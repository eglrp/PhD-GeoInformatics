/******************************************************************************
 * $Id: CrossCalibration.cpp $
 *
 * Project:  Cross calibration
 * Purpose:  Radiometric calibration of VHR/aerial imagery to well known 
 *           BRDF calibrated reference like MODIS, at lower res.
 * Author:   Dugal Harris <dugalh@gmail.com>
 ******************************************************************************
 * Copyright (c) 2014 Dugal Harris
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

 /******************************************************************************
 * Some Notes:
 * On array and disk storage: 
 * I is for row, J is for col - I have tried to standardise on this
 * OpenCV and Buf3d API generally take indices in (row, col) / (i, j) order but...
 * GDAL (eg RasterIO) takes indices in (x, y) / (col, row) order 
 * RasterIO returns data in col, row, channel order which can be stored in a 3dBuf, so...
 * 3dBuf storage is also in col, row, channel order i.e. along the col dimension first, then ... etc i.e. 3dBuf can point to RasterIO data w/o issue
 * OpenCV array storage is in pixel(channel), col, row order (for 2D mats - I think 3D mats will be col, row, channel)
 * The storage of RasterIO, OpenCV, 3dBuf is called "row major order" as the data is stored first down the col dimension, this is a row hence "row major"
 * Matlab uses column major order 
 * So note that the OpenCV API generally takes indices in the opposite order to storage
 * GTiff storage on disk depends on interleave and tiling settings but I suspect is also col, row order.  Pixel interleave, which seems to be GTiff default,
 * means channel dimension first.  Data is stored tile-by-tile so that individual tiles are contiguous and can be read from disk as such.  When there is no
 * tiling specified, the tile can still be thought of to be the size of 1 row (and gdal_info reports the block/tile size as such).  
 * A row tile may not be ideal for display or compression purposes where a 2d ROI is what is actually been accessed and operated on.  

 ******************************************************************************/

 //

#include "targetver.h"

#include <tchar.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <omp.h>

#include "CrossCalibration.h"

using namespace std;

#define XCALIB_DEBUG 0
#define MAX_PATH 1024
#define SEAMLINE_COVERAGE_FIX 1					// excludes partially covered ref pixels
#define SEAMLINE_EXTRAP_FIX 0					// erodes away one boundary pixel to exclude extrapolated params, SEAMLINE_COVERAGE_FIX must be 1
#define DO_COMPRESS_OUTPUT 1
#define BIGTIFF 1
#define COVERAGE_PORTION 0.95
#define SLIDE_WIN_CENTER 0

int defaultWinSize[2] = {1, 1};
int gdalwarp(int argc, char ** argv);
int gdal_translate(int argc, char ** argv);

//-------------------------------------------------------------------------------------------------
//Parses argString and calls gdalwarp backend
//-------------------------------------------------------------------------------------------------
//TODO: I don't think we actually need gdalwarp - gdaltranslate can resample.  gdalwarp is for reprojecting / mosaicing.  
//TODO: Also, we may not need gdal_translate as it essentially seems to be a call to CopyDataset (or something similar...)
int GdalWarpWrapper(const string& argString)
{
	int numFields = (int)count(argString.begin(), argString.end(), '~');
	char** argv = new char*[numFields];
	int res = 0;
	stringstream ss(argString);
	string item;
	int i = 0;
	while (getline(ss, item, '~'))
	{
		argv[i] = new char[MAX_PATH];
		strcpy_s(argv[i++], MAX_PATH, item.c_str());
	}
	
	res = gdalwarp(i, argv); //GDAL utility function
	
	for (int j = 0; j < numFields; j++)
	{
		if (argv[j] != NULL)
			delete[] argv[j];
		argv[j] = NULL;
	}
	delete[] argv;
	return res;
}

//-------------------------------------------------------------------------------------------------
//Parses argString and calls gdal_translate backend
//-------------------------------------------------------------------------------------------------
int GdalTranslateWrapper(const string& argString)
{
	int numFields = (int)count(argString.begin(), argString.end(), '~');
	char** argv = new char*[numFields];
	int res = 0;
	stringstream ss(argString);
	string item;
	int i = 0;
	while (getline(ss, item, '~'))
	{
		argv[i] = new char[MAX_PATH];
		strcpy_s(argv[i++], MAX_PATH, item.c_str());
	}

	res = gdal_translate(i, argv); //GDAL utility function

	for (int j = 0; j < numFields; j++)
	{
		if (argv[j] != NULL)
			delete[] argv[j];
		argv[j] = NULL;
	}
	delete[] argv;
	return res;
}

//-------------------------------------------------------------------------------------------------
//Error checking (projections, resolutions, spatial coverage)
//-------------------------------------------------------------------------------------------------
void ErrorCheckDatasets(GDALDataset* srcDataSet, GDALDataset* refDataSet)
{
	//check band count
	if (srcDataSet->GetRasterCount() != refDataSet->GetRasterCount())
		throw string("source and reference images don't have the same number of bands");

	//compare co-ord systems
	OGRSpatialReference srcRef, refRef;
	const char* srcProjectionRef = srcDataSet->GetProjectionRef();
	const char* refProjectionRef = refDataSet->GetProjectionRef();

	OGRErr err = srcRef.importFromWkt((char**)&srcProjectionRef);
	err = refRef.importFromWkt((char**)&refProjectionRef);
	if (srcRef.IsSame(&refRef) == 0)
		throw string("source and reference co-ord systems are not the same");

	//check srcDataSet lies inside refDataSet
	double refGeoTransform[6], srcGeoTransform[6];
	refDataSet->GetGeoTransform(refGeoTransform);
	srcDataSet->GetGeoTransform(srcGeoTransform);

	double srcX[2], srcY[2], refX[2], refY[2];
	GDALApplyGeoTransform(srcGeoTransform, 0, 0, &srcX[0], &srcY[0]); //upper left
	GDALApplyGeoTransform(srcGeoTransform, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterYSize(), &srcX[1], &srcY[1]); // bottom right
	GDALApplyGeoTransform(refGeoTransform, 0, 0, &refX[0], &refY[0]);
	GDALApplyGeoTransform(refGeoTransform, refDataSet->GetRasterXSize(), refDataSet->GetRasterYSize(), &refX[1], &refY[1]);

	OGRCoordinateTransformation* srcTForm = OGRCreateCoordinateTransformation(&srcRef, srcRef.CloneGeogCS());
	OGRCoordinateTransformation* refTform = OGRCreateCoordinateTransformation(&refRef, refRef.CloneGeogCS());

	srcTForm->Transform(2, srcX, srcY); //to lat lon (not necessary)
	refTform->Transform(2, refX, refY);

	delete srcTForm;
	delete refTform;

	double refXMin = *std::min_element(refX, refX+2), refYMin = *std::min_element(refY, refY+2);
	double refXMax = *std::max_element(refX, refX+2), refYMax = *std::max_element(refY, refY+2);

	if (!((srcX[0] >= refXMin) && (srcX[1] >= refXMin) && (srcX[0] <= refXMax) && (srcX[1] <= refXMax) &&
		(srcY[0] >= refYMin) && (srcY[1] >= refYMin) && (srcY[0] <= refYMax) && (srcY[1] <= refYMax)))
		throw string("source image does not lie inside the reference image"); 

	//check resolutions
	if (!(abs(srcGeoTransform[1]) <= abs(refGeoTransform[1]) && abs(srcGeoTransform[5]) <= abs(refGeoTransform[5])))
		throw string("source image resolution should be finer than reference image resolution"); 
}

//-------------------------------------------------------------------------------------------------
//Cross calibrates srcFileName to refFileName
//-------------------------------------------------------------------------------------------------
int CrossCalib(const string& refFileName, const string& srcFileName_, const int* winSize = defaultWinSize, const int modelForm = ModelForms::GAIN_ONLY)
{
	string srcFileName = srcFileName_;
	char* argv[20] = {};
	char gdalString[MAX_PATH];
	stringstream ss(gdalString);
	string item;
	int i = 0;
	int res = 0;
	GDALDataset* refDataSet = NULL;
	//int winSize[2] = {1, 1}; //TO DO: make configurable
	
#if 0 //XCALIB_DEBUG 
	//--------------------------------------------------------------------------------------------------------------
	//Start by downsampling the src to 2.5m pixels to speed up init testing
	//-------------------------------------------------------------------------------------------------

	string srcInitFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_INIT.tif";
	
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~\"0\"~-dstnodata~\"0\"~-wm~1024~-r~bilinear~-tap~-tr~%1.1f~%1.1f~%s~%s~", 2.5, 2.5, 
		srcFileName.c_str(), srcInitFileName.c_str());
	res = GdalWarpWrapper(gdalString);

	srcFileName = srcInitFileName;
#endif

	//--------------------------------------------------------------------------------------------------------------
	//Downsample the source file to the ref file resolution & grid
	//-------------------------------------------------------------------------------------------------

	double refGeoTransform[6], srcGeoTransform[6];
	refDataSet = (GDALDataset*)GDALOpen(refFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (refDataSet == NULL)
		throw string("Could not open: " + refFileName);
	refDataSet->GetGeoTransform(refGeoTransform);

	GDALDataset* srcDataSet = (GDALDataset*)GDALOpen(srcFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcDataSet == NULL)
		throw string("Could not open: " + srcFileName);

	ErrorCheckDatasets(srcDataSet, refDataSet);


	srcDataSet->GetGeoTransform(srcGeoTransform);
	int srcXSize = srcDataSet->GetRasterXSize();
	int srcYSize = srcDataSet->GetRasterYSize();
	GDALClose(srcDataSet);

	std::cout << "------------------------------------------------------------" << endl;
	std::cout << "Processing: " << srcFileName << endl;
	cout << "Reference file: " << srcFileName << endl;
	cout << "Window size: " << winSize[0] << ", " << winSize[1] << endl;
	cout << "Model: " << ModelFormsToString(modelForm) << endl;
	std::cout << "------------------------------------------------------------" << endl;

#if SEAMLINE_COVERAGE_FIX   // TO DO: somehow change seamline fix to improve processing Eg only downsample/upsample fully covered ref areas to reduce processing by figuring 
					// out the full coverage mask, then update nodata in source im with this.  OR: avoid upsampling the DS mask - just make one rather than being lazy 
	std::cout << "Extracting mask from: " << srcFileName << endl << endl;
#if XCALIB_DEBUG
	string srcMaskFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_MASK.tif"; //make separate files for each source file
#else
	string srcMaskFileName = string(CPLGetCurrentDir()) + "\\SourceMask.tif";
#endif

	//gdal_translate - b mask - a_nodata none - ot Float32 3321b_3172_12_0415_rgbn.tif 3321b_3172_12_0415_nodata.tif
	sprintf_s(gdalString, MAX_PATH, "gdal_translate~-b~mask~-a_nodata~none~-ot~Float32~%s~%s~", srcFileName.c_str(), srcMaskFileName.c_str());
	res = GdalTranslateWrapper(gdalString);
	if (res != 0)
		throw string("Could not extract mask " + srcFileName);

	std::cout << "Downsampling mask: " << srcMaskFileName << endl << endl;

#if XCALIB_DEBUG
	string srcMaskDsFileName = srcMaskFileName.substr(0, srcMaskFileName.length() - 4) + "_DS.tif"; //make separate files for each source file
#else
	string srcMaskDsFileName = string(CPLGetCurrentDir()) + "\\SourceMaskDs.tif"; //reuse the same file
#endif

	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~none~-dstnodata~none~-r~average~-tap~-tr~%d~%d~%s~%s~", (int)abs(refGeoTransform[1]), (int)abs(refGeoTransform[5]),
		srcMaskFileName.c_str(), srcMaskDsFileName.c_str());
	res = GdalWarpWrapper(gdalString);

	if (res != 0)
		throw string("Could not warp " + srcMaskFileName);

#endif //SEAMLINE_COVERAGE_FIX

	std::cout << "Downsampling: " << srcFileName << endl << endl;

#if XCALIB_DEBUG
	string srcDsFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_DS.tif"; //make separate files for each source file
#else
	//string srcDsFileName = string(CPLGetCurrentDir()) + "\\SourceDs.tif"; //reuse the same file
	string srcDsFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_DS.tif"; //make separate files for each source file
#endif
	//

	//TODO if the raster already has average(/cubicspline? - to test) overviews, it seems this downsampling can go much faster
	/*if (modelForm == ModelForms::OFFSET_ONLY)
		sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-overwrite~-srcnodata~\"0\"~-dstnodata~\"0\"~-r~average~-tap~-tr~%d~%d~%s~%s~", (int)abs(refGeoTransform[1]), (int)abs(refGeoTransform[5]),
			srcFileName.c_str(), srcDsFileName.c_str());
	else*/
		sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-overwrite~-srcnodata~\"0\"~-dstnodata~\"0\"~-r~cubicspline~-tap~-tr~%d~%d~%s~%s~", (int)abs(refGeoTransform[1]), (int)abs(refGeoTransform[5]), 
			srcFileName.c_str(), srcDsFileName.c_str());
	res = GdalWarpWrapper(gdalString);

	if (res != 0)
		throw string("Could not warp " + srcFileName);

	//--------------------------------------------------------------------------------------------------------------
	//Find the ROI in the ref file that corresponds to the source file
	//-------------------------------------------------------------------------------------------------

	GDALDataset* srcDsDataSet = (GDALDataset*)GDALOpen(srcDsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcDsDataSet == NULL)
		throw string("Could not open: " + srcDsFileName);
	double srcDsGeoTransform[6];
	double refInvGeoTransform[6];

	srcDsDataSet->GetGeoTransform(srcDsGeoTransform);
	//refDataSet->GetGeoTransform(refGeoTransform);
	GDALInvGeoTransform(refGeoTransform, refInvGeoTransform);

	double ulX, ulY, brX, brY; //geog co-ords of source
	double refUlI, refUlJ, refBrI, refBrJ; //pixel co-ords in ref

	GDALApplyGeoTransform(srcDsGeoTransform, 0, 0, &ulX, &ulY);
	GDALApplyGeoTransform(srcDsGeoTransform, srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(), &brX, &brY);

	GDALApplyGeoTransform(refInvGeoTransform, ulX, ulY, &refUlJ, &refUlI);
	GDALApplyGeoTransform(refInvGeoTransform, brX, brY, &refBrJ, &refBrI);

	//--------------------------------------------------------------------------------------------------------------
	//Find the calibration params
	//-------------------------------------------------------------------------------------------------

	int nParamBands = srcDsDataSet->GetRasterCount();
	if (modelForm == ModelForms::GAIN_AND_OFFSET)  // store both gain and offset in same raster (gain first, then offset)
		nParamBands *= 2;
	//TODO remove all Buf3d and use only cv::Mat - will need to be a 3D Mat though which is not working in this version
	Buf3d<float> srcDsData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());
	Buf3d<float> refSubData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());
	Buf3d<float> paramDsData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), nParamBands);

#if XCALIB_DEBUG
	string paramDsFileName = srcDsFileName.substr(0, srcDsFileName.length() - 4) + "_PARAMS.tif"; //make separate files for each source file
#else
	//string paramDsFileName = string(CPLGetCurrentDir()) + "\\ParamDs.tif"; //reuse the same file
	string paramDsFileName = srcDsFileName.substr(0, srcDsFileName.length() - 4) + "_PARAMS.tif"; //HACK for now - these files are tiny
#endif

	std::cout << "Calculating calibration params (" << paramDsFileName << ")" << endl << endl;

	GDALDataset* paramDsDataSet = srcDsDataSet->GetDriver()->Create(paramDsFileName.c_str(), srcDsDataSet->GetRasterXSize(), 
		srcDsDataSet->GetRasterYSize(), nParamBands, GDALDataType::GDT_Float32, NULL);

	if (paramDsDataSet == NULL)
		throw string("Could not create: " + paramDsFileName);

	paramDsDataSet->SetGeoTransform(srcDsGeoTransform);
	paramDsDataSet->SetProjection(srcDsDataSet->GetProjectionRef());	

	CPLErr err = srcDsDataSet->RasterIO(GDALRWFlag::GF_Read, 0, 0, srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(),
		srcDsData.Buf(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float32, 
		srcDsDataSet->GetRasterCount(), NULL, 0, 0, 0);
	
	if (err != CPLErr::CE_None)
		throw string("srcDsDataSet->RasterIO failure");

	err = refDataSet->RasterIO(GDALRWFlag::GF_Read, (int)refUlJ, (int)refUlI, srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(),
		refSubData.Buf(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float32, 
		srcDsDataSet->GetRasterCount(), NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("refDataSet->RasterIO failure");


#if SEAMLINE_COVERAGE_FIX
	//read the ds mask
	GDALDataset* srcMaskDsDataSet = (GDALDataset*)GDALOpen(srcMaskDsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcMaskDsDataSet == NULL)
		throw string("Could not open: " + srcMaskDsFileName);
	Buf3d<float> srcMaskDsData(srcMaskDsDataSet->GetRasterYSize(), srcMaskDsDataSet->GetRasterXSize(), srcMaskDsDataSet->GetRasterCount());

	err = srcMaskDsDataSet->RasterIO(GDALRWFlag::GF_Read, 0, 0, srcMaskDsDataSet->GetRasterXSize(), srcMaskDsDataSet->GetRasterYSize(),
		srcMaskDsData.Buf(), srcMaskDsDataSet->GetRasterXSize(), srcMaskDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float32,
		srcMaskDsDataSet->GetRasterCount(), NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("srcMaskDsData->RasterIO failure");
#if SEAMLINE_EXTRAP_FIX  //TODO change the coverage mask to a vector cutline that can be shrunk inwards for the extrap fix.  then we don't have to create & store us and ds mask images which is clumsy, space and time consuming
	//make eroded mask for application to upsampled images
	//string winTitle = "mask";
	cv::Mat cvMask(srcMaskDsDataSet->GetRasterYSize(), srcMaskDsDataSet->GetRasterXSize(), CV_32FC1, (void*)srcMaskDsData.Buf());
	//cv::namedWindow(winTitle, cv::WINDOW_NORMAL | CV_WINDOW_KEEPRATIO | CV_GUI_EXPANDED);
	//cv::imshow(winTitle, cvMask / 255);
	//cv::waitKey(0);
	cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));
	cv::Mat cvMaskErode;
	//string winTitle2 = "maskErode";
	//cv::namedWindow(winTitle2, cv::WINDOW_NORMAL | CV_WINDOW_KEEPRATIO | CV_GUI_EXPANDED);
	cv::erode(cvMask >= 255, cvMaskErode, kernel);
	//cv::imshow(winTitle2, cvMaskErode);
	//cv::waitKey(0);
	//cv::destroyAllWindows();

#if XCALIB_DEBUG
	string erodeMaskDsFileName = srcMaskDsFileName.substr(0, srcMaskDsFileName.length() - 4) + "_ERODE.tif"; //make separate files for each source file
#else
	string erodeMaskDsFileName = string(CPLGetCurrentDir()) + "\\ErodeMaskDs.tif"; //reuse the same file
#endif

	GDALDataset* erodeMaskDsDataSet = srcDsDataSet->GetDriver()->Create(erodeMaskDsFileName.c_str(), srcMaskDsDataSet->GetRasterXSize(),
		srcMaskDsDataSet->GetRasterYSize(), 1, GDALDataType::GDT_Byte, NULL);

	if (erodeMaskDsDataSet == NULL)
		throw string("Could not create: " + paramDsFileName);

	erodeMaskDsDataSet->SetGeoTransform(srcDsGeoTransform);
	erodeMaskDsDataSet->SetProjection(srcDsDataSet->GetProjectionRef());

	err = erodeMaskDsDataSet->RasterIO(GDALRWFlag::GF_Write, 0, 0, erodeMaskDsDataSet->GetRasterXSize(), erodeMaskDsDataSet->GetRasterYSize(),
		cvMaskErode.data, erodeMaskDsDataSet->GetRasterXSize(), erodeMaskDsDataSet->GetRasterYSize(), GDALDataType::GDT_Byte,
		1, NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("erodeMaskDsDataSet->RasterIO failure");

	erodeMaskDsDataSet->FlushCache();
	GDALClose(erodeMaskDsDataSet);

#else //#if SEAMLINE_EXTRAP_FIX
	string erodeMaskDsFileName = srcMaskDsFileName;  // don't do erosion
#endif  //#if SEAMLINE_EXTRAP_FIX
	GDALClose(srcMaskDsDataSet);
	cv::Mat& srcMaskDsMat = srcMaskDsData.ToMat(0);
#else 
	cv::Mat_<unsigned char> srcMaskDsMat(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), (unsigned char)255);
#endif  //#if SEAMLINE_COVERAGE_FIX

	cv::Range winRanges[2];

	std::vector<cv::Mat_<float>>& refSubBands = refSubData.ToMatVec();
	std::vector<cv::Mat_<float>>& srcDsBands = srcDsData.ToMatVec();
	//cv::Mat onesVec;  //(int rows, int cols, int type, const Scalar& s)
	//cv::Mat srcDsConcatMat;

	int nRows = srcDsDataSet->GetRasterYSize(), nCols = srcDsDataSet->GetRasterXSize();
	
	std::vector<float> imageGain(srcDsDataSet->GetRasterCount(), 1.);
	std::vector<float> imageOffset(srcDsDataSet->GetRasterCount(), 0.);
	cout << "Image (gain, offset): ";
	for (int k = 0; k < srcDsDataSet->GetRasterCount(); k++)
	{
#if FALSE   // this approach has a problem
		// the offset src can go close to zero which means the ref/(src+offset) gains can grow large and change +ve to -ve
		// this makes a mess of things when upsampling and these wildly varying gains must be interpolated between.
		// also, I don't think it is physically valid.  It should be more like the min(src) should be offset to be similar to
		// min(ref) which the alternative below seems to achieve 
		cv::Scalar srcMean, srcStd, refSubMean, refSubStd;
		cv::meanStdDev(srcDsBands[k], srcMean, srcStd, srcMaskDsMat >= COVERAGE_PORTION*255.0);
		cv::meanStdDev(refSubBands[k], refSubMean, refSubStd, srcMaskDsMat >= COVERAGE_PORTION*255.0);
		if (modelForm == ModelForms::IMAGE_GAIN_AND_OFFSET)  // find per band gains for IMAGE_GAIN_AND_OFFSET model
		{
			//cv::Scalar srcMean, srcStd, refSubMean, refSubStd;
			//cv::meanStdDev(srcDsBands[k], srcMean, srcStd, srcMaskDsMat >= 0.95*255.0);
			//cv::meanStdDev(refSubBands[k], refSubMean, refSubStd, srcMaskDsMat >= 0.95*255.0);
			//bandScaleForOffsetModel[k] = refSubStd.val[0] / srcStd.val[0];
			imageGain[k] = refSubStd.val[0] / srcStd.val[0];
		}
		else if (modelForm == ModelForms::GAIN_AND_IMAGE_OFFSET)  // find per band offsets for GAIN_AND_IMAGE_OFFSET model
		{
			cv::meanStdDev(srcDsBands[k] * refSubStd.val[0] / srcStd.val[0], srcMean, srcStd, srcMaskDsMat >= COVERAGE_PORTION*255.0);
			imageOffset[k] = refSubMean.val[0] - srcMean.val[0];
		}
#else
		cv::Mat coverageMask = srcMaskDsMat >= COVERAGE_PORTION*255.0;
		cv::Mat onesVec(nCols*nRows, 1, CV_32FC1, cv::Scalar(1.));  //(int rows, int cols, int type, const Scalar& s)
		cv::Mat srcDsConcatMat(nCols*nRows, 2, CV_32FC1, cv::Scalar(1.));  //(int rows, int cols, int type, const Scalar& s)
																	//srcDsConcatMat = cv::Mat(n*m, 2, CV_32FC1, cv::Scalar(0.));
		cv::Mat lsSoln;
		cv::hconcat(srcDsBands[k].clone().setTo(0, ~coverageMask).reshape(1, nCols*nRows), onesVec, srcDsConcatMat);
		cv::solve(srcDsConcatMat, refSubBands[k].clone().setTo(0, ~coverageMask).reshape(1, nCols*nRows), lsSoln, cv::DECOMP_SVD);
		//cv::solve(srcDsConcatMat, refSubBands[k].clone().setTo(0, ~coverageMask).reshape(1, nCols*nRows), lsSoln, cv::DECOMP_SVD);
		if (modelForm == ModelForms::IMAGE_GAIN_AND_OFFSET)  // find per band gains for IMAGE_GAIN_AND_OFFSET model
			imageGain[k] = 1./lsSoln.at<float>(0);
		else if (modelForm == ModelForms::GAIN_AND_IMAGE_OFFSET)  // find per band offsets for GAIN_AND_IMAGE_OFFSET model
			imageOffset[k] = lsSoln.at<float>(1)/lsSoln.at<float>(0);
#endif
		cout << "(" << imageGain[k] << ", " << imageOffset[k] << "), ";
	}
	cout << endl;

	int winCenterOffset[2] = {(winSize[0] - 1)/2, (winSize[1] - 1)/2};
	//TODO: As an alternative to the GAIN_AND_IMAGE_OFFSET model above, we could have different window sizes for gain and offset 
	//params, use the setTo(0) & onesVec in conjunction with Mx+C model to implement this.  This way we could have a larger offset window
	//and smaller gain window.  Or perhaps just set offset window to full image by default.

	for (int k = 0; k < srcDsDataSet->GetRasterCount(); k++)
	{
		err = paramDsDataSet->GetRasterBand(k + 1)->SetNoDataValue(0);
		cv::Mat& srcDsBand = srcDsBands[k];
		cv::Mat& refSubBand = refSubBands[k];
		//TODO (add option to) make i,j the window center not UL cnr - that way we will have more valid xcalib data
		//TODO make notes on SLIDE_WIN_CENTER behaviour, disk storage ordering, compression bottlenecks
#if SLIDE_WIN_CENTER	
		for (int j = 0; j < nCols; j++)
		{
			winRanges[1] = cv::Range(max<int>(j - winCenterOffset[1], 0), min<int>(nCols, j + winSize[1] - winCenterOffset[1]));
			for (int i = 0; i < nRows; i++)
			{
				winRanges[0] = cv::Range(max<int>(i - winCenterOffset[0], 0), min<int>(nRows, i + winSize[0] - winCenterOffset[0]));
				// extract sliding window ROI - (i, j) is center
				const cv::Mat srcDsWin = srcDsBand(winRanges);
				const cv::Mat refSubWin = refSubBand(winRanges);
				const cv::Mat srcDsMaskDsWin = srcMaskDsMat(winRanges);
				const cv::Mat::MSize coveredWinSize = srcDsWin.size;

				cv::Mat coverageMask = srcDsMaskDsWin >= (COVERAGE_PORTION*255.0);
				// convert sliding win data to col vectors etc suitable for LS
				// find LS param estimates
				//int valid = (cv::compare( refSubWin == 0);
				if (modelForm == ModelForms::GAIN_ONLY || modelForm == ModelForms::GAIN_AND_IMAGE_OFFSET)
				{
					float gain = cv::mean(refSubWin / (srcDsWin + imageOffset[k]), coverageMask).val[0];
					paramDsData(i, j, k) = gain;
					if (abs(gain) > 10)
					{
						cout << "coverageMask" << coverageMask << endl;
						cout << "refSubWin" << refSubWin << endl;
						cout << "srcDsWin" << srcDsWin << endl;
					}
				}
				else if (modelForm == ModelForms::GAIN_AND_OFFSET)
				{
					cv::Mat lsSoln;
					float coverage = cv::sum(coverageMask).val[0];


					//#if XCALIB_DEBUG
					//#endif
					if (coverage >= 255. * 2)  //check we have at least 2 pieces of valid data for our 2 params
					{
						cv::Mat onesVec = cv::Mat(coveredWinSize[0] * coveredWinSize[1], 1, CV_32FC1, cv::Scalar(1.));  //(int rows, int cols, int type, const Scalar& s)
						cv::Mat srcDsConcatMat = cv::Mat(coveredWinSize[0] * coveredWinSize[1], 2, CV_32FC1, cv::Scalar(0.));
						cv::Mat srcVec = srcDsWin.clone().setTo(0, ~coverageMask).reshape(1, coveredWinSize[0] * coveredWinSize[1]);
						cv::Mat refVec = refSubWin.clone().setTo(0, ~coverageMask).reshape(1, coveredWinSize[0] * coveredWinSize[1]);
						cv::Mat coveredOnesVec = onesVec.setTo(0, ~coverageMask.reshape(1, coveredWinSize[0] * coveredWinSize[1]));
						/*cout << "(i, j)" << i << ", " << j << endl;
						cout << "coverageMask:" << coverageMask << endl;
						cout << "onesVec:" << onesVec << endl;
						cout << "coveredOnesVec:" << coveredOnesVec << endl;
						cout << "srcDsConcatMat:" << srcDsConcatMat << endl;
						cout << "srcDsWin:" << srcDsWin << endl;
						cout << "srcVec:" << srcVec << endl;
						cout << "refSubWin:" << refSubWin << endl;
						cout << "refVec:" << refVec << endl;*/
						//if (coverage < winSize[0] * winSize[0] * 255)
						{
							cv::hconcat(srcVec, coveredOnesVec, srcDsConcatMat);
							cv::solve(srcDsConcatMat, refVec, lsSoln, cv::DECOMP_QR);
						}
						//cout << "lsSoln:" << lsSoln << endl;
						//else
						//{
						//	cv::Mat srcVec = srcDsWin.clone().reshape(1, srcDsWin.rows*srcDsWin.cols);
						//	cv::hconcat(srcVec, onesVec, srcDsConcatMat);
						//	cv::solve(srcDsConcatMat, refSubWin.clone().reshape(1, refSubWin.rows*refSubWin.cols), lsSoln, cv::DECOMP_SVD);
						//}

					}
					else
						lsSoln = cv::Mat(2, 1, CV_32F, cv::Scalar::all(0));  // if we dont have enough coverage, set params to zero 

					paramDsData(i, j, k) = lsSoln.at<float>(0, 0);
					paramDsData(i, j, k + srcDsDataSet->GetRasterCount()) = lsSoln.at<float>(1, 0);
				}
				else if (modelForm == ModelForms::OFFSET_ONLY || modelForm == ModelForms::IMAGE_GAIN_AND_OFFSET)
				{
					float offset = 0.;
					if (cv::sum(coverageMask).val[0] > 0.)
						offset = cv::mean(refSubWin - srcDsWin * imageGain[k], coverageMask).val[0];
					paramDsData(i, j, k) = offset;
				}

				// place param estimates at center of sliding win
			}
		}
#else  //SLIDE_WIN_ULCNR

		cv::Mat onesVec = cv::Mat(winSize[0] * winSize[1], 1, CV_32FC1, cv::Scalar(1.));  //(int rows, int cols, int type, const Scalar& s)
		cv::Mat srcDsConcatMat = cv::Mat(winSize[0] * winSize[1], 2, CV_32FC1, cv::Scalar(0.));

		for (int j = 0; j < nCols - winSize[1] + 1; j++)
		{
			//int colIdx[2] = {max<int>(j-(winSize[1]-1)/2, 0), min<int>(m-1, j+(winSize[1]-1)/2)};
			winRanges[1] = cv::Range(j, j + winSize[1]);
			for (int i = 0; i < nRows - winSize[0] + 1; i++)
			{
				//int rowIdx[2] = {max<int>(i-(winSize[0]-1)/2, 0), min<int>(n-1, i+(winSize[0]-1)/2)};
				winRanges[0] = cv::Range(i, i + winSize[0]);
				// extract sliding window ROI - (i, j) is UL cnr
				const cv::Mat srcDsWin = srcDsBand(winRanges);
				const cv::Mat refSubWin = refSubBand(winRanges);
				const cv::Mat srcDsMaskDsWin = srcMaskDsMat(winRanges);
				cv::Mat coverageMask = srcDsMaskDsWin >= COVERAGE_PORTION * 255.0;

				// convert sliding win data to col vectors etc suitable for LS
				// find LS param estimates
				//int valid = (cv::compare( refSubWin == 0);
				if (modelForm == ModelForms::GAIN_ONLY || modelForm == ModelForms::GAIN_AND_IMAGE_OFFSET)
				{
					float gain = cv::mean(refSubWin / (srcDsWin + imageOffset[k]), coverageMask).val[0];
					paramDsData(i + winCenterOffset[0], j + winCenterOffset[1], k) = gain;
				}
				else if (modelForm == ModelForms::GAIN_AND_OFFSET)
				{
					cv::Mat lsSoln;
					float coverage = cv::sum(coverageMask).val[0];
					if (coverage >= 255. * 2)  //check we have at least 2 pieces of valid data for our 2 params
					{
						//cv::Mat onesVec = cv::Mat(winSize[0] * winSize[1], 1, CV_32FC1, cv::Scalar(1.));  //(int rows, int cols, int type, const Scalar& s)
						//cv::Mat srcDsConcatMat = cv::Mat(winSize[0] * winSize[1], 2, CV_32FC1, cv::Scalar(0.));
						onesVec.setTo(1.);   //NB
						cv::Mat srcVec = srcDsWin.clone().setTo(0, ~coverageMask).reshape(1, winSize[0] * winSize[1]);
						cv::Mat refVec = refSubWin.clone().setTo(0, ~coverageMask).reshape(1, winSize[0] * winSize[1]);
						cv::Mat coveredOnesVec = onesVec.setTo(0, ~coverageMask.reshape(1, winSize[0] * winSize[1]));

						/*cout << "(i, j)" << i << ", " << j << endl;
						cout << "coverageMask:" << coverageMask << endl;
						cout << "onesVec:" << onesVec << endl;
						cout << "coveredOnesVec:" << coveredOnesVec << endl;
						cout << "srcDsConcatMat:" << srcDsConcatMat << endl;
						cout << "srcDsWin:" << srcDsWin << endl;
						cout << "srcVec:" << srcVec << endl;
						cout << "refSubWin:" << refSubWin << endl;
						cout << "refVec:" << refVec << endl;*/

						cv::hconcat(srcVec, coveredOnesVec, srcDsConcatMat);
						cv::solve(srcDsConcatMat, refVec, lsSoln, cv::DECOMP_SVD);

						//cout << "lsSoln:" << lsSoln << endl;


						//if (coverage < winSize[0] * winSize[0] * 255)
						//{
						//	cv::hconcat(srcDsWin.clone().setTo(0, ~coverageMask).reshape(1, winSize[0] * winSize[0]),
						//		onesVec.setTo(0, ~coverageMask.reshape(1, winSize[0] * winSize[0])), srcDsConcatMat);
						//	cv::solve(srcDsConcatMat, refSubWin.clone().setTo(0, ~coverageMask).reshape(1, winSize[0] * winSize[0]), lsSoln, cv::DECOMP_SVD);
						//}
						//else
						//{
						//	cv::hconcat(srcDsWin.clone().reshape(1, winSize[0] * winSize[0]), onesVec, srcDsConcatMat);
						//	cv::solve(srcDsConcatMat, refSubWin.clone().reshape(1, winSize[0] * winSize[0]), lsSoln, cv::DECOMP_SVD);
						//}
					}
					else
						lsSoln = cv::Mat(2, 1, CV_32F, cv::Scalar::all(0));  // if we dont have enough coverage, set params to zero 

					paramDsData(i + winCenterOffset[0], j + winCenterOffset[1], k) = lsSoln.at<float>(0, 0);
					paramDsData(i + winCenterOffset[0], j + winCenterOffset[1], k + srcDsDataSet->GetRasterCount()) = lsSoln.at<float>(1, 0);
				}
				else if (modelForm == ModelForms::OFFSET_ONLY || modelForm == ModelForms::IMAGE_GAIN_AND_OFFSET)
				{
					//double offset = cv::mean(refSubWin - srcDsWin * bandScaleForOffsetModel[k], coverageMask).val[0];
					//paramDsData(i + (winSize[0] - 1) / 2, j + (winSize[1] - 1) / 2, k) = offset;
					float offset = cv::mean(refSubWin - srcDsWin * imageGain[k], coverageMask).val[0];
					paramDsData(i + winCenterOffset[0], j + winCenterOffset[1], k) = offset;
				}
			}
		}
#endif
	}

	//write out params to file
	err = paramDsDataSet->RasterIO(GDALRWFlag::GF_Write, 0, 0, paramDsDataSet->GetRasterXSize(), paramDsDataSet->GetRasterYSize(),
		paramDsData.Buf(), paramDsDataSet->GetRasterXSize(), paramDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float32, 
		nParamBands, NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("paramDsDataSet->RasterIO failure");

	paramDsDataSet->FlushCache();
	GDALClose(paramDsDataSet);

	GDALClose(srcDsDataSet);
	GDALClose(refDataSet);
	//--------------------------------------------------------------------------------------------------------------
	//Upsample param to source resolution
	//-------------------------------------------------------------------------------------------------

#if XCALIB_DEBUG
	string paramUsFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_US_PARAM.tif";
#else
	string paramUsFileName = string(CPLGetCurrentDir()) + "\\ParamUs.tif";
#endif
	std::cout << "Interpolating calibration params (" << paramUsFileName << ")" << endl << endl;

	//TODO: check if upsampling faster with square tile size i.e. tiles are not rows
	//TODO: use new gdal warp code as it may be faster or better threading
#if BIGTIFF
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-co~BIGTIFF=YES~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), paramDsFileName.c_str(), paramUsFileName.c_str());
#else
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), paramDsFileName.c_str(), paramUsFileName.c_str());
#endif

	res = GdalWarpWrapper(gdalString);
	if (res != 0)
		throw string("Could not warp " + paramDsFileName);


	//-------------------------------------------------------------------------------------------------
	//Apply umsampled params to source image
	//-------------------------------------------------------------------------------------------------

	string calibFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_XCALIB.tif";
	std::cout << "Applying params (" << calibFileName << ")" << endl << endl;

	srcDataSet = (GDALDataset*)GDALOpen(srcFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcDataSet == NULL)
		throw string("Could not open: " + srcFileName);

	GDALDataset* paramDataSet = (GDALDataset*)GDALOpen(paramUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (paramDataSet == NULL)
		throw string("Could not open: " + paramUsFileName);
	char **papszOptions = NULL;
#if DO_COMPRESS_OUTPUT
	// make block size same as src?  to speed writing?
	// the bottleneck in the xcalib creation below is mostly in the deflate compression, then in disk access - threading helps a lot
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");  //LZW faster than DEFLATE
	papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "2");
#endif
	//it seems tiff creates a default block size of a row though for large tiffs (if TILED=NO or not specified I think).
	//the above options seem to force this default behaviour which is what we want for now 
	//we may be better of copying the original dataset though as in gdal_translate for both param and xcalib files. 
	//this way, they all have same block size which we need for efficiency below
	//even if the tiff is not specifically tiled, it is "tiled" by row by default
	int xBlockSize = 0, yBlockSize = 0;
	srcDataSet->GetRasterBand(1)->GetBlockSize(&xBlockSize, &yBlockSize);
	papszOptions = CSLSetNameValue(papszOptions, "BLOCKXSIZE", std::to_string(xBlockSize).c_str());
	papszOptions = CSLSetNameValue(papszOptions, "BLOCKYSIZE", std::to_string(yBlockSize).c_str());
	if (yBlockSize > 1)
		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");  //only if src dataset appears to be tiled
	papszOptions = CSLSetNameValue(papszOptions, "NUM_THREADS", "ALL_CPUS");
#if BIGTIFF
	papszOptions = CSLSetNameValue(papszOptions, "BIGTIFF", "YES");
#endif

	GDALDataset* calibDataSet = srcDataSet->GetDriver()->Create(calibFileName.c_str(), srcDataSet->GetRasterXSize(),
		srcDataSet->GetRasterYSize(), 4, GDALDataType::GDT_Int16, papszOptions);  //int16 to allow -ve values for offset models gone wrong
	if (calibDataSet == NULL)
		throw string("Could not create: " + calibFileName);
	//err = srcDataSet->SetMetadata(pMetaData, "TIFFTAG_GDAL_METADATA");
	//store param info in the file - use "gdalinfo -mdd "TIFFTAG_GDAL_METADATA"" to see this 
	std::string tag = "Radiom. Homogenised - Win: (" + std::to_string(winSize[0]) + "," + std::to_string(winSize[1]) +
		") Model: " + ModelFormsToString(modelForm) +" Ref: " + refFileName.substr(refFileName.find_last_of("/\\") + 1);
	err = calibDataSet->SetMetadataItem("IMAGEDESCRIPTION", tag.c_str());
	err = calibDataSet->SetMetadataItem("TIFFTAG_IMAGEDESCRIPTION", tag.c_str(), "TIFFTAG_GDAL_METADATA");
	calibDataSet->SetDescription(tag.c_str());
	//CSLDestroy(pMetaData);

	calibDataSet->SetGeoTransform(srcGeoTransform);
	calibDataSet->SetProjection(srcDataSet->GetProjectionRef());
	res = calibDataSet->SetMetadataItem("BitsPerSample", "12", NULL);

	//find ROI of source image in upsampled param image
	double paramInvGeoTransform[6], paramGeoTransform[6];
	paramDataSet->GetGeoTransform(paramGeoTransform);
	GDALInvGeoTransform(paramGeoTransform, paramInvGeoTransform);

	double paramUlI, paramUlJ, paramBrI, paramBrJ;  

	GDALApplyGeoTransform(srcGeoTransform, 0, 0, &ulX, &ulY);
	GDALApplyGeoTransform(srcGeoTransform, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterYSize(), &brX, &brY);

	GDALApplyGeoTransform(paramInvGeoTransform, ulX, ulY, &paramUlJ, &paramUlI);
	GDALApplyGeoTransform(paramInvGeoTransform, brX, brY, &paramBrJ, &paramBrI);

	int nSrcBands = srcDataSet->GetRasterCount();
	nRows = srcDataSet->GetRasterYSize(); nCols = srcDataSet->GetRasterXSize();

	cv::Mat_<float> srcData(nSrcBands, nCols);
	cv::Mat_<float> paramSubData(nParamBands, nCols);
	cv::Mat_<float> calibData(nSrcBands, nCols);
	cv::Mat_<float> imageOffsetMat(nSrcBands, nCols, 0.);
	for (int k = 0; k < nSrcBands; k++)
		imageOffsetMat.row(k).setTo(imageOffset[k]);
	cv::Mat_<float> imageGainMat(nSrcBands, nCols, 1.);
	for (int k = 0; k < nSrcBands; k++)
		imageGainMat.row(k).setTo(imageGain[k]);

	int prog = 0;
	int updateProgNum = 0;

#if SEAMLINE_EXTRAP_FIX
	//upsample eroded mask to source resolution
#if XCALIB_DEBUG
	string erodeMaskUsFileName = srcMaskFileName.substr(0, srcFileName.length() - 4) + "_MASK_US_ERODE.tif";
#else
	string erodeMaskUsFileName = string(CPLGetCurrentDir()) + "\\ErodeMaskUs.tif";
#endif
	std::cout << "Upsampling eroded mask (" << erodeMaskUsFileName << ")" << endl << endl;
#if BIGTIFF
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-co~BIGTIFF=YES~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~None~-wm~2048~-r~bilinear~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), erodeMaskDsFileName.c_str(), erodeMaskUsFileName.c_str());
#else
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~None~-wm~2048~-r~bilinear~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), erodeMaskDsFileName.c_str(), erodeMaskUsFileName.c_str());
#endif

	res = GdalWarpWrapper(gdalString);
	if (res != 0)
		throw string("Could not warp " + erodeMaskDsFileName);

	GDALDataset* erodeMaskDataSet = (GDALDataset*)GDALOpen(erodeMaskUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (erodeMaskDataSet == NULL)
		throw string("Could not open: " + erodeMaskUsFileName);

	//Buf3d<unsigned char> erodeMaskSubData(1, srcDataSet->GetRasterXSize(), 1);
	cv::Mat_<unsigned char> erodeMaskSubData(1, nCols, 0.);

	/*GDALDataset* erodeMaskDataSet = (GDALDataset*)GDALOpen(erodeMaskUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (erodeMaskDataSet == NULL)
	throw string("Could not open: " + erodeMaskUsFileName);*/
#endif
	//TODO: consider rewriting to make sure it is by block (xcalib block size is 256x256), src and param blocks are rows
	//TODO would gdal_calc do this any faster?  if so, how
	//process by row (usually same size as block so should be fast)
	for (int i = 0; i < nRows; i++)
	{
		//TODO: check the implications of storage pixel vs band ordering i.e. in terms of RasterIO speed i.e. is the data stored on disk in the same way we are reading it?
#pragma omp parallel sections  // do disk read/writes in parallel in case files are on different disks
		{
#pragma omp section
			{
				err = paramDataSet->RasterIO(GDALRWFlag::GF_Read, (int)paramUlJ, (int)paramUlI + i, nCols, 1,
					paramSubData.data, nCols, 1, GDALDataType::GDT_Float32,
					nParamBands, NULL, 0, 0, 0);
				if (err != CPLErr::CE_None)
					throw string("paramDataSet->RasterIO failure");
			}
#pragma omp section
			{
				err = srcDataSet->RasterIO(GDALRWFlag::GF_Read, (int)0, (int)i, nCols, 1,
					srcData.data, nCols, 1, GDALDataType::GDT_Float32,
					nSrcBands, NULL, 0, 0, 0);
				if (err != CPLErr::CE_None)
					throw string("srcDataSet->RasterIO failure");
			}
#if SEAMLINE_EXTRAP_FIX  
#pragma omp section
			{
				err = erodeMaskDataSet->RasterIO(GDALRWFlag::GF_Read, (int)paramUlJ, (int)paramUlI + i, nCols, 1,
					erodeMaskSubData.data, nCols, 1, GDALDataType::GDT_Byte,
					1, NULL, 0, 0, 0); //nLineSpace is wrong I think but not used (?)
				if (err != CPLErr::CE_None)
					throw string("erodeMaskDataSet->RasterIO failure");
			}
#endif
#pragma omp section
			if (i > 0) // write calibData from the previous loop in parallel with above reads 
			{
				err = calibDataSet->RasterIO(GDALRWFlag::GF_Write, (int)0, (int)i, nCols, 1,
					calibData.data, nCols, 1, GDALDataType::GDT_Float32,
					nSrcBands, NULL, 0, 0, 0);

				if (err != CPLErr::CE_None)
					throw string("calibDataSet->RasterIO failure");
			}
		}
#pragma omp barrier
		{
			//cv::Mat nodataMask = cv::Mat_<float>(paramSubData(cv::Range(0, nSrcBands), cv::Range::all()) != 0. & srcData != 0) / 255.;  
			if (modelForm == ModelForms::GAIN_ONLY || modelForm == ModelForms::GAIN_AND_IMAGE_OFFSET)
				calibData = paramSubData.mul(srcData + imageOffsetMat); // +imageOffset[k]);
			else if (modelForm == ModelForms::GAIN_AND_OFFSET)
			{
				cv::Mat_<float> M = paramSubData(cv::Range(0, nSrcBands), cv::Range::all());
				cv::Mat_<float> C = paramSubData(cv::Range(nSrcBands, nParamBands), cv::Range::all());
				//calibData = (M.mul(srcData) + C).mul(nodataMask);
				calibData = M.mul(srcData) + C;
			}
			else if (modelForm == ModelForms::OFFSET_ONLY || modelForm == ModelForms::IMAGE_GAIN_AND_OFFSET)
			{
				calibData = srcData.mul(imageGainMat) + paramSubData;
				calibData.setTo(0, paramSubData == 0.);
			}
#if TRUE //SLIDE_WIN_CENTER  this can also be necessary for SLIDE_WIN_CENTER = 0 eg partial pixels have param vals because of win size
			//cv::Mat nodataMask = cv::Mat_<float>(srcData == 0) / 255.;  // necessary for SLIDE_WIN_CENTER which usually has no nodata in paramDs
			//calibData = calibData.mul(nodataMask);
			calibData.setTo(0, srcData == 0);  // necessary for SLIDE_WIN_CENTER which usually has no nodata in paramDs
#endif
#if SEAMLINE_EXTRAP_FIX
			cv::Mat erodeMask = erodeMaskSubData == 0;
			for (int k = 0; k < nSrcBands; k++)
				calibData.row(k).setTo(0, erodeMask);
#endif
		}

		if ((40*(i+1)) / nRows > prog)
		{
			prog = (40*(i+1)) / nRows;
			if (2.5*prog > updateProgNum)
			{
				std::cout << updateProgNum;
				updateProgNum += 10;
			}
			else
				std::cout << ".";
		}
	}
	//write out calibData from final loop
	{
		err = calibDataSet->RasterIO(GDALRWFlag::GF_Write, (int)0, (int)nRows-1, nCols, 1,
			calibData.data, nCols, 1, GDALDataType::GDT_Float32,
			nSrcBands, NULL, 0, 0, 0);

		if (err != CPLErr::CE_None)
			throw string("calibDataSet->RasterIO failure");
	}
	std::cout << endl << endl;
		
	for (int k = 0; k < nSrcBands; k++)
	{
		GDALRasterBand* calibBand = calibDataSet->GetRasterBand(k+1);
		res = calibBand->SetNoDataValue(0);
	}

	calibDataSet->FlushCache();
	GDALClose(calibDataSet);
	GDALClose(srcDataSet);
	GDALClose(paramDataSet);
#if SEAMLINE_EXTRAP_FIX
	GDALClose(erodeMaskDataSet);
#endif
	CSLDestroy(papszOptions);
	//CSLDestroy(pMetaData);

	return res;
}

//-------------------------------------------------------------------------------------------------
// Usage
//-------------------------------------------------------------------------------------------------
void Usage()
{
	std::cout << endl << "Usage: CrossCalibration.exe [-o] [-w height width] [-p num params] [ref file name] [src file name]" << endl << endl;
	std::cout << "Source and reference files must be tiffs in the same co-ordinates and with matching bands. Reference image must contain the source area." << endl << endl;
	std::cout << "[ref file name]: name of a calibrated reference file" << endl;
	std::cout << "[src file name]: name of the file to be calibrated" << endl;
	std::cout << "-o: overwrite calibrated file if it exists otherwise exit (optional - default off)" << endl;
	std::cout << "-w: specify the height and width of the sliding window in reference pixels (optional - default 1 1)" << endl;
	std::cout << "-p: specify the linear model form (1 = gain only, 2 = gain and offset, 3 = offset only)" << endl;
	exit(-1);
}

//-------------------------------------------------------------------------------------------------
// Main
//-------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	int res = 0;
	bool doOverwrite = false;
	int winSize[2] = {1,1};
    char* pSrcFileName = NULL;
    char* pRefFilename = NULL;
	int modelForm = ModelForms::GAIN_ONLY;

	//parse arguments
	if (argc < 3)
		Usage();

	for (int i = 1; i < argc; i++)
	{
		if (_strcmpi(argv[i], "-o") == 0)
			doOverwrite = true;
		else if (_strcmpi(argv[i], "-w") == 0 && i < argc-2)
		{
            winSize[0] = atoi(argv[++i]);
            winSize[1] = atoi(argv[++i]);
		}
		else if (_strcmpi(argv[i], "-p") == 0 && i < argc - 1)
		{
			modelForm = atoi(argv[++i]);
		}
		else if (argv[i][0] == '-')
			Usage();
		else if (i == argc-2)
			pRefFilename = argv[i];
		else if (i == argc-1)
			pSrcFileName = argv[i];
	}

    if (pRefFilename == NULL || pSrcFileName == NULL)
        Usage();

	string srcFileName(pSrcFileName);
	string refFileName(pRefFilename);

	if (!doOverwrite) //check if out calib file exists
	{
#if XCALIB_DEBUG
		string calibFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_INIT_XCALIB.tif";
#else
		string calibFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_XCALIB.tif";
#endif
		//char tmp[MAX_PATH];
		//strcpy(tmp, calibFileName.c_str());
		if (CPLCheckForFile((char*)calibFileName.c_str(), NULL) == 1)
		{
			std::cout << calibFileName  <<  " exists.  Exiting" << endl;
			return -1;
		}
	}

	GDALAllRegister();
#if XCALIB_DEBUG 
	//CPLSetConfigOption("CPL_DEBUG", "ON");  //hack
#else
	CPLSetConfigOption("CPL_DEBUG", "OFF");
#endif
	CPLSetConfigOption("GDAL_CACHEMAX", "2048");//?
	try
	{
		res = CrossCalib(refFileName, srcFileName, winSize, modelForm);
	}
	catch(string ex) //catch exceptions thrown by CrossCalib
	{
		cerr << "Error: " << ex << endl;
	}

#if XCALIB_DEBUG 
	/*CPLErrorReset();
	const char  *pszDebug = CPLGetConfigOption("CPL_DEBUG", NULL);
	if (pszDebug && (EQUAL(pszDebug,"ON") || EQUAL(pszDebug,"")))
	{  
		GDALDumpOpenDatasets(stderr); //access violation?
		CPLDumpSharedList(NULL);
	}*/
#endif
	GDALDestroyDriverManager();

	return res;
}

