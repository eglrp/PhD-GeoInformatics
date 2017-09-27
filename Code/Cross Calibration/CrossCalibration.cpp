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


#include "CrossCalibration.h"

using namespace std;

#define XCALIB_DEBUG 1
#define MAX_PATH 1024
#define SEAMLINE_FIX 1
#define DO_COMPRESS_OUTPUT 1
#define BIGTIFF 1

int defaultWinSize[2] = {1, 1};
int gdalwarp(int argc, char ** argv);
int gdal_translate(int argc, char ** argv);

//-------------------------------------------------------------------------------------------------
//Parses argString and calls gdalwarp backend
//-------------------------------------------------------------------------------------------------
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
int CrossCalib(const string& refFileName, const string& srcFileName_, const int* winSize = defaultWinSize)
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
	
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-overwrite~-srcnodata~\"0\"~-dstnodata~\"0\"~-wm~1024~-r~bilinear~-tap~-tr~%1.1f~%1.1f~%s~%s~", 2.5, 2.5, 
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

	/*CPLErr err = srcDataSet->GetRasterBand(1)->CreateMaskBand(GMF_NODATA | GMF_PER_DATASET);
	if (err != CPLErr::CE_None)
		throw string("srcDsDataSet->CreateMaskBand failure");
	GDALRasterBand* maskBand = srcDataSet->GetRasterBand(1)->GetMaskBand();
	if (maskBand == NULL)
		throw string("Could not GetMaskBand");
	//maskBand->RasterIO()*/

	srcDataSet->GetGeoTransform(srcGeoTransform);
	int srcXSize = srcDataSet->GetRasterXSize();
	int srcYSize = srcDataSet->GetRasterYSize();
	GDALClose(srcDataSet);

	cout << "------------------------------------------------------------" << endl;
	cout << "Processing: " << srcFileName << endl << endl;
#if SEAMLINE_FIX
	cout << "Extracting mask from: " << srcFileName << endl << endl;
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

	cout << "Downsampling mask: " << srcMaskFileName << endl << endl;

#if XCALIB_DEBUG
	string srcMaskDsFileName = srcMaskFileName.substr(0, srcMaskFileName.length() - 4) + "_DS.tif"; //make separate files for each source file
#else
	string srcMaskDsFileName = string(CPLGetCurrentDir()) + "\\SourceMaskDs.tif"; //reuse the same file
#endif

	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-overwrite~-srcnodata~none~-dstnodata~none~-r~average~-tap~-tr~%d~%d~%s~%s~", (int)abs(refGeoTransform[1]), (int)abs(refGeoTransform[5]),
		srcMaskFileName.c_str(), srcMaskDsFileName.c_str());
	res = GdalWarpWrapper(gdalString);

	if (res != 0)
		throw string("Could not warp " + srcMaskFileName);

#endif //SEAMLINE_FIX

	cout << "Downsampling: " << srcFileName << endl << endl;

#if XCALIB_DEBUG
	string srcDsFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_DS.tif"; //make separate files for each source file
#else
	string srcDsFileName = string(CPLGetCurrentDir()) + "\\SourceDs.tif"; //reuse the same file
#endif
	//
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-co~\"BIGTIFF=YES\"~-overwrite~-srcnodata~\"0\"~-dstnodata~\"0\"~-r~cubicspline~-tap~-tr~%d~%d~%s~%s~", (int)abs(refGeoTransform[1]), (int)abs(refGeoTransform[5]), 
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

	GDALApplyGeoTransform(refInvGeoTransform, ulX, ulY, &refUlI, &refUlJ);
	GDALApplyGeoTransform(refInvGeoTransform, brX, brY, &refBrI, &refBrJ);

	//--------------------------------------------------------------------------------------------------------------
	//Find the calibration gains
	//-------------------------------------------------------------------------------------------------

	Buf3d<double> srcDsData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());
	Buf3d<double> refSubData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());
	Buf3d<double> gainDsData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());

#if XCALIB_DEBUG
	string gainDsFileName = srcDsFileName.substr(0, srcDsFileName.length() - 4) + "_GAIN.tif"; //make separate files for each source file
#else
	string gainDsFileName = string(CPLGetCurrentDir()) + "\\GainDs.tif"; //reuse the same file
#endif

	cout << "Calculating calibration gains (" << gainDsFileName << ")" << endl << endl;

	GDALDataset* gainDsDataSet = srcDsDataSet->GetDriver()->Create(gainDsFileName.c_str(), srcDsDataSet->GetRasterXSize(), 
		srcDsDataSet->GetRasterYSize(), 4, GDALDataType::GDT_Float64, NULL);

	if (gainDsDataSet == NULL)
		throw string("Could not create: " + gainDsFileName);

	gainDsDataSet->SetGeoTransform(srcDsGeoTransform);
	gainDsDataSet->SetProjection(srcDsDataSet->GetProjectionRef());	

	CPLErr err = srcDsDataSet->RasterIO(GDALRWFlag::GF_Read, 0, 0, srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(),
		srcDsData.Buf(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float64, 
		srcDsDataSet->GetRasterCount(), NULL, 0, 0, 0);
	
	if (err != CPLErr::CE_None)
		throw string("srcDsDataSet->RasterIO failure");

	err = refDataSet->RasterIO(GDALRWFlag::GF_Read, (int)refUlI, (int)refUlJ, srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(),
		refSubData.Buf(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float64, 
		srcDsDataSet->GetRasterCount(), NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("refDataSet->RasterIO failure");

#if SEAMLINE_FIX
	//read the ds mask
	GDALDataset* srcMaskDsDataSet = (GDALDataset*)GDALOpen(srcMaskDsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcMaskDsDataSet == NULL)
		throw string("Could not open: " + srcMaskDsFileName);
	Buf3d<double> srcMaskDsData(srcDsDataSet->GetRasterYSize(), srcDsDataSet->GetRasterXSize(), srcDsDataSet->GetRasterCount());

	err = srcMaskDsDataSet->RasterIO(GDALRWFlag::GF_Read, 0, 0, srcMaskDsDataSet->GetRasterXSize(), srcMaskDsDataSet->GetRasterYSize(),
		srcMaskDsData.Buf(), srcMaskDsDataSet->GetRasterXSize(), srcMaskDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float64,
		srcMaskDsDataSet->GetRasterCount(), NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("srcMaskDsData->RasterIO failure");
	//make eroded mask for application to upsampled images
	//string winTitle = "mask";
	cv::Mat cvMask(srcMaskDsDataSet->GetRasterYSize(), srcMaskDsDataSet->GetRasterXSize(), CV_64FC1, (void*)srcMaskDsData.Buf());
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
		throw string("Could not create: " + gainDsFileName);

	erodeMaskDsDataSet->SetGeoTransform(srcDsGeoTransform);
	erodeMaskDsDataSet->SetProjection(srcDsDataSet->GetProjectionRef());

	err = erodeMaskDsDataSet->RasterIO(GDALRWFlag::GF_Write, 0, 0, erodeMaskDsDataSet->GetRasterXSize(), erodeMaskDsDataSet->GetRasterYSize(),
		cvMaskErode.data, erodeMaskDsDataSet->GetRasterXSize(), erodeMaskDsDataSet->GetRasterYSize(), GDALDataType::GDT_Byte,
		1, NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("erodeMaskDsDataSet->RasterIO failure");

	erodeMaskDsDataSet->FlushCache();
	GDALClose(erodeMaskDsDataSet);
	GDALClose(srcMaskDsDataSet);

#endif //#if SEAMLINE_FIX


	//gainDsDataSet->GetRasterBand(1)->Mas

	int n = srcDsDataSet->GetRasterYSize(), m = srcDsDataSet->GetRasterXSize();
	for (int k = 0; k < srcDsDataSet->GetRasterCount(); k++)
	{
		err = gainDsDataSet->GetRasterBand(k+1)->SetNoDataValue(0);
		for (int j = 0; j < m; j++)
		{
			//average gain over the sliding window
			int colIdx[2] = {max<int>(j-(winSize[1]-1)/2, 0), min<int>(m-1, j+(winSize[1]-1)/2)};
			for (int i = 0; i < n; i++) 			//omit the border pixels so that we are only using ds pixels fully covered by src pixels
			{
				int rowIdx[2] = {max<int>(i-(winSize[0]-1)/2, 0), min<int>(n-1, i+(winSize[0]-1)/2)};
				double sum = 0;
				int count = 0;
				for (int jw = colIdx[0]; jw <= colIdx[1]; jw++)
				{
					for (int iw = rowIdx[0]; iw <= rowIdx[1]; iw++)
					{
						if (srcDsData(iw, jw, k) != 0 && refSubData(iw, jw, k) != 0) //ignore nodata
						{
							sum += (double)refSubData(iw, jw, k)/(double)srcDsData(iw, jw, k);
							count++;
						}
					}
				}
#if SEAMLINE_FIX
				if (count > 0 && srcMaskDsData(i, j, 0) >= 0.95*255.0)
#else
				if (count > 0)
#endif
					gainDsData(i, j, k) = sum / (double)count;
				else
					gainDsData(i, j, k) = 0;
			}
		}
	}

	//write out gains to file
	err = gainDsDataSet->RasterIO(GDALRWFlag::GF_Write, 0, 0, gainDsDataSet->GetRasterXSize(), gainDsDataSet->GetRasterYSize(),
		gainDsData.Buf(), gainDsDataSet->GetRasterXSize(), gainDsDataSet->GetRasterYSize(), GDALDataType::GDT_Float64, 
		srcDsDataSet->GetRasterCount(), NULL, 0, 0, 0);

	if (err != CPLErr::CE_None)
		throw string("gainDsDataSet->RasterIO failure");

	gainDsDataSet->FlushCache();
	GDALClose(gainDsDataSet);

	GDALClose(srcDsDataSet);
	GDALClose(refDataSet);
	//--------------------------------------------------------------------------------------------------------------
	//Upsample gain to source resolution
	//-------------------------------------------------------------------------------------------------

#if XCALIB_DEBUG
	string gainUsFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_US_GAIN.tif";
#else
	string gainUsFileName = string(CPLGetCurrentDir()) + "\\GainUs.tif";
#endif
	cout << "Interpolating calibration gains (" << gainUsFileName << ")" << endl << endl;

#if BIGTIFF
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-co~\"BIGTIFF=YES\"~-multi~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), gainDsFileName.c_str(), gainUsFileName.c_str());
#else
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), gainDsFileName.c_str(), gainUsFileName.c_str());
#endif

	res = GdalWarpWrapper(gdalString);
	if (res != 0)
		throw string("Could not warp " + gainDsFileName);


	//-------------------------------------------------------------------------------------------------
	//Apply umsampled gains to source image
	//-------------------------------------------------------------------------------------------------

	string calibFileName = srcFileName.substr(0, srcFileName.length() - 4) + "_XCALIB.tif";
	cout << "Applying gains (" << calibFileName << ")" << endl << endl;

	srcDataSet = (GDALDataset*)GDALOpen(srcFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (srcDataSet == NULL)
		throw string("Could not open: " + srcFileName);

	GDALDataset* gainDataSet = (GDALDataset*)GDALOpen(gainUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (gainDataSet == NULL)
		throw string("Could not open: " + gainUsFileName);
	char **papszOptions = NULL;
#if DO_COMPRESS_OUTPUT
	papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
	papszOptions = CSLSetNameValue(papszOptions, "PREDICTOR", "2");
#endif
#if BIGTIFF
	papszOptions = CSLSetNameValue(papszOptions, "BIGTIFF", "YES");
#endif
	GDALDataset* calibDataSet = srcDataSet->GetDriver()->Create(calibFileName.c_str(), srcDataSet->GetRasterXSize(),
		srcDataSet->GetRasterYSize(), 4, GDALDataType::GDT_UInt16, papszOptions);

	if (calibDataSet == NULL)
		throw string("Could not create: " + calibFileName);

	calibDataSet->SetGeoTransform(srcGeoTransform);
	calibDataSet->SetProjection(srcDataSet->GetProjectionRef());
	res = calibDataSet->SetMetadataItem("BitsPerSample", "12", NULL);

	//find ROI of source image in upsampled gain image
	double gainInvGeoTransform[6], gainGeoTransform[6];
	gainDataSet->GetGeoTransform(gainGeoTransform);
	GDALInvGeoTransform(gainGeoTransform, gainInvGeoTransform);

	double gainUlI, gainUlJ, gainBrI, gainBrJ;

	GDALApplyGeoTransform(srcGeoTransform, 0, 0, &ulX, &ulY);
	GDALApplyGeoTransform(srcGeoTransform, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterYSize(), &brX, &brY);

	GDALApplyGeoTransform(gainInvGeoTransform, ulX, ulY, &gainUlI, &gainUlJ);
	GDALApplyGeoTransform(gainInvGeoTransform, brX, brY, &gainBrI, &gainBrJ);

	Buf3d<double> srcData(1, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterCount());
	Buf3d<double> gainSubData(1, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterCount());
	Buf3d<double> calibData(1, srcDataSet->GetRasterXSize(), srcDataSet->GetRasterCount());
	int prog = 0;

#if SEAMLINE_FIX
	//upsample eroded mask to source resolution
#if XCALIB_DEBUG
	string erodeMaskUsFileName = srcMaskFileName.substr(0, srcFileName.length() - 4) + "_MASK_US_ERODE.tif";
#else
	string erodeMaskUsFileName = string(CPLGetCurrentDir()) + "\\ErodeMaskUs.tif";
#endif
	cout << "Upsampling eroded mask (" << erodeMaskUsFileName << ")" << endl << endl;
#if BIGTIFF
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-co~\"BIGTIFF=YES\"~-multi~-overwrite~-srcnodata~0~-dstnodata~None~-wm~2048~-r~bilinear~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), erodeMaskDsFileName.c_str(), erodeMaskUsFileName.c_str());
#else
	sprintf_s(gdalString, MAX_PATH, "gdalwarp~-multi~-overwrite~-srcnodata~0~-dstnodata~None~-wm~2048~-r~bilinear~-tr~%f~%f~%s~%s~",
		fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), erodeMaskDsFileName.c_str(), erodeMaskUsFileName.c_str());
#endif

	res = GdalWarpWrapper(gdalString);
	if (res != 0)
		throw string("Could not warp " + erodeMaskDsFileName);

	GDALDataset* erodeMaskDataSet = (GDALDataset*)GDALOpen(erodeMaskUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (erodeMaskDataSet == NULL)
		throw string("Could not open: " + erodeMaskUsFileName);

	Buf3d<unsigned char> erodeMaskSubData(1, srcDataSet->GetRasterXSize(), 1);
	/*GDALDataset* erodeMaskDataSet = (GDALDataset*)GDALOpen(erodeMaskUsFileName.c_str(), GDALAccess::GA_ReadOnly);
	if (erodeMaskDataSet == NULL)
	throw string("Could not open: " + erodeMaskUsFileName);*/
#endif

	//process by row (usually same size as block so should be fast)
	for (int j = 0; j < srcDataSet->GetRasterYSize(); j++)
	{
		err = gainDataSet->RasterIO(GDALRWFlag::GF_Read, (int)gainUlI, (int)gainUlJ+j, srcDataSet->GetRasterXSize(), 1,
			gainSubData.Buf(), srcDataSet->GetRasterXSize(), 1, GDALDataType::GDT_Float64, 
			srcDataSet->GetRasterCount(), NULL, 0, 0, 0);
		if (err != CPLErr::CE_None)
			throw string("gainDataSet->RasterIO failure");
			
		err = srcDataSet->RasterIO(GDALRWFlag::GF_Read, (int)0, (int)j, srcDataSet->GetRasterXSize(), 1,
			srcData.Buf(), srcDataSet->GetRasterXSize(), 1, GDALDataType::GDT_Float64, 
			srcDataSet->GetRasterCount(), NULL, 0, 0, 0);
		if (err != CPLErr::CE_None)
			throw string("srcDataSet->RasterIO failure");

#if SEAMLINE_FIX
		err = erodeMaskDataSet->RasterIO(GDALRWFlag::GF_Read, (int)gainUlI, (int)gainUlJ + j, srcDataSet->GetRasterXSize(), 1,
			erodeMaskSubData.Buf(), srcDataSet->GetRasterXSize(), 1, GDALDataType::GDT_Byte,
			1, NULL, 0, 0, 0); //nLineSpace is wrong I think but not used (?)
		if (err != CPLErr::CE_None)
			throw string("erodeMaskDataSet->RasterIO failure");
#endif
		for (int k = 0; k < srcDataSet->GetRasterCount(); k++)
		{
			for (int i = 0; i < srcDataSet->GetRasterXSize(); i++)
			{
#if SEAMLINE_FIX
				if (erodeMaskSubData(0, i, 0) > (unsigned char)0.95*255.0)
					calibData(0, i, k) = gainSubData(0, i, k) * srcData(0, i, k);
				else
#endif
					calibData(0, i, k) = 0;
			}
		}
		err = calibDataSet->RasterIO(GDALRWFlag::GF_Write, (int)0, (int)j, srcDataSet->GetRasterXSize(), 1,
			calibData.Buf(), srcDataSet->GetRasterXSize(), 1, GDALDataType::GDT_Float64, 
			srcDataSet->GetRasterCount(), NULL, 0, 0, 0);
		
		if (err != CPLErr::CE_None)
			throw string("calibDataSet->RasterIO failure");

		if ((100*j) / srcDataSet->GetRasterYSize() > prog)
		{
			prog = (100*j) / srcDataSet->GetRasterYSize();
			cout << ".";
		}
			
	}
	cout << endl << endl;
		
	for (int k = 0; k < calibDataSet->GetRasterCount(); k++)
	{
		GDALRasterBand* calibBand = calibDataSet->GetRasterBand(k+1);
		res = calibBand->SetNoDataValue(0);
	}

	calibDataSet->FlushCache();
	GDALClose(calibDataSet);
	GDALClose(srcDataSet);
	GDALClose(gainDataSet);
#if SEAMLINE_FIX
	GDALClose(erodeMaskDataSet);
#endif
	return res;
}

//-------------------------------------------------------------------------------------------------
// Usage
//-------------------------------------------------------------------------------------------------
void Usage()
{
	cout << endl << "Usage: CrossCalibration.exe [-o] [-w width height] [ref file name] [src file name]" << endl << endl;
	cout << "Source and reference files must be tiffs in the same co-ordinates and with matching bands. Reference image must contain the source area." << endl << endl;
	cout << "[ref file name]: name of a calibrated reference file" << endl;
	cout << "[src file name]: name of the file to be calibrated" << endl;
	cout << "-o: overwrite calibrated file if it exists otherwise exit (optional - default off)" << endl;
	cout << "-w: specify the width and height of the sliding window in reference pixels (optional - default 1 1)" << endl;
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

	//parse arguments
	if (argc < 3)
		Usage();

	for (int i = 1; i < argc; i++)
	{
		if (EQUAL(argv[i], "-o"))
			doOverwrite = true;
		else if (EQUAL(argv[i], "-w") && i < argc-2)
		{
            winSize[0] = atoi(argv[++i]);
            winSize[1] = atoi(argv[++i]);
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
			cout << calibFileName  <<  " exists.  Exiting" << endl;
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
		res = CrossCalib(refFileName, srcFileName, winSize);
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

