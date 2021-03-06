/**
* A simple abstraction layer which uses GDAL to open
* image formats not supported by OpenCV.
*
* @author Marvin Smith
* @date 21 March 2013
*/

// OpenCV Dependencies

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

// GDAL Dependencies
#include <cpl_conv.h>
#include <gdal_priv.h>

// Boost Dependencies
//#include <boost/filesystem.hpp>

// STL Dependencies
#include <string>
#include <iostream>

using namespace cv;
using namespace std;

#if 1
/**
* Compute the scale factor for the conversion between image depths
*
* @param[in] gdalDepth - GDAL Depth Type
* @param[in] opencvDepth - OpenCV Depth Type
* @return Scale factor between depths
*/
double gdal2OpenCVScale( const int& gdalDepth, const int& opencvDepth ){

	if( opencvDepth == CV_8U && gdalDepth == GDT_Byte ) return 1;
	if( opencvDepth == CV_8U && gdalDepth == GDT_Int16 ) return 1/16.0;
	if( opencvDepth == CV_8U && gdalDepth == GDT_UInt16 ) return 1/16.0;
	if( opencvDepth == CV_16U && gdalDepth == GDT_Byte ) return 16;
	if( opencvDepth == CV_16U && gdalDepth == GDT_Int16 ) return 16;
	if( opencvDepth == CV_16U && gdalDepth == GDT_UInt16 ) return 1;
	if( opencvDepth == CV_16S && gdalDepth == GDT_Byte ) return 16;
	if( opencvDepth == CV_16S && gdalDepth == GDT_Int16 ) return 16;
	if( opencvDepth == CV_16S && gdalDepth == GDT_UInt16 ) return 1;

	throw string("Error: Unknown OpenCV Type or Unknown GDAL Type");
}

/**
* Convert and set the value from the gdal raster to the opencv image
*
* @param[in] val Pixel value from GDAL dataset
* @param[in] image OpenCV Image to load
* @param[in] point Pixel to set the value to
* @param[in] scale Scale factor to convert with
*/
void GDAL2OpenCV( double const& val, Mat& image, Point const& point, const double& scale ){

	// get the required depth
	int cvDepth = image.depth();

	// get the value for the image
	double pixVal = val * scale;

	// place the value
	if( image.depth() == CV_8U ){
		image.at<uchar>(point) = (uchar)pixVal;
	}
	else if( image.depth() == CV_16U ){
		image.at<ushort>(point) = (ushort)pixVal;
	}
	else if( image.depth() == CV_16S ){
		image.at<short>(point) = (short)pixVal;
	}
	else if( image.depth() == CV_32S ){
		image.at<int>(point) = (int)pixVal;
	}
	else if( image.depth() == CV_32F ){
		image.at<float>(point) = (float)pixVal;
	}
	else if( image.depth() == CV_64F ){
		image.at<double>(point) = (double)pixVal;
	}
	else{
		throw string("Error: Unknown Depth");
	}

}

/**
* Convert Type 2 Depth
*
* @param[in] type OpenCV Type to convert
* @return Associated depth
*/
int cvType2Depth( const int& type ){

	switch( type ){

	case CV_8UC1:
		return CV_8U;
	case CV_8UC2:
		return CV_8U;
	case CV_8UC3:
		return CV_8U;
	case CV_16UC1:
		return CV_16U;
	case CV_16UC2:
		return CV_16U;
	case CV_16UC3:
		return CV_16U;
	case CV_16SC1:
		return CV_16S;
	case CV_16SC2:
		return CV_16S;
	case CV_16SC3:
		return CV_16S;
	case CV_32SC1:
		return CV_32S;
	case CV_32SC2:
		return CV_32S;
	case CV_32SC3:
		return CV_32S;
	case CV_32FC1:
		return CV_32F;
	case CV_32FC2:
		return CV_32F;
	case CV_32FC3:
		return CV_32F;
	case CV_64FC1:
		return CV_64F;
	case CV_64FC2:
		return CV_64F;
	case CV_64FC3:
		return CV_64F;
	default:
		return -1;
	}
}


/**
* Extract the number of channels from the type
*
* @param[in] type
* @return Number of channels
*/
int cvType2Channels( const int& type ){

	switch( type ){

	case CV_8UC1:
	case CV_16UC1:
	case CV_16SC1:
	case CV_32SC1:
	case CV_32FC1:
	case CV_64FC1:
		return 1;
	case CV_8UC2:
	case CV_16UC2:
	case CV_16SC2:
	case CV_32SC2:
	case CV_32FC2:
	case CV_64FC2:
		return 2;
	case CV_8UC3:
	case CV_16UC3:
	case CV_16SC3:
	case CV_32SC3:
	case CV_32FC3:
	case CV_64FC3:
		return 3;

	default:
		return 0;
	}
}

/**
* Convert Depth and Channels into Type
*
* @param[in] Depth OpenCV Depth value
* @param[in] Channels Number of channels
* @return OpenCV Type
*/
int cvDepthChannel2Type( const int Depth, const int Channels ){

	if( Depth == CV_8U && Channels == 3 ) return CV_8UC3;
	if( Depth == CV_8U && Channels == 2 ) return CV_8UC2;
	if( Depth == CV_8U && Channels == 1 ) return CV_8UC1;
	if( Depth == CV_16U && Channels == 3 ) return CV_16UC3;
	if( Depth == CV_16U && Channels == 2 ) return CV_16UC2;
	if( Depth == CV_16U && Channels == 1 ) return CV_16UC1;
	if( Depth == CV_16S && Channels == 3 ) return CV_16SC3;
	if( Depth == CV_16S && Channels == 2 ) return CV_16SC2;
	if( Depth == CV_16S && Channels == 1 ) return CV_16SC1;
	if( Depth == CV_32S && Channels == 3 ) return CV_32SC3;
	if( Depth == CV_32S && Channels == 2 ) return CV_32SC2;
	if( Depth == CV_32S && Channels == 1 ) return CV_32SC1;
	if( Depth == CV_64F && Channels == 1 ) return CV_64FC1;
	//if( Depth == CV_16U && Channels == 4 ) return CV_16UC(4); //dh

	throw string("Error: combo not supported");

	return 0;
}


/**
* Checks if the OpenCV Image Type is supported for this.
*
* 2 channel images are not supported only because I am not sure what
* types of images these are.
*/
bool validOpenCVImageType( const int& imageType )
{
	//int ngiType = CV_16U

	switch( imageType )
	{

	case CV_8UC1:
	case CV_8UC3:
	case CV_16UC1:
	case CV_16UC2:
	case CV_16UC3:
	case CV_16SC1:
	case CV_16SC2:
	case CV_16SC3:
	case CV_32SC1:
	case CV_32SC2:
	case CV_32SC3:
	case CV_32FC1:
	case CV_32FC2:
	case CV_32FC3:
	case CV_64FC1:
	case CV_64FC2:
	case CV_64FC3:
		//case CV_16UC(4): //dh
		return true;
	default:
		return false;
	}
}

int CvDataType(GDALDataset* ds)
{
	switch (ds->GetRasterBand(0)->GetRasterDataType()) //NB: assume bands are all of the same type
	{
	case GDALDataType::GDT_Byte:
		return CV_8UC(ds->GetRasterCount());
	case GDALDataType::GDT_Int16:
		return CV_16SC(ds->GetRasterCount());
	case GDALDataType::GDT_UInt16:
		return CV_16UC(ds->GetRasterCount());
	case GDALDataType::GDT_Int32:
		return CV_32SC(ds->GetRasterCount());
	case GDALDataType::GDT_Float32:
		return CV_32FC(ds->GetRasterCount());
	case GDALDataType::GDT_Float64:
		return CV_64FC(ds->GetRasterCount());
	default:
		throw string("Error: Unsupported data type");
	}
}

int CvDataType(GDALRasterBand* band)
{
	switch (band->GetRasterDataType()) //NB: assume bands are all of the same type
	{
	case GDALDataType::GDT_Byte:
		return CV_8UC1;
	case GDALDataType::GDT_Int16:
		return CV_16SC1;
	case GDALDataType::GDT_UInt16:
		return CV_16UC1;
	case GDALDataType::GDT_Int32:
		return CV_32SC1;
	case GDALDataType::GDT_Float32:
		return CV_32FC1;
	case GDALDataType::GDT_Float64:
		return CV_64FC1;
	default:
		throw string("Error: Unsupported data type");
	}
}


vector<Mat> GdalImread(const string& filename)
{
	// ensure the file exists
	if (CPLCheckForFile((char*)filename.c_str(), NULL) == 0)
		throw string( "Error: File ") + filename + string(" does not exist" );

	// load the dataset
	GDALDataset* dataset = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly);

	// if the dataset returns null, then the file does not exist, or there was a read error
	if (dataset == NULL)
		throw string("Error: GDAL Dataset returned null from read");

	// check if pixel data even exists
	if (dataset->GetRasterCount() <= 0) 
		throw string("Error: File does not contain pixel data");

	// get raster image size
	Size imgSize(dataset->GetRasterXSize(), dataset->GetRasterYSize());
	
	//Reading on block boundaries is more efficient
	int nBands = dataset->GetRasterCount();
	vector<Mat> res(nBands);
	int cacheVal = GDALGetCacheMax();
	GDALSetCacheMax(1024); //for fast reading and GDALCLose (in bytes)
	for (int iBand = 0; iBand < nBands; iBand++) // < nBands; iBand++)
	{
		cout << "Reading band " << iBand << endl;
		int nXBlocks, nYBlocks, nXBlockSize, nYBlockSize;

		GDALRasterBand* band = dataset->GetRasterBand(iBand+1);
		//int nXSize = band->GetXSize(), nYSize = band->GetYSize();
		band->GetBlockSize(&nXBlockSize, &nYBlockSize);
		nXBlocks = (band->GetXSize() + nXBlockSize - 1) / nXBlockSize;
		nYBlocks = (band->GetYSize() + nYBlockSize - 1) / nYBlockSize;
		
		GDALDataType gdalDataType = band->GetRasterDataType(); //NB: assume bands all have same data type
		
		res[iBand] = Mat(band->GetYSize(), band->GetXSize(), CvDataType(band));
		for (int iYBlock = 0; iYBlock < nYBlocks; iYBlock++)
		{
			for (int iXBlock = 0; iXBlock < nXBlocks; iXBlock++)
			{
				if (nYBlockSize == 1) //read row by row
				{
					void* ptr = res[iBand].ptr(iYBlock);
					CPLErr err = band->ReadBlock(iXBlock, iYBlock, ptr);
					if (err != CPLErr::CE_None)
					{
						throw string("Error: RasterIO failure");
					}
				}
				else //make roi header for each row
				{
					Mat roi(res[iBand], cv::Rect(iXBlock*nXBlockSize, iYBlock*nYBlockSize, nXBlockSize, nYBlockSize));
					CPLErr err = band->ReadBlock(iXBlock, iYBlock, roi.data);
					if (err != CPLErr::CE_None)
					{
						throw string("Error: RasterIO failure");
					}
				}
			}
		}
	}
	GDALSetCacheMax(cacheVal);

	GDALClose(dataset);
	return res;
}

/**
* Read an image format using GDAL as the toolset.
*
* @param[in] filename Image filename
* @param[in] flags Type of image you want to read
*/
Mat imread_geo( const string& filename, const int& imageType ){

	// parse the flags to determine which color type we want to read
	if( validOpenCVImageType(imageType) == false ) throw string("Invalid Image Type");

	// ensure the file exists
	/*if( boost::filesystem::exists( filename ) == false ){
	throw string( "Error: File ") + filename + string(" does not exist" );
	}*/
	if (CPLCheckForFile((char*)filename.c_str(), NULL) == 0)
	{
		//cout << calibFileName << " exists.  Exiting" << endl;
		throw string( "Error: File ") + filename + string(" does not exist" );
	}

	// register the gdal driver
	GDALAllRegister();

	// load the dataset
	GDALDataset* dataset = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly);

	// if the dataset returns null, then the file does not exist, or there was a read error
	if( dataset == NULL ){ throw string("Error: GDAL Dataset returned null from read"); }

	// check if pixel data even exists
	if( dataset->GetRasterCount() <= 0 ) throw string("Error: File does not contain pixel data");

	//get the driver infomation
	GDALDriver* driver = dataset->GetDriver();

	// get raster image size
	Size imgSize( dataset->GetRasterXSize(), dataset->GetRasterYSize() );

	// create mats for each raster layer
	//Mat cvImage
	//dataset->Get

	//Mat bob(imgSize, CvDataType(dataset));

	vector<Mat> layers( dataset->GetRasterCount() );
	for( size_t i=0; i<layers.size(); i++ )
		layers[i] = Mat( imgSize, cvDepthChannel2Type(cvType2Depth( imageType), 1));

	// iterate through each band
	for (int i = 0; i < dataset->GetRasterCount(); i++) {

		//create the band object
		GDALRasterBand *band = dataset->GetRasterBand(i + 1);

		// get the gdal band type
		int datatype = band->GetRasterDataType();

		// compute the scale factor
		double scale = gdal2OpenCVScale( datatype, layers[i].depth());

		//read each image row
		for ( int r = 0; r < layers[i].rows; r++) {

			float* pafScanline;
			pafScanline = (float*) CPLMalloc(sizeof (float) *layers[i].cols);
			band->RasterIO(GF_Read, 0, r, layers[i].cols, 1, pafScanline, layers[i].cols, 1, GDT_Float32, 0, 0);

			// iterate through each column
			for ( int c = 0; c < layers[i].cols; c++) {
				GDAL2OpenCV( pafScanline[c], layers[i], Point(c,r), scale );
			}

		}

	}

	// close our GDAL Dataset and clean up
	GDALClose(dataset);

	//merge channels into single image
	Mat output;
	if( layers.size() > 1 )
		cv::merge( layers, output );
	else
		output = layers[0].clone();

	// do channel conversions if necessary
	if( cvType2Channels(imageType) == output.channels() )
		return output;

	if( cvType2Channels(imageType) == 1 && output.channels() == 3 ){
		cvtColor( output, output, CV_BGR2GRAY );
		return output;
	}
	if( cvType2Channels(imageType) == 3 && output.channels() == 1 ){
		cvtColor( output, output, CV_GRAY2BGR);
		return output;
	}

	// Return the output image
	return output;
}

#endif

/**
* Test Program
*/
int main( int argc, char* argv[] )
{
	cv::Mat_<float> tmp(3, 4, 0.);
	cout << tmp << endl;

	if( argc < 2 )
	{
		cerr << "Error: must provide one image name" << endl;
		cerr << "usage: " << argv[0] << " [input image] [optional: output image]" << endl;
		return 1;
	}



	try

	{
		// register the gdal driver
		GDALAllRegister();
		CPLSetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "FALSE");
		CPLSetConfigOption("CPL_DEBUG", "ON");
		//CPLSetConfigOption("GDAL_CACHEMAX", "256"); //large vals here make GDALCLose v slow - rather use GDALSetCacheMax64(GIntBig nNewSize) programmatically where it is needed i.e. for writing

		// Read the image
		//Mat image = cv::imread(argv[1], CV_LOAD_IMAGE_COLOR); //-1 for as is
		//Mat image = imread_geo(argv[1], CV_16UC1);
		vector<Mat> images = GdalImread(argv[1]);

		// Save the image
		if(argc > 2)
			imwrite(argv[2], images[0]);

		int ch = images[0].channels();
		vector<Mat> rgbV(3);
		rgbV[2] = images[3];
		rgbV[1] = images[0];
		rgbV[0] = images[1];
		Mat rgb;
		cv::merge(rgbV, rgb);

		namedWindow("IMAGE", cv::WINDOW_NORMAL|CV_WINDOW_KEEPRATIO|CV_GUI_EXPANDED);
		imshow("IMAGE", rgb*16);
		waitKey(60000);
	}
	catch( string e )
	{
		cout << e << endl;
	}

	CPLErrorReset();
	GDALDumpOpenDatasets(stderr);
	GDALDestroyDriverManager();

	return 0;
}

