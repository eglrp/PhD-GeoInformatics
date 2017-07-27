#include "SpekBoomClassifier.h"

// OpenCV Dependencies
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>
#include <ogr_api.h>
#include <ogr_geometry.h>

#define LOG 0
#define MAX_PATH 1024
#define DO_OUT_BOUND 0

SpekboomClassifier::SpekboomClassifier(ClassifierType _clfrType = ClassifierType::TypeRandomTrees)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "FALSE");
	CPLSetConfigOption("CPL_DEBUG", "ON");
	OGRRegisterAll();
#if LOG
	this->blockSize = cv::Size(128, 128);
#else
	this->blockSize = cv::Size(128, 128);
	//this->blockSize = cv::Size(512, 512);
#endif

	this->entropyFilt = new EntropyFilt(cv::Size(5, 5), this->blockSize);
	this->clfrType = _clfrType;

	//CPLSetConfigOption("GDAL_CACHEMAX", "256"); //large vals here make GDALCLose v slow - rather use GDALSetCacheMax64(GIntBig nNewSize) programmatically where it is needed i.e. for writing
}

SpekboomClassifier::~SpekboomClassifier()
{
	delete this->entropyFilt;
	OGRCleanupAll();
	CPLErrorReset();
	GDALDumpOpenDatasets(stderr);
	GDALDestroyDriverManager();
	destroyAllWindows();
}

void SpekboomClassifier::Configure(const string& configFileName)
{
	//clfrType = ClassifierType::TypeRandomTrees;

	cv::FileStorage fs;

	fs.open(configFileName, FileStorage::READ);
	//this->svm.read(fs, "svm");
	//fs["svm"] >> this->svm;
	//this->svm.lo
	fs["rgbPca"] >> this->rgbPca;
	fs["rggPca"] >> this->rggPca;
	fs["scale"] >> this->scale;

	//train knn
	Mat scaleRot = Mat::diag(Mat(this->scale, cv::Rect(0, 0, this->scale.cols, 1)));
	Mat scaleOffset = Mat(this->scale, cv::Rect(0, 1, this->scale.cols, 1));

	/*Mat tr, trLab;
	fs["tr"] >> tr;
	fs["trLab"] >> trLab;*/

	//tr = tr*scaleRot + cv::repeat(scaleOffset, tr.rows, 1);
	//bool res = this->kNearest.train(tr, trLab);


	//CPLGetFilename(configFileName)
	//CPLGetPath(configFileName);
	string svmConfigFileName = (configFileName.substr(0, configFileName.length() - 5) + "Svm.yaml");
	this->svm.load(svmConfigFileName.c_str(), "my_svm");
	this->randomForest.load((configFileName.substr(0, configFileName.length() - 5) + "RTreesSmall.yaml").c_str(), "my_random_trees");
	this->normalBayes.load((configFileName.substr(0, configFileName.length() - 5) + "NormalBayes.yaml").c_str(), "my_nb");
	this->decisionTree.load((configFileName.substr(0, configFileName.length() - 5) + "DTree.yaml").c_str(), "my_tree");

#if LOG
	cout << "rgbPca: " << this->rgbPca << endl << endl;
	cout << "rggPca: " << this->rggPca << endl << endl;
	cout << "scale: " << this->scale << endl << endl;
#endif
	//fs.
}

cv::Mat ConfusionMatrix(cv::Mat& responses, cv::Mat& labels)
{
	double minVal, maxVal;
	//assume labels start at 1 and end at max
	cv::minMaxLoc(labels, &minVal, &maxVal);
	int nc = maxVal;
	cv::Mat res(nc, nc, CV_32F);
	res = 0;

	for (int i = 0; i < responses.rows; i++)
	{
		res.at<float>(labels.at<float>(i, 0)-1, responses.at<float>(i, 0)-1)++;
	}
	cout << "Confusion Matrix: " << endl << res << endl;
	cv::Mat classSizes;
	cv::reduce(res, classSizes, 1, CV_REDUCE_SUM);
	res = res/cv::repeat(classSizes, 1, nc);
	cout << "Normalised Confusion Matrix: " << endl << res << endl;
	cout << "Prior independent error: " << endl << 1.0 - cv::sum(res.diag(0)).val[0]/nc << endl << endl;
	return res;
}

void SpekboomClassifier::Train(const string& dataFileName)
{
	cv::FileStorage fs;

	fs.open(dataFileName, FileStorage::READ);

	//train knn
	Mat scaleRot = Mat::diag(Mat(this->scale, cv::Rect(0, 0, this->scale.cols, 1)));
	Mat scaleOffset = Mat(this->scale, cv::Rect(0, 1, this->scale.cols, 1));

	Mat tr, trLab;
	fs["tr"] >> tr;
	fs["trLab"] >> trLab;
	fs.release();

	bool res = false;

	switch (this->clfrType)
	{
	case ClassifierType::TypeRandomTrees:
		{
			//flatFeatIm = flatFeatIm*scaleRot + cv::repeat(scaleOffset, flatFeatIm.rows, 1);
			CvRTParams rtParams;
			float priors[3] = {0.2,0.6,0.2}; //{1/3, 1/3, 1/3};
			//rtParams = CvRTParams(7, (int)0.01 * tr.rows, 0.05, true, 3, (const float*)&priors, true, 4, 50, 0.05, CV_TERMCRIT_ITER); 
			rtParams = CvRTParams(10, (int)0.01 * tr.rows, 0.025, false, 3, (const float*)&priors, false, 4, 10, 0.025, CV_TERMCRIT_ITER); 
			/*rtParams.max_depth = 7;
			rtParams.nactive_vars = 4;
			rtParams.calc_var_importance = false;
			rtParams.regression_accuracy = 0.05;*/

			this->randomForest.train(tr, CV_ROW_SAMPLE, trLab, cv::Mat(), cv::Mat(), cv::Mat(), cv::Mat(), rtParams);
			cv::Mat responses(tr.rows, 1, CV_32F);
			float* outPtr = responses.ptr<float>(0);
			for (int i = 0; i < tr.rows; i++)
				outPtr[i] = this->randomForest.predict(tr.row(i));
			
			ConfusionMatrix(responses, trLab);
		}
		//float error = this->randomForest.calc_error(tr, CV_TRAIN_ERROR, trLab);
		break;
	case ClassifierType::TypeSVM:
		{
			tr = tr*scaleRot + cv::repeat(scaleOffset, tr.rows, 1);

			cv::SVMParams svmParams;
			svmParams.svm_type = CvSVM::C_SVC;
			svmParams.kernel_type = CvSVM::RBF;
			svmParams.gamma = 5;
			svmParams.C = 1;
			double cw[3] = {1, 3, 1};
			cv::Mat cwm(1, 3, CV_64F, cw);
			CvMat cwm_ = cwm;
			svmParams.class_weights = &cwm_;
			res = this->svm.train(tr, trLab, cv::Mat(), cv::Mat(), svmParams);

			cv::Mat responses(tr.rows, 1, CV_32F);
			this->svm.predict(tr, responses);
			ConfusionMatrix(responses, trLab);
		}
		break;
	case ClassifierType::TypeKNearest:
		{
			tr = tr*scaleRot + cv::repeat(scaleOffset, tr.rows, 1);
			res = this->kNearest.train(tr, trLab);

			cv::Mat responses(tr.rows, 1, CV_32F);
			//this->kNearest.
			Mat neighbourResponses, dists;
			this->kNearest.find_nearest(tr, 5, responses, neighbourResponses, dists);
			ConfusionMatrix(responses, trLab);
		}
	case ClassifierType::TypeNormalBayes:
		{
			res = this->normalBayes.train(tr, trLab);

			cv::Mat responses(tr.rows, 1, CV_32F);
			this->normalBayes.predict(tr, &responses);
			ConfusionMatrix(responses, trLab);
		}
	case ClassifierType::TypeDecisionTree:
		{
/*[err, cerr, nlabOut] = prcrossval(subData(:, feats), opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 5, 1);*/

			float priors[3] = {0.25, 0.5, 0.25};
			CvDTreeParams dtreeParams;
			dtreeParams.max_depth = 12;
			dtreeParams.use_1se_rule = false;
			dtreeParams.use_surrogates = false;
			dtreeParams.cv_folds = 5;
			dtreeParams.truncate_pruned_tree = true;
			dtreeParams.priors = priors;
			dtreeParams.min_sample_count = 20;

			res = this->decisionTree.train(tr, CV_ROW_SAMPLE, trLab, cv::Mat(), cv::Mat(), cv::Mat(), cv::Mat(), dtreeParams);

			cv::Mat responses(tr.rows, 1, CV_32F);
			float* outPtr = responses.ptr<float>(0);
			for (int i = 0; i < tr.rows; i++)
				outPtr[i] = this->decisionTree.predict(tr.row(i))->value;
			ConfusionMatrix(responses, trLab);
			//dtree->predict(sample,mask)->value
		}
		break;
	}
}

void SpekboomClassifier::ShowImageBands(const vector<Mat>& bands, const string& winTitle = "ShowImageBands", const bool doWait = true)
{
	namedWindow(winTitle, cv::WINDOW_NORMAL|CV_WINDOW_KEEPRATIO|CV_GUI_EXPANDED);
	for (int i = 0; i < bands.size(); i++)
	{
		double min, max;
		cv::minMaxLoc(bands[i], &min, &max);		

		imshow(winTitle, SpekboomClassifier::StretchImage(bands[i]));
		if (doWait)
			waitKey(0);
	}
	if (doWait)
		destroyWindow(winTitle);
}

Mat SpekboomClassifier::StretchImage(const Mat& image)
{
	vector<Mat> bands;
	cv::split(image, bands);
	vector<Mat> stretchBands(bands.size());
	double min, max;
	for (int i = 0; i < bands.size(); i++)
	{
		cv::minMaxLoc(bands[i], &min, &max);
		stretchBands[i] = (bands[i]-min)/(max-min);
	}
	Mat res;
	cv::merge(stretchBands, res);
	return res;
}

void SpekboomClassifier::ShowImage(const Mat& image, const string& winTitle = "ShowImageBands", const bool doWait = true)
{
	vector<Mat> bands;
	if (image.channels() == 3) //display as rgb image rather than splitting into bands
	{
		bands.resize(1);
		bands[0] = image;
	}
	else
		cv::split(image, bands);

	ShowImageBands(bands, winTitle, doWait);
}


void SpekboomClassifier::Classify(const string& inputFileName, const string& outputFileName)
{
	vector<Mat> bands = SpekboomClassifier::GdalImread(inputFileName);
	Mat image;
	cv::merge(bands, image);

	//cv::Size blockSize(128, 128);
	//cv::Size numBlocks(image.cols / blockSize.width, image.rows / blockSize.height);
#if LOG
	cv::FileStorage fs;
	fs.open("D:/results.yaml.gz", cv::FileStorage::WRITE);
#endif

    char **papszOptions = NULL;
	papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "LZW");
	papszOptions = CSLSetNameValue( papszOptions, "TILED", "YES");

	int cacheVal = GDALGetCacheMax();
	GDALSetCacheMax(1024*1024*1024); //for fast reading and GDALCLose (in bytes)

	GDALDataset* inputDataSet = (GDALDataset*) GDALOpen(inputFileName.c_str(), GA_ReadOnly);
	GDALDataset* outputDataSet = inputDataSet->GetDriver()->Create(outputFileName.c_str(), inputDataSet->GetRasterXSize(), 
		inputDataSet->GetRasterYSize(), 1, GDALDataType::GDT_Byte, papszOptions);

	if (outputDataSet == NULL)
		throw string("Error: could not create: ") + outputFileName;

	double inputGeoTransform[6];
	inputDataSet->GetGeoTransform(inputGeoTransform);
	outputDataSet->SetGeoTransform(inputGeoTransform);
	outputDataSet->SetProjection(inputDataSet->GetProjectionRef());
	//res = calibDataSet->SetMetadataItem("BitsPerSample", "12", NULL);
	GDALRasterBand* outputBand = outputDataSet->GetRasterBand(1);
	int res = outputBand->SetNoDataValue(0);

#if DO_OUT_BOUND
	string outputBoundFileName = outputFileName.substr(0, outputFileName.length()-4) + "_BND.tif";

	GDALDataset* outputBoundDataSet = inputDataSet->GetDriver()->Create(outputBoundFileName.c_str(), inputDataSet->GetRasterXSize(), 
		inputDataSet->GetRasterYSize(), 1, GDALDataType::GDT_Byte, papszOptions);
	if (outputBoundDataSet == NULL)
		throw string("Error: could not create: ") + outputBoundFileName;

	outputBoundDataSet->SetGeoTransform(inputGeoTransform);
	outputBoundDataSet->SetProjection(inputDataSet->GetProjectionRef());
	GDALRasterBand* outputClosedBand = outputBoundDataSet->GetRasterBand(1);
	res = outputClosedBand->SetNoDataValue(0);
#endif
	GDALClose(inputDataSet);
	CSLDestroy(papszOptions);

	Mat scaleRot = Mat::diag(Mat(this->scale, cv::Rect(0, 0, this->scale.cols, 1)));
	Mat scaleOffset = Mat(this->scale, cv::Rect(0, 1, this->scale.cols, 1));

	bool log = false;
	int blockTtl = (1 + image.cols * image.rows)/(blockSize.width * blockSize.height);
	int blockProgressUpdate = 0;
	int blockCount = 0;
	int64 startTick = 0, processTicks[10]={};

	Mat strEl = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(-1, -1));
	Mat strEl5 = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(-1, -1));

#if DO_OUT_BOUND
	//Mat strElBig = getStructuringElement(MORPH_RECT, Size(32, 32), Point(-1, -1));
#endif

#if LOG
	for (int x = this->blockSize.width*2*4; x < image.cols; x += blockSize.width)
#else
	//namedWindow("Classification");
	//for (int x = this->blockSize.width*2*4; x < image.cols; x += blockSize.width)	
	for (int x = 0; x < image.cols; x += blockSize.width)
#endif
	{
		//TO DO: change back to 0
#if LOG
		for (int y = this->blockSize.height*2*4; y < image.rows; y += blockSize.height)
#else
		//for (int y = this->blockSize.height*2*4; y < image.rows; y += blockSize.height)
		for (int y = 0; y < image.rows; y += blockSize.height)
#endif
		{
			//TO DO: partial block for edges
			//TO DO: convert types to single/double
				// Compute the portion of the block that is valid
				// for partial edge blocks.
			int nXValid, nYValid;
			Size validBlockSize;
			int iXBlock = x/blockSize.width, iYBlock = y/blockSize.height;

			if ((iXBlock+1) * blockSize.width > image.cols)
				validBlockSize.width = image.cols - iXBlock * blockSize.width;
			else
				validBlockSize.width = blockSize.width;

			if ((iYBlock+1) * blockSize.height > image.rows)
				validBlockSize.height = image.rows - iYBlock * blockSize.height;
			else
				validBlockSize.height = blockSize.height;
			startTick = cv::getTickCount();
			Mat block(image, cv::Rect(x, y, validBlockSize.width, validBlockSize.height));
			vector<Mat> noDataChannels;
			cv::split(block==0, noDataChannels);
			Mat noDataMask; //(block.rows, block.cols, CV_8U), noDataMaskB;
			noDataMask = noDataChannels[0];
			for (int c = 1; c < noDataChannels.size(); c++)
				noDataMask = noDataMask & noDataChannels[c];
			//cv::compare(block, cv::Scalar(0), noDataMask, CV_CMP_EQ); // (block == 0); //cv::Scalar(0,0,0,0);
			//cv::reduce(block, noDataMaskB, 2, CV_RE

			//noDataMask.convertTo(noDataMaskB, CV_8U);
			//int type = noDataMaskB.type();
			//type = noDataMaskB.channels();

			Mat blockSingle;
			block.convertTo(blockSingle, CV_32F);

			vector<Mat> feats = this->ExtractFeatures(blockSingle/(4096*4));
			Mat featIm;
			cv::merge(feats, featIm);
			
			//NB scalem([], 'variance')
			Mat flatFeatIm = featIm.reshape(1, block.rows*block.cols);

			Mat flatOutIm(1, flatFeatIm.rows, CV_32FC1);
			float* outPtr = flatOutIm.ptr<float>(0);
			processTicks[0] += (cv::getTickCount() - startTick);
			startTick = cv::getTickCount();
			switch (this->clfrType)
			{
			case ClassifierType::TypeSVM:
				flatFeatIm = flatFeatIm*scaleRot + cv::repeat(scaleOffset, flatFeatIm.rows, 1);
				this->svm.predict(flatFeatIm, flatOutIm);
				break;

			case ClassifierType::TypeRandomTrees:
				for (int i = 0; i < flatFeatIm.rows; i++)
					outPtr[i] = this->randomForest.predict(flatFeatIm.row(i));
				break;
				
			case ClassifierType::TypeNormalBayes:
				//flatFeatIm = flatFeatIm*scaleRot + cv::repeat(scaleOffset, flatFeatIm.rows, 1);
				this->normalBayes.predict(flatFeatIm, &flatOutIm);
				break;
			case ClassifierType::TypeDecisionTree:
				for (int i = 0; i < flatFeatIm.rows; i++)
					outPtr[i] = this->decisionTree.predict(flatFeatIm.row(i))->value;
				break;
			case ClassifierType::TypeKNearest:
				{
					Mat neighbourResponses, dists;
					flatFeatIm = flatFeatIm*scaleRot + cv::repeat(scaleOffset, flatFeatIm.rows, 1);
					this->kNearest.find_nearest(flatFeatIm, 5, flatOutIm, neighbourResponses, dists);
				}
				break;
			default:
				throw(string("Error: unknown classifier type"));
				break;
			}
#if LOG
			cout << ".";
#endif
			Mat outIm = (flatOutIm.reshape(1, block.rows)); //outIm must be continuous
			outIm.setTo(0, noDataMask);
			
			//swap sb and tree
			Mat sbMask = outIm==2;
			outIm = 0;
			outIm.setTo(255, sbMask);
			//outIm.setTo(2, outIm==3);
			//outIm.setTo(3, sbMask);

			//default border behaviour should be to not apply morphology there
			//ShowImage(outIm, "outIm", false);
			cv::morphologyEx(outIm, outIm, CV_MOP_OPEN, strEl);
			cv::morphologyEx(outIm, outIm, CV_MOP_CLOSE, strEl);
			//ShowImage(outIm, "outIm Morph", true);
#if DO_OUT_BOUND
			//ShowImage(outIm, "outIm", false);
			Mat outBoundIm; // = outIm.clone();
			Mat outDilateIm;
			cv::morphologyEx(outIm, outDilateIm, CV_MOP_DILATE, strEl5);
			outBoundIm =  cv::abs(outDilateIm - outIm);
			//ShowImage(outBoundIm, "outBoundIm", true);

			CPLErr e = outputClosedBand->RasterIO(GDALRWFlag::GF_Write, x, y, validBlockSize.width, validBlockSize.height, (float*)outBoundIm.data, 
				validBlockSize.width, validBlockSize.height, GDALDataType::GDT_Float32, 0, 0); //sizeof(float), sizeof(float)*blockSize.width);
			if (e != 0)
				throw string("Error: outputClosedBand->RasterIO");

#endif
			//ShowImage(outIm, "outIm Morph", true);
			//ShowImage(noDataMask/255, "noDataMask", true);
			//outIm.setTo(0, noDataMaskB); TO DO combine 4 ch to 1
			//outIm.
			processTicks[1] += (cv::getTickCount() - startTick);
			startTick = cv::getTickCount();
			//we can't write in blocks as they usuallu consist of 1 row and we need at least winSize rows for entropyfilt
			CPLErr cplErr = outputBand->RasterIO(GDALRWFlag::GF_Write, x, y, validBlockSize.width, validBlockSize.height, (float*)outIm.data, 
				validBlockSize.width, validBlockSize.height, GDALDataType::GDT_Float32, 0, 0); //sizeof(float), sizeof(float)*blockSize.width);
			//outputBand->FlushCache();
			processTicks[2] += (cv::getTickCount() - startTick);

			if (cplErr != 0)
				throw string("Error: outputBand->RasterIO");

			if ((100*++blockCount)/blockTtl > blockProgressUpdate)
			{
				blockProgressUpdate += 10;
#if LOG
				cout << endl << blockProgressUpdate << "%, ";

				cout << "Process Times" << endl;
				cout << " ExtractFeatures: " << (1000*processTicks[0]/blockCount)/cv::getTickFrequency() << "ms" << endl;
				cout << " Classify: " << (1000*processTicks[1]/blockCount)/cv::getTickFrequency() << "ms" << endl;
				cout << " Write: " << (1000*processTicks[2]/blockCount)/cv::getTickFrequency() << "ms" << endl;
#else
				cout << blockProgressUpdate << "%, ";
#endif
				/*processTicks[0] = 0;
				processTicks[1] = 0;
				processTicks[2] = 0;*/
			}


#if LOG
			//TO DO: recheck feats against matlab - some look suspect
			//fs << "entropy" << feats[4];
			//ShowImageBands(feats, "Feats");

			vector<Mat> blockBands;
			cv::split(blockSingle, blockBands);
			vector<Mat> cirV(3);
			cirV[2] = blockBands[3];
			cirV[1] = blockBands[0];
			cirV[0] = blockBands[1];
			cv::Mat cir;
			cv::merge(cirV, cir);
			ShowImage(cir*16, "CIR", false);
			ShowImage(outIm/3, "Classification", true);

			cout << "scaleRot: " << scaleRot << endl << endl;
			cout << "scaleOffset: " << scaleOffset << endl << endl;

			if (blockCount > 1)
			{
				fs.release();
 				break;
			}
#else
			//cv::imshow("Classification", outIm/3);
#endif
		}
#if LOG
		break;
#endif
	}
	cout << endl << endl;
#if LOG
	fs.release();
#else
	//cv::destroyWindow("Classification");
#endif
	outputDataSet->FlushCache();
	GDALSetCacheMax(cacheVal); //for fast reading and GDALCLose (in bytes)
	
		/*OGRSFDriver *poDriver  = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
		OGRDataSource* ds = poDriver->CreateDataSource("d:/test.shp", NULL);
		OGRLayer* poLayer = ds->CreateLayer("Spekboom", NULL, wkbMultiPolygon, NULL);
		OGRFieldDefn oField("DF", OFTInteger);
		OGRErr oer = poLayer->CreateField(&oField);
		CPLErr er = GDALPolygonize(outputBand, NULL, poLayer, 0, NULL, NULL, NULL);
		poLayer->SyncToDisk();
		OGRDataSource::DestroyDataSource(ds);*/

	//GDALPolygonize(outputDataSet, NULL, NULL, 2, NULL, NULL);
	//OGRGeometry geom;
	//OGRPolygon poly = 
	//OGRGeometryH
	GDALClose(outputDataSet);
#if DO_OUT_BOUND
	outputBoundDataSet->FlushCache();
	GDALClose(outputBoundDataSet);
#endif
}

vector<Mat> SpekboomClassifier::ExtractFeatures(const Mat& image)
{
	vector<Mat> bands;
	cv::split(image, bands); //does this copy the data?

	Mat flatRgbImage = image.reshape(1, image.rows*image.cols);
	Mat sumChannels;
	cv::reduce(flatRgbImage, sumChannels, 1, CV_REDUCE_SUM, -1);
	Mat den = cv::repeat(sumChannels, 1, image.channels());
	Mat flatRggImage = flatRgbImage/den;
	int imageDims[] = {image.rows, image.cols};
	Mat rgG = flatRggImage.reshape(image.channels(), image.rows);
	vector<Mat> rggBands;
	cv::split(rgG, rggBands);
	//ShowImageBands(rggBands);

	Mat ndvi = (bands[3] - bands[0])/(bands[0] + bands[3]);

	Mat rgbPcaRot(this->rgbPca, cv::Rect(0, 0, image.channels(), image.channels()));
	Mat rgbPcaOffset(this->rgbPca, cv::Rect(0, image.channels(), image.channels(), 1));
	Mat rggPcaRot(this->rggPca, cv::Rect(0, 0, image.channels(), image.channels()));
	Mat rggPcaOffset(this->rggPca, cv::Rect(0, image.channels(), image.channels(), 1));

#if LOG
	/*cout << "rgbPcaRot: " << rgbPcaRot << endl << endl;
	cout << "rgbPcaOffset: " << rgbPcaOffset << endl << endl;
	cout << "rggPcaRot: " << rggPcaRot << endl << endl;
	cout << "rggPcaOffset: " << rggPcaOffset << endl << endl;*/
	/*Mat tst(flatRgbImage, cv::Rect(0,0,1,flatRgbImage.rows));
	Mat tst_ = tst.clone().reshape(1, image.rows);
	ShowImage(tst_, "TEST");*/
#endif

	Mat rgbPcaProj = flatRgbImage*rgbPcaRot + cv::repeat(rgbPcaOffset, flatRgbImage.rows, 1);
	Mat rggPcaProj = flatRggImage*rggPcaRot + cv::repeat(rggPcaOffset, flatRggImage.rows, 1);

 	cv::reduce(flatRgbImage, sumChannels, 1, CV_REDUCE_AVG, -1);
	sumChannels = sumChannels.reshape(1, image.rows);

	rgbPcaProj = rgbPcaProj.reshape(image.channels(), image.rows);
	rggPcaProj = rggPcaProj.reshape(image.channels(), image.rows);
	vector<Mat> rgbPcaBands;
	cv::split(rgbPcaProj, rgbPcaBands);
	vector<Mat> rggPcaBands;
	cv::split(rggPcaProj, rggPcaBands);

	Mat pc1, entropy(image.rows, image.cols, CV_32FC1);
	rgbPcaBands[0].convertTo(pc1, CV_8UC1, 256);
	if (image.rows == this->blockSize.height && image.cols == this->blockSize.width)
		entropy = this->entropyFilt->Execute(pc1);
	else
		entropy = EntropyFiltFn(pc1); //partial block

#if LOG
	//ShowImage(sumChannels, "sumChannels");
	//ShowImage(rgbPcaProj, "rgbPcaProj");
	//ShowImage(rggPcaProj, "rggPcaProj");
	//ShowImage(entropy, "entropy", false);
#endif

	vector<Mat> feats(6);
	feats[0] = ndvi;
	feats[1] = rgbPcaBands[0];
	feats[2] = rggPcaBands[1];
	feats[3] = entropy;
	feats[4] = rggBands[2];
	feats[5] = rggBands[1];

	return feats;
}


vector<Mat> SpekboomClassifier::GdalImread(const string& filename)
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
#if LOG
	cout << "Reading " << filename;
#endif
	for (int iBand = 0; iBand < nBands; iBand++) // < nBands; iBand++)
	{
#if LOG
		cout << ", Band " << iBand;
#endif
		int nXBlocks, nYBlocks, nXBlockSize, nYBlockSize;

		GDALRasterBand* band = dataset->GetRasterBand(iBand+1);
		//int nXSize = band->GetXSize(), nYSize = band->GetYSize();
		band->GetBlockSize(&nXBlockSize, &nYBlockSize);

		if (nYBlockSize > 1 && nXBlockSize < imgSize.width)
			throw string( "Error: blocks will not be continous");

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
				else //make roi header for each block
				{
					int nXValid, nYValid;
					Size validBlockSize;
					if (((iXBlock+1) * nXBlockSize > res[iBand].cols) || ((iYBlock+1) * nYBlockSize > res[iBand].rows)) //can't use ReadBlock for partial blocks
					{
						validBlockSize.width = res[iBand].cols - iXBlock * nXBlockSize;
						validBlockSize.height = res[iBand].rows - iYBlock * nYBlockSize;
						Mat roi(res[iBand], cv::Rect(iXBlock*nXBlockSize, iYBlock*nYBlockSize,  validBlockSize.width, validBlockSize.height));
						CPLErr err = band->RasterIO(GDALRWFlag::GF_Read, iXBlock * nXBlockSize, iYBlock * nYBlockSize, validBlockSize.width, validBlockSize.height, roi.data, 
							validBlockSize.width, validBlockSize.height, gdalDataType, 0, 0);
					}
					else
					{
						validBlockSize.width = nXBlockSize;
						validBlockSize.height = nYBlockSize;
						Mat roi(res[iBand], cv::Rect(iXBlock*nXBlockSize, iYBlock*nYBlockSize,  validBlockSize.width, validBlockSize.height));
						CPLErr err = band->ReadBlock(iXBlock, iYBlock, roi.data);
						if (err != CPLErr::CE_None)
						{
							throw string("Error: RasterIO failure");
						}
					}
				}
			}
		}
	}
	GDALSetCacheMax(cacheVal);
#if LOG
	cout << endl;
#endif

	GDALClose(dataset);
	return res;
}


int SpekboomClassifier::CvDataType(GDALDataset* ds)
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

int SpekboomClassifier::CvDataType(GDALRasterBand* band)
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

/**
* Test Program
*/
int main( int argc, char* argv[] )
{
	if( argc < 3 )
	{
		cerr << "Error: must provide two image names" << endl;
		cerr << "usage: " << argv[0] << " [input image] [output image] [-o] [-t 1-4]" << endl;
		cerr << "       -o: overwrite output image if it exists" << endl;
		cerr << "       -t: 0 = Rand Forest, 1 = KNN, 2 = SVM, 3 = Normal Bayes, 4 = Decision Tree  (default)" << endl;
		return 1;
	}

	bool doOverwrite = false;
	ClassifierType clfrType = ClassifierType::TypeDecisionTree;

	for (int i=3; i<argc; i++)
	{
		if (string(argv[i]).compare("-o") == 0)
		{
			doOverwrite = true;
			cout << "Will overwrite " << argv[2] << endl;
		}
		else if (string(argv[i]).compare("-t") == 0)
		{
			i++;
			clfrType = (ClassifierType)atoi(argv[i]);
			cout << "Using Classifier: " << clfrType << endl;
		}
	}

	if (!doOverwrite)
	{
		if (CPLCheckForFile(argv[2], NULL) == 1)
		{
			cout << string(argv[2]) << " exists.  Exiting" << endl;
			return 1;
		}
	}

	try
	{
		/*int data[] = {1, 2, 3, 4, 5, 6};
		cv::Mat x = cv::Mat(2,3,CV_32SC1, data);
		cout << "x: " << x << endl << endl;
		cout << "x.reshape(1,6): " << x.reshape(1,6) << endl << endl;
		cout << "x.reshape(1,1): " << x.reshape(1,1) << endl << endl;
		Mat x_ = x.t();
		cout << "x.t().reshape(1,6): " << x_.reshape(1,1) << endl << endl;
		return 0;*/
		SpekboomClassifier clfr(clfrType);
		//clfr.Configure("G:/MSc GeoInformatics/Data/NGI/My Rectified/Ground Truth Images/config.yaml");
		//clfr.Train("G:/MSc GeoInformatics/Data/NGI/My Rectified/Ground Truth Images/trDataDTree.yaml");
		char path[MAX_PATH];
		CPLGetExecPath(path, MAX_PATH);
		
		clfr.Configure(string(CPLGetDirname(path)) + "/config.yaml");
		clfr.Train(string(CPLGetDirname(path)) + "/trDataDTree.yaml");
		

		// Read the image
		//Mat image = cv::imread(argv[1], CV_LOAD_IMAGE_COLOR); //-1 for as is
		//Mat image = imread_geo(argv[1], CV_16UC1);
		clfr.Classify(argv[1], argv[2]);
		return 0;

		vector<Mat> images = SpekboomClassifier::GdalImread(argv[1]);

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
	catch (string e)
	{
		cout << e << endl;
	}


	return 0;
}

