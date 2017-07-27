#pragma once
#include <string>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/ml/ml.hpp>

#include <cpl_conv.h>
#include <gdal_priv.h>
#include "EntropyFilt.h"

using namespace cv;
using namespace std;

enum ClassifierType
{
	TypeRandomTrees = 0,
	TypeKNearest = 1,
	TypeSVM = 2,
	TypeNormalBayes = 3,
	TypeDecisionTree = 4
};

class SpekboomClassifier
{
public:
	SpekboomClassifier(ClassifierType clfrType);
	~SpekboomClassifier();
	void Configure(const string& configFileName);
	void Classify(const string& inputFileName, const string& outputFileName);
	vector<Mat> ExtractFeatures(const Mat& image);
	void Train(const string& dataFileName);

	static vector<Mat> GdalImread(const string& fileName);
	static int CvDataType(GDALDataset* ds);
	static int CvDataType(GDALRasterBand* band);
	static void ShowImageBands(const vector<Mat>& bands, const string& winTitle, const bool doWait);
	static void ShowImage(const Mat& image, const string& winTitle, const bool doWait);
	static Mat StretchImage(const Mat& bands);

private:
	//string InputFileName;
	//vector<Mat> bands;
	//string OutputFileName;
	Mat rgbPca, rggPca, scale;
	SVM svm;
	RandomTrees randomForest;
	KNearest kNearest;
	NormalBayesClassifier normalBayes;
	DecisionTree decisionTree;
	EntropyFilt* entropyFilt;	
	cv::Size blockSize;
	ClassifierType clfrType;
};

