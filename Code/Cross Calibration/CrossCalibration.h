/******************************************************************************
 * $Id: CrossCalibration.h $
 *
 * Project:  Cross calibration
 * Purpose:  Utility buffer classes
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
#pragma once
#include <opencv2/core/core.hpp>

template <typename T>
class A2D 
{
    T *m_buf;
    size_t m_n; 
    size_t m_m; 
public:
    A2D(T *buf, const size_t &n, const size_t &m)
        : m_buf(buf), m_n(n), m_m(m) { }
    ~A2D() { }

    T& operator()(const size_t &i, const size_t &j)
    {
        return *(this->m_buf + i * this->m_m + j);
    }
};

template <typename T>
class Buf3d 
{
private:
    T* m_buf;
    size_t m_n; //rows
    size_t m_m; //cols
    size_t m_o; //bands 

public:
	// rows, cols, bands
    Buf3d(const size_t &n, const size_t &m, const size_t &o)
        : m_n(n), m_m(m), m_o(o) 
	{ 
		this->m_buf = new T[n*m*o];
		std::memset(this->m_buf, 0, n*m*o*sizeof(T));
	}
    ~Buf3d() 
	{ 
		delete[] (this->m_buf);
	}

	// row, col, band
    T& operator()(const size_t &i, const size_t &j, const size_t &k)  
    {
		return *(this->m_buf + j + (i * this->m_m) + (k * this->m_m * this->m_n));
    }

	T*& Buf()
	{
		return m_buf;
	}

	// open cv Mat data ordering is channel (i.e. a pixel is stored as all channels in the image), col, row while gdal (geotiff) gives col, row, channel
	// this means we cannot simply wrap a multiband gdal raster with a opencv Mat header
	// so this function creates a vector of Mats that point to each band of the gdal raster
	// (newer versions of opencv support n-dimensional Mats - I would think that a 3 dim Mat should then follow the gdal ordering, but the 2.4.8 opencv does
	// not seem to support this properly or at all and I don't particularly want to upgrade wo good reason)
	std::vector<cv::Mat_<T>> ToMatVec()   
	{
		std::vector<cv::Mat_<T>> bandVec(this->m_o);
		for (int band = 0; band < this->m_o; band++)
		{
			bandVec[band] = this->ToMat(band);  //cv::Mat_<T>(this->m_n, this->m_m, this->m_buf + (band * this->m_m * this->m_n));
			// C++: Mat::Mat(int rows, int cols, int type, void* data, size_t step=AUTO_STEP)
			// crude switch to support templates - m
			//if (sizeof(T) == 8)
			//	bandVec[band] = cv::Mat(this->m_n, this->m_m, CV_64F, this->m_buf + (band * this->m_m * this->m_n));
			//else if (sizeof(T) == 4)
			//	bandVec[band] = cv::Mat(this->m_n, this->m_m, CV_32F, this->m_buf + (band * this->m_m * this->m_n));
			//else if (sizeof(T) == 2)
			//	bandVec[band] = cv::Mat(this->m_n, this->m_m, CV_16U, this->m_buf + (band * this->m_m * this->m_n));
			//else if (sizeof(T) == 1)
			//	bandVec[band] = cv::Mat(this->m_n, this->m_m, CV_8U, this->m_buf + (band * this->m_m * this->m_n));
			//else
			//	throw new std::exception("Unsupported data type for conversion cv::Mat");
		}
		return bandVec;
	}

	cv::Mat_<T> ToMat(int band)
	{
		//Mat_(int _rows, int _cols, _Tp* _data, size_t _step = AUTO_STEP);
		cv::Mat_<T> bandMat = cv::Mat_<T>(this->m_n, this->m_m, this->m_buf + (band * this->m_m * this->m_n));
		return bandMat;
	}
};

enum ModelForms
{
	GAIN_ONLY = 1,
	GAIN_AND_OFFSET = 2,
	OFFSET_ONLY = 3,
	GAIN_AND_IMAGE_OFFSET = 4,
	IMAGE_GAIN_AND_OFFSET = 5,
	//GAIN_AND_OFFSET_IND_WIN = 6
};