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
    size_t m_n;
    size_t m_m;
    size_t m_o;

public:
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

    T& operator()(const size_t &i, const size_t &j, const size_t &k)
    {
		return *(this->m_buf + j + (i * this->m_m) + (k * this->m_m * this->m_n));
    }

	T*& Buf()
	{
		return m_buf;
	}
};

