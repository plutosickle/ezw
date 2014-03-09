/*
 *  wavelet.hpp
 *  ezw
 *
 *  Created by Tim Thirion on 4/10/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __wavelet_hpp__
#define __wavelet_hpp__

#include "matrix.hpp"

//------------------------------------------------------------------------------
// Haar
//------------------------------------------------------------------------------
template<typename Real>
void Haar(Matrix<Real> &H, const int size) {
	assert(H.Rows() == size);
	assert(H.Columns() == size);
	
	H.SetIdentity();
	Matrix<Real> A(size, size);
	Matrix<Real> R(size, size);
	
	int i(size);
	int half(0);
	while (i >= 2) {
		A.SetZero();
		half = i >> 1;
		for (int j=0; j < half; ++j) {
			// Low pass filter
			A(j, 2*j)	= Real(0.707106781186547);
			A(j, 2*j+1)	= Real(0.707106781186547);
			
			// High pass filter
			A(j+half, 2*j)	 = Real(0.707106781186547);
			A(j+half, 2*j+1) = Real(-0.707106781186547);
		}
		for (int j=i; j<size; ++j) {
			A(j,j) = Real(1);
		}
		Multiply(R, A, H);
		H.Set(R);
		i = half;
	}
}

#endif

