/*
 *  Scan.h
 *  ezw
 *
 *  Created by Tim Thirion on 4/10/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __scan_hpp__
#define __scan_hpp__

#include <vector>
#include "matrix.hpp"

typedef std::pair<int, int>	Index;

//------------------------------------------------------------------------------
// MortonScan
//------------------------------------------------------------------------------
void MortonScan(Matrix<int> &R, int size) {
	assert(size > 1);
	
	if (size == 2) {
		R.Resize(2,2);
		R(0,0) = 0;
		R(0,1) = 1;
		R(1,0) = 2;
		R(1,1) = 3;
	} else {
		int half = size >> 1;
		int halfSquared = half*half;
		Matrix<int> B(half, half);
		MortonScan(B, half);
		
		R.Resize(size, size);
		R.SetSubmatrix(0, 0, B);
		Add(B, B, halfSquared);
		R.SetSubmatrix(0, half, B);
		Add(B, B, halfSquared);
		R.SetSubmatrix(half, 0, B);
		Add(B, B, halfSquared);
		R.SetSubmatrix(half, half, B);
	}
}

//------------------------------------------------------------------------------
// VectorizeScan
//------------------------------------------------------------------------------
template<typename Real>
void VectorizeScan(std::vector<Index>& scan,
				   const Matrix<Real> &M)
{
	Index index(0,0);
	scan.reserve(M.Rows() * M.Columns());
	for (unsigned i=0; i<scan.capacity(); ++i) {
		scan.push_back(index);
	}
	for (int i=0; i<M.Rows(); ++i) {
		for (int j=0; j<M.Columns(); ++j) {
			scan[M(i,j)].first	= i;
			scan[M(i,j)].second	= j;
		}
	}
}

#endif

