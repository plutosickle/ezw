/*
 *  math.hpp
 *  ezw
 *
 *  Created by Tim Thirion on 4/7/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __math_hpp__
#define __math_hpp__

#define INFINITY 1e20

template<typename Real>
Real Abs(const Real x) {
	return x > Real(0) ? x : -x;
}

bool IsPowerOfTwo(const int n) {
	return !(n & (n-1));
}

#endif
