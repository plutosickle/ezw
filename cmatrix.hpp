/*
 *  cmatrix.hpp
 *  ezw
 *
 *  Created by Tim Thirion on 9/28/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __cmatrix_hpp__
#define __cmatrix_hpp__

#include "matrix.hpp"

struct CMatrix
{
    Matrix<unsigned>    M;
    float               scale;
    float               bias;
    float               value;
    int                 qbits;
    int                 rows;
    int                 columns;
};

#endif
