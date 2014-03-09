/*
 *  Matrix.h
 *  ezw
 *
 *  Created by Tim Thirion on 4/7/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __matrix_hpp__
#define __matrix_hpp__

#include <iostream>
#include <iomanip>
#include <cassert>

template<typename Real>
class Matrix
{
public:
    Matrix()
    : R(0), C(0), array(0)
    {
    }

    Matrix(int rows, int columns=1)
    : R(rows), C(columns)
    {
        array = new Real[R*C];
        for (int i=0; i<R*C; ++i) {
            array[i] = Real(0);
        }
    }

    const Matrix<Real>& operator=(const Matrix<Real> &M) {
        assert(R == M.Rows());
        assert(C == M.Columns());

        for (int i=0; i<M.Rows()*M.Columns(); ++i) {
            array[i] = M[i];
        }
        return *this;
    }

    ~Matrix() {
        delete [] array;
    }

    Real& operator[](const int i) {
        assert(i < R*C);
        return array[i];
    }

    const Real& operator[](const int i) const {
        assert(i < R*C);
        return array[i];
    }

    Real& operator()(const int r, const int c) {
        assert(r >= 0 && r < R);
        assert(c >= 0 && c < C);
        return array[r*C+c];
    }

    const Real& operator()(const int r, const int c) const {
        assert(r >= 0 && r < R);
        assert(c >= 0 && c < C);
        return array[r*C+c];
    }

    int Rows() { return R; }
    const int Rows() const { return R; }
    int Columns() { return C; }
    const int Columns() const { return C; }

    void Resize(const int rows, const int columns) {
        delete [] array;
        R = rows;
        C = columns;
        array = new Real[R*C];
    }

    void Set(const Matrix<float> &M) {
        assert(R == M.Rows());
        assert(C == M.Columns());
        for (int i=0; i<M.Rows()*M.Columns(); ++i) {
            array[i] = M[i];
        }
    }

    void SetSubmatrix(int i, int j, const Matrix<Real> &M) {
        assert(i+M.Rows() <= R);
        assert(j+M.Columns() <= C);
        for (int r=0; r<M.Rows(); ++r) {
            for (int c=0; c<M.Columns(); ++c) {
                array[(r+i)*C+(c+j)] = M(r,c);
            }
        }
    }

    void SetZero() {
        for (int r=0; r<R; ++r) {
            for (int c=0; c<C; ++c) {
                array[r*C+c] = Real(0);
            }
        }
    }

    void SetIdentity() {
        assert(R == C);
        for (int r=0; r<R; ++r) {
            for (int c=0; c<C; ++c) {
                if (r == c) {
                    array[r*C+c] = Real(1);
                } else {
                    array[r*C+c] = Real(0);
                }
            }
        }

        for (int i=0; i<R; ++i) {
            array[i*C+i] = Real(1);
        }
    }

private:
    int R, C;
    Real *array;
};

template<typename Real>
void Scale(Matrix<Real> &R,
           const Real s,
           const Matrix<Real>& A)
{
    assert(R.Rows() == A.Rows());
    assert(R.Columns() == A.Columns());

    for (int i=0; i<A.Rows()*A.Columns(); ++i) {
        R[i] = s * A[i];
    }
}

template<typename Real>
void Add(Matrix<Real> &R,
         const Matrix<Real> &A,
         const Real s)
{
    assert(R.Rows() == A.Rows());
    assert(R.Columns() == R.Rows());

    for (int i=0; i<A.Rows()*A.Columns(); ++i) {
        R[i] = A[i] + s;
    }
}

template<typename Real>
void Add(Matrix<Real> &R,
         const Matrix<Real> &A,
         const Matrix<Real> &B)
{
    assert(R.Rows() == A.Rows() == B.Rows());
    assert(R.Columns() == A.Columns() == B.Columns());

    for (int i=0; i<A.Rows()*A.Columns(); ++i) {
        R[i] = A[i] + B[i];
    }
}

//------------------------------------------------------------------------------
// Subtract
//------------------------------------------------------------------------------
template<typename Real>
void Subtract(Matrix<Real> &R,
              const Matrix<Real> &A,
              const Matrix<Real> &B)
{
    assert(R.Rows() == A.Rows());
    assert(A.Rows() == B.Rows());
    assert(R.Columns() == A.Columns());
    assert(A.Columns() == B.Columns());

    for (int i=0; i<A.Rows()*A.Columns(); ++i) {
        R[i] = A[i] - B[i];
    }
}

//------------------------------------------------------------------------------
// Multiply
//------------------------------------------------------------------------------
template<typename Real>
void Multiply(Matrix<Real> &R,
              const Matrix<Real> &A,
              const Matrix<Real> &B)
{
    assert(A.Columns() == B.Rows());
    assert(R.Rows() == A.Rows());
    assert(R.Columns() == B.Columns());

    for (int i=0; i<A.Rows(); ++i) {
        for (int j=0; j<B.Columns(); ++j) {
            Real sum(0);
            for (int k=0; k < A.Columns(); ++k) {
                sum += A(i, k) * B(k, j);
            }
            R(i, j) = sum;
        }
    }
}

//------------------------------------------------------------------------------
// Transpose
//------------------------------------------------------------------------------
template<typename Real>
void Transpose(Matrix<Real> &R, const Matrix<Real> &A) {
    assert(R.Rows() == A.Rows());
    assert(R.Columns() == A.Columns());

    for (int r=0; r<A.Rows(); ++r) {
        for (int c=0; c<A.Columns(); ++c) {
            R(c, r) = A(r, c);
        }
    }
}

//------------------------------------------------------------------------------
// operator<<
//------------------------------------------------------------------------------
template<typename Real>
std::ostream& operator<<(std::ostream& os, const Matrix<Real>& M) {
    for (int r=0; r<M.Rows(); ++r) {
        for (int c=0; c<M.Columns(); ++c) {
            os << std::setw(8);
            os << M(r, c);
            if (c == M.Columns()-1) {
                os << "\n";
            } else {
                os << ", ";
            }
        }
    }
    return os;
}

void Cast(Matrix<float> &M, const Matrix<unsigned> &U) {
    M.Resize(U.Rows(), U.Columns());

    for (int r=0; r<U.Rows(); ++r) {
        for (int c=0; c<U.Columns(); ++c) {
            M(r,c)=float(U(r,c));
        }
    }
}

#endif

