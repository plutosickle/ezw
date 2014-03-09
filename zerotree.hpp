/*
 *  Zerotree.h
 *  ezw
 *
 *  Created by Tim Thirion on 4/10/06.
 *  Copyright 2006 Tim Thirion. All rights reserved.
 *
 */
#ifndef __Zerotree_h__
#define __Zerotree_h__

#include <fstream>
#include "scan.hpp"
#include "math.hpp"
#include "wavelet.hpp"
#include "cmatrix.hpp"

//------------------------------------------------------------------------------
// HasSignificantChildren
//------------------------------------------------------------------------------
template<typename Real>
bool HasSignificantChildren(const Matrix<Real> &W,
                            const std::vector<Index> &scan,
                            const int index,
                            const Real threshold)
{
    assert(W.Rows() == W.Columns());

    int size = scan.size();
    bool result = false;
    for (int i=4*index; i<4*(index+1); ++i) {
        if (i<size) {
            int r = scan[i].first;
            int c = scan[i].second;
            if (Abs(W(r, c)) >= threshold) {
                return true;
            } else {
                result |= HasSignificantChildren(W, scan, i, threshold);
            }
        }
    }
    return result;
}

//------------------------------------------------------------------------------
// SetZerotree
//------------------------------------------------------------------------------
void SetZerotreeChildren(Matrix<int> &R,
                         const std::vector<Index> &scan,
                         const int index)
{
    assert(index > 0);
    assert(R.Rows() == R.Columns());

    int size = scan.size();
    for (int i=4*index; i<4*(index+1); ++i) {
        if (i < size) {
            R(scan[i].first, scan[i].second) = 1;
            SetZerotreeChildren(R, scan, i);
        }
    }
}

//------------------------------------------------------------------------------
// DominantPass
//------------------------------------------------------------------------------
template<typename Real>
void DominantPass(std::vector<char>& symbols,
                  std::vector<Real>& values,
                  std::vector<Index>& positions,
                  Matrix<Real> &W,
                  const std::vector<Index>& scan,
                  const Real threshold)
{
    // R is a reference matrix.
    // A non-zero value indicates that the coefficient is the descendant of a
    // zerotree root, or a zerotree root itself.
    Matrix<int> R(W.Rows(), W.Columns());

    // First encode the DC coefficient.
    if (Abs(W(0,0)) >= threshold) {
        values.push_back(Abs(W(0,0)));
        positions.push_back(scan[0]);
        if (W(0,0) > Real(0)) {
            symbols.push_back('p');
        } else {
            symbols.push_back('n');
        }
        W(0,0) = Real(0);
    } else {
        symbols.push_back('z');
    }

    // Encode the next three coefficients. These have no parent.
    for (int i=1; i<4; ++i) {
        int r = scan[i].first;
        int c = scan[i].second;

        if (Abs(W(r, c)) >= threshold) {
            values.push_back(Abs(W(r,c)));
            positions.push_back(scan[i]);
            if (W(r,c) > Real(0)) {
                symbols.push_back('p');
            } else {
                symbols.push_back('n');
            }
            W(r,c) = Real(0);
        } else {
            // The second, third, and fourth elements have no parent.
            // Check their descendants.
            if (HasSignificantChildren(W, scan, i, threshold)) {
                // Coefficient has significant descendants.
                // Encode as an isolated zero.
                symbols.push_back('z');
            } else {
                // Coefficient has no significant descendants.
                // Encode as the root of a zerotree.
                symbols.push_back('t');
                R(r,c) = 1;
                SetZerotreeChildren(R, scan, i);
            }
        }
    }

    // Encode the remaining coefficients.
    for (int i=4; i<scan.size(); ++i) {
        int r = scan[i].first;
        int c = scan[i].second;

        if (Abs(W(r,c)) >= threshold) {
            values.push_back(Abs(W(r,c)));
            positions.push_back(scan[i]);
            if (W(r,c) > Real(0)) {
                symbols.push_back('p');
            } else {
                symbols.push_back('n');
            }
            W(r,c) = Real(0);
        } else if (!R(r,c)) {
            if (HasSignificantChildren(W, scan, i, threshold)) {
                // Encode as an isolated zero.
                symbols.push_back('z');
            } else {
                // Encode as the root of a zerotree.
                symbols.push_back('t');
                R(r,c) = 1;
                SetZerotreeChildren(R, scan, i);
            }
        }
    }
}

//------------------------------------------------------------------------------
// SubordinatePass
//------------------------------------------------------------------------------
template<typename Real>
void SubordinatePass(std::vector<bool> &bits,
                     std::vector<Real> &values,
                     const Real threshold)
{
    for (unsigned i=0; i<values.size(); ++i) {
        Real prevThreshold = Real(2.0f)*threshold;
        if (values[i] >= prevThreshold) {
            values[i] -= prevThreshold;
        }
        if (values[i] >= threshold) {
            bits.push_back(true);
        } else {
            bits.push_back(false);
        }
    }
}

//------------------------------------------------------------------------------
// FindThreshold
//------------------------------------------------------------------------------
template<typename Real>
void FindThreshold(Real &threshold,
                   int &lastLevel,
                   const Matrix<Real> &W)
{
    // Find the maximum
    Real maximum(-INFINITY);
    for (int i=0; i<W.Rows()*W.Columns(); ++i) {
        if (Abs(W[i]) > maximum) {
            maximum = Abs(W[i]);
        }
    }

    // Find the powers of two that bookend the maximum
    for (int i=0; i<20; ++i) {
        Real p = Real(1 << i);
        Real p2 = Real(1 << (i+1));
        if (maximum > p && maximum < p2) {
            threshold = p;
            lastLevel = i;
        }
    }
}

//------------------------------------------------------------------------------
// ZerotreeEncode
//------------------------------------------------------------------------------
template<typename Real>
bool ZerotreeEncode(const std::string& filename,
                    const Matrix<Real>& M)
{
    assert(M.Rows() == M.Columns());
    assert(IsPowerOfTwo(M.Rows()));

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

    if (!file.is_open()) {
        std::cout << "Could not open file: " << filename << std::endl;
        return false;
    }

    int size = M.Rows();

    // Wavelet-transform M
    Matrix<Real> W(size, size);
    {
        Matrix<Real> H(size, size);
        Haar(H, size);

        Matrix<Real> Ht(size, size);
        Transpose(Ht, H);

        Matrix<Real> T(size, size);
        Multiply(T, H, M);
        Multiply(W, T, Ht);
    }
    W.Set(M);

    // Write the size of the matrix to file.
    file.write(reinterpret_cast<char*>(&size), sizeof(int));

    // Determine the threshold
    Real threshold(0);
    int numPasses(0);
    FindThreshold(threshold, numPasses, W);

    // Write the number of passes done.
    int passCount(0);
    file.write(reinterpret_cast<char*>(&numPasses), sizeof(int));
    numPasses++;

    unsigned byteCount = 8;
    // Write the threshold.
    file.write(reinterpret_cast<char*>(&threshold), sizeof(Real));
    byteCount += sizeof(Real);

    // Determine the scan order
    Matrix<int> S;
    MortonScan(S, size);
    std::vector<Index> scan;
    VectorizeScan(scan, S);

    // Encode
    //std::cout << "\nEncoding...\n";
    std::vector<char>   symbols;
    std::vector<Real>   values;
    std::vector<Index>  positions;
    std::vector<bool>   bits;
    for (;;) {
        symbols.clear();
        bits.clear();
        positions.clear();

        passCount++;
        if (passCount == numPasses) {
            break;
        }

        DominantPass(symbols, values, positions, W, scan, threshold);
        threshold *= Real(0.5);
        SubordinatePass(bits, values, threshold);

        // Write the symbol stream to the file
        int symbolCount = int(symbols.size());
        file.write(reinterpret_cast<char*>(&symbolCount), sizeof(int));
        byteCount += 4;

        int bytePosition = 0;
        unsigned char byte = 0;
        for (unsigned i=0; i<symbols.size(); ++i) {
            switch (symbols[i]) {
                case 'p':
                    // p encodes as 0.
                    bytePosition += 1;
                    break;
                case 't':
                    // t encodes as 10.
                    bytePosition += 2;
            }
        }

        symbolCount = 0;
        byte = 0;
        for (unsigned i=0; i<symbols.size(); ++i) {
            switch (symbols[i]) {
                case 't':
                    // t encodes as 0.
                    // Do nothing.
                    break;
                case 'n':
                    // n encodes as 1.
                    byte += 1;
                    break;
                case 'p':
                    // p encodes as 2.
                    byte += 2;
                    break;
                case 'z':
                    // z encodes as 3.
                    byte += 3;
                    break;
            }
            symbolCount++;
            if (symbolCount == 4) {
                file.write(reinterpret_cast<char*>(&byte),
                           sizeof(unsigned char));
                byteCount++;
                symbolCount = 0;
                byte = 0;
            } else {
                byte <<= 2;
            }
        }
        if (symbolCount > 0) {
            file.write(reinterpret_cast<char*>(&byte),
                       sizeof(unsigned char));
            byteCount++;
        }

        // Write the subordinate bits to file.
        int bitCount = int(bits.size());
        file.write(reinterpret_cast<char*>(&bitCount), sizeof(int));
        byteCount += 4;

        byte = 0;
        symbolCount = 0;
        for (unsigned i=0; i<bits.size(); ++i) {
            if (bits[i]) {
                byte += 1;
            }

            if (symbolCount < 7) {
                byte <<= 1;
                symbolCount++;
            } else {
                file.write(reinterpret_cast<char*>(&byte),
                           sizeof(unsigned char));
                byteCount++;
                symbolCount = 0;
                byte = 0;
            }
        }
        if (symbolCount > 0) {
            byte <<= (7-symbolCount);
            file.write(reinterpret_cast<char*>(&byte),
                       sizeof(unsigned char));
            byteCount++;
        }

        if (0)
        {
            // Echo results
            std::cout << "\tSymbols: ";
            for (unsigned i=0; i<symbols.size(); ++i) {
                std::cout << symbols[i];
            }
            std::cout << std::endl;
            std::cout << "\tBits: ";
            for (unsigned i=0; i<bits.size(); ++i) {
                if (bits[i]) {
                    std::cout << 1;
                } else {
                    std::cout << 0;
                }
            }
            std::cout << std::endl;
        }
    }
    //std::cout << byteCount << " bytes written." << std::endl;
    return true;
}

//------------------------------------------------------------------------------
// ZerotreeDecode
//------------------------------------------------------------------------------
template<typename Real>
bool ZerotreeDecode(Matrix<Real> &M, const std::string& filename)
{
    std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cout << "Could not open file: " << filename << std::endl;
        return false;
    }
    //std::cout << "Decoding...\n";

    int size = 0;
    file.read(reinterpret_cast<char*>(&size), sizeof(int));

    // Resize M
    M.Resize(size, size);

    // Reference matrix
    Matrix<int> R(size, size);

    // Determine the scan order
    Matrix<int> S;
    MortonScan(S, size);
    std::vector<Index> scan;
    VectorizeScan(scan, S);

    // Read the number of passes done.
    int passes = 0;
    file.read(reinterpret_cast<char*>(&passes), sizeof(int));

    // Read the initial threshold.
    Real threshold(0);
    file.read(reinterpret_cast<char*>(&threshold), sizeof(Real));

    // Push-back all significant values for each pass, so we can refine
    // them in each subsequent pass.
    std::vector<Index> significants;

    // Read information per-pass
    int bytesRead = 12;
    for (int i=0; i<passes; ++i) {
        unsigned char byte;
        std::vector<char> symbols;

        // Read the t, p, n, and z symbols.
        int symbolCount = 0;
        file.read(reinterpret_cast<char*>(&symbolCount), sizeof(int));
        bytesRead += 4;
        while (symbolCount > 0) {
            file.read(reinterpret_cast<char*>(&byte),
                      sizeof(unsigned char));
            bytesRead++;

            for (int j=3; j>=0; --j) {
                int p = 3 & (byte >> 2*j);
                switch (p) {
                    case 0:
                        symbols.push_back('t');
                        break;
                    case 1:
                        symbols.push_back('n');
                        break;
                    case 2:
                        symbols.push_back('p');
                        break;
                    case 3:
                        symbols.push_back('z');
                        break;
                }
            }
            symbolCount -= 4;
        }

        // Read the subordinate bits.
        symbolCount = 0;
        byte = 0;
        file.read(reinterpret_cast<char*>(&symbolCount), sizeof(int));
        bytesRead += 4;

        std::vector<bool> bits;
        while (symbolCount > 0) {
            file.read(reinterpret_cast<char*>(&byte),
                      sizeof(unsigned char));
            bytesRead += 1;

            int min = 0;
            if (symbolCount < 8) {
                min = 8-symbolCount;
            }

            for (int j=7; j>=min; --j) {
                if (byte & (1<<j)) {
                    bits.push_back(true);
                } else {
                    bits.push_back(false);
                }
            }
            symbolCount -= 8;
        }

        if (0)
        {
            // Echo results
            std::cout << "\tSymbols: ";
            for (unsigned j=0; j<symbols.size(); ++j) {
                std::cout << symbols[j];
            }
            std::cout << std::endl;
            std::cout << "\tBits: ";
            for (unsigned j=0; j<bits.size(); ++j) {
                if (bits[j]) {
                    std::cout << 1;
                } else {
                    std::cout << 0;
                }
            }
            std::cout << std::endl;
        }

        // Rebuild the matrix.
        // Dominant pass
        unsigned j=0;
        unsigned k=0;
        while (j<symbols.size()) {
            int r = scan[k].first;
            int c = scan[k].second;

            if (!R(r,c)) {
                if (symbols[j] == 'p') {
                    if (i == passes-1) {
                        M(r,c) = threshold;
                    } else {
                        M(r,c) = Real(1.5) * threshold;
                    }
                    significants.push_back(scan[k]);
                } else if (symbols[j] == 'n') {
                    if (i == passes-1) {
                        M(r,c) = -threshold;
                    } else {
                        M(r,c) = -Real(1.5) * threshold;
                    }
                    significants.push_back(scan[k]);
                } else if (symbols[j] == 't') {
                    SetZerotreeChildren(R, scan, k);
                }
                R(r,c) = 1;
                j++;
            }
            k++;
        }

        // Subordinate pass
        threshold *= Real(0.5f);
        if (i < passes-1) {
            for (unsigned j=0; j<bits.size(); ++j) {
                int r = significants[j].first;
                int c = significants[j].second;

                Real value = Real(0.5) * threshold;
                if (bits[j] == 1) {
                    if (M(r,c) > Real(0.0)) {
                        M(r,c) += value;
                    } else {
                        M(r,c) -= value;
                    }
                } else {
                    if (M(r,c) > Real(0.0)) {
                        M(r,c) -= value;
                    } else {
                        M(r,c) += value;
                    }
                }
            }
        } else {
            //std::cout << "Normalizing...\n";
            for (unsigned i=0; i<M.Rows()*M.Columns(); ++i) {
                if (int(M[i]) % 2 != 0) {
                    if (M[i] < 0.0) {
                        M[i] += 1.0;
                    } else {
                        M[i] -= 1.0;
                    }
                }
            }
        }

        // Clear reference matrix and vectors for next pass.
        R.SetZero();
        symbols.clear();
        bits.clear();
    }

    //std::cout << bytesRead << " bytes read." << std::endl;
    return true;
}



//------------------------------------------------------------------------------
// ZerotreeEncode
//------------------------------------------------------------------------------
bool ZerotreeEncode(const std::string& filename,
                    const CMatrix& C)
{
    assert(C.M.Rows() == C.M.Columns());
    assert(IsPowerOfTwo(C.M.Rows()));

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

    if (!file.is_open()) {
        std::cout << "Could not open file: " << filename << std::endl;
        return false;
    }

    int size = C.M.Rows();

    Matrix<float> M;
    Cast(M, C.M);

    // Wavelet-transform M
    Matrix<float> W(size, size);
    {
        Matrix<float> H(size, size);
        Haar(H, size);

        Matrix<float> Ht(size, size);
        Transpose(Ht, H);

        Matrix<float> T(size, size);
        Multiply(T, H, M);
        Multiply(W, T, Ht);
    }

    // Write the information of CMatrix to file.
    float scale = C.scale;
    float bias = C.bias;
    float value = C.bias;
    int qbits = C.qbits;
    int rows = C.rows;
    int columns = C.columns;
    file.write(reinterpret_cast<char*>(&scale), sizeof(float));
    file.write(reinterpret_cast<char*>(&bias), sizeof(float));
    file.write(reinterpret_cast<char*>(&value), sizeof(float));
    file.write(reinterpret_cast<char*>(&qbits), sizeof(int));
    file.write(reinterpret_cast<char*>(&rows), sizeof(int));
    file.write(reinterpret_cast<char*>(&columns), sizeof(int));

    // Determine the threshold
    float threshold(0);
    int numPasses(0);
    FindThreshold(threshold, numPasses, W);

    // Write the number of passes done.
    int passCount(0);
    file.write(reinterpret_cast<char*>(&numPasses), sizeof(int));
    numPasses++;

    unsigned byteCount = 8;

    // Write the threshold.
    file.write(reinterpret_cast<char*>(&threshold), sizeof(float));
    byteCount += sizeof(float);

    // Determine the scan order
    Matrix<int> S;
    MortonScan(S, size);
    std::vector<Index> scan;
    VectorizeScan(scan, S);

    // Encode
    std::cout << "\nEncoding...\n";
    std::vector<char>   symbols;
    std::vector<float>  values;
    std::vector<Index>  positions;
    std::vector<bool>   bits;

    std::cout << "Total passes: " << numPasses << std::endl;
    for (;;) {
        symbols.clear();
        bits.clear();
        positions.clear();

        passCount++;
        std::cout << "Pass: " << passCount << std::endl;
        if (passCount == numPasses) {
            break;
        }

        DominantPass(symbols, values, positions, W, scan, threshold);
        threshold *= 0.5f;
        SubordinatePass(bits, values, threshold);

        // Write the symbol stream to the file
        int symbolCount = int(symbols.size());
        file.write(reinterpret_cast<char*>(&symbolCount), sizeof(int));
        byteCount += 4;

        symbolCount = 0;
        unsigned char byte = 0;
        for (unsigned i=0; i<symbols.size(); ++i) {
            switch (symbols[i]) {
                case 't':
                    // t encodes as 0.
                    // Do nothing.
                    break;
                case 'n':
                    // n encodes as 1.
                    byte += 1;
                    break;
                case 'p':
                    // p encodes as 2.
                    byte += 2;
                    break;
                case 'z':
                    // z encodes as 3.
                    byte += 3;
                    break;
            }
            symbolCount++;
            if (symbolCount == 4) {
                file.write(reinterpret_cast<char*>(&byte),
                           sizeof(unsigned char));
                byteCount++;
                symbolCount = 0;
                byte = 0;
            } else {
                byte <<= 2;
            }
        }
        if (symbolCount > 0) {
            file.write(reinterpret_cast<char*>(&byte),
                       sizeof(unsigned char));
            byteCount++;
        }

        // Write the subordinate bits to file.
        int bitCount = int(bits.size());
        file.write(reinterpret_cast<char*>(&bitCount), sizeof(int));
        byteCount += 4;

        byte = 0;
        symbolCount = 0;
        for (unsigned i=0; i<bits.size(); ++i) {
            if (bits[i]) {
                byte += 1;
            }

            if (symbolCount < 7) {
                byte <<= 1;
                symbolCount++;
            } else {
                file.write(reinterpret_cast<char*>(&byte),
                           sizeof(unsigned char));
                byteCount++;
                symbolCount = 0;
                byte = 0;
            }
        }
        if (symbolCount > 0) {
            byte <<= (7-symbolCount);
            file.write(reinterpret_cast<char*>(&byte),
                       sizeof(unsigned char));
            byteCount++;
        }

        // Echo results
        //std::cout << "\tSymbols: ";
        //for (unsigned i=0; i<symbols.size(); ++i) {
        //  std::cout << symbols[i];
        //}
        //std::cout << std::endl;
        //std::cout << "\tBits: ";
        //for (unsigned i=0; i<bits.size(); ++i) {
        //  if (bits[i]) {
        //      std::cout << 1;
        //  } else {
        //      std::cout << 0;
        //  }
        //}
        //std::cout << std::endl;
    }
    //std::cout << byteCount << " bytes written." << std::endl;
    return true;
}

#endif

