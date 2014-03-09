//------------------------------------------------------------------------------
// Embedded Zerotree Wavelet Compressor
// Tim Thirion
// thirion@cs.unc.edu
// Spring 2006
//------------------------------------------------------------------------------
#define USE_EXAMPLE_MATRIX 0
#include "zerotree.hpp"
#include "cmatrix.hpp"

typedef Matrix<float> Matrixf;

const int SIZE = 8;

//------------------------------------------------------------------------------
// ReadDTED
//------------------------------------------------------------------------------
template<typename Real>
bool ReadDTED(Matrix<Real>& M, const std::string &filename)
{
    std::fstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    int rows, columns;
    file >> rows >> columns;

    M.Resize(rows, columns);
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            int height;
            file >> height;
            M(i, j) = float(height);
            if (M(i, j) < Real(0)) {
                M(i, j) = Real(0);
            }
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// ReadMatrixFromFile
//------------------------------------------------------------------------------
template<typename Real>
bool ReadMatrixFromFile(Matrix<Real>& M, const std::string &filename)
{
    std::fstream file(filename.c_str(), std::ios::in);
    if (!file.is_open()) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    int rows, columns;
    file >> rows >> columns;

    M.Resize(rows, columns);
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            float height;
            file >> height;
            M(i, j) = height;
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// WriteMatrixToFile
//------------------------------------------------------------------------------
template<typename Real>
bool WriteMatrixToFile(const Matrix<Real> &M, const std::string& filename) {
    std::ofstream file(filename.c_str());

    if (!file.is_open()) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    file << M.Rows() << "\n";
    file << M.Columns() << "\n";
    for (int r=0; r<M.Rows(); ++r) {
        for (int c=0; c<M.Columns(); ++c) {
            file << M(r,c) << " ";
        }
        file << "\n";
    }
    return true;
}

//------------------------------------------------------------------------------
// ReadCMatrixFromFile
//------------------------------------------------------------------------------
bool ReadCMatrixFromFile(CMatrix &C, const std::string &filename)
{
    std::fstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    file >> C.scale;
    file >> C.bias;
    file >> C.value;
    file >> C.qbits;
    file >> C.rows;
    file >> C.columns;

    std::cout << C.scale << "\n";
    std::cout << C.bias << "\n";
    std::cout << C.value << "\n";
    std::cout << C.qbits << "\n";
    std::cout << C.rows << "\n";
    std::cout << C.columns << "\n";

    // FIX! - assume #rows is a power of two and > than #columns
    //C.M.Resize(rows, columns);
    C.M.Resize(C.rows, C.rows);
    for (int i=0; i<C.rows; ++i) {
        for (int j=0; j<C.columns; ++j) {
            unsigned height;
            file >> height;
            C.M(i, j) = height;
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// WriteCMatrixToFile
//------------------------------------------------------------------------------
bool WriteCMatrixToFile(CMatrix &C, const std::string &filename)
{
    std::fstream file(filename.c_str());

    if (!file.is_open()) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    file << C.scale << "\n";
    file << C.bias << "\n";
    file << C.value << "\n";
    file << C.qbits << "\n";
    file << C.M.Rows() << "\n";
    file << C.M.Columns() << "\n";
    for (int r=0; r<C.M.Rows(); ++r) {
        for (int c=0; c<C.M.Columns(); ++c) {
            file << C.M(r,c) << " ";
        }
        file << "\n";
    }
    return true;
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main (int argc, char * const argv[]) {
    std::cout.setf(std::ios::right);
    std::cout << std::setprecision(4) << std::fixed << std::showpoint;
    std::cout << std::setfill(' ');

    if (argc == 4 && !strncmp(argv[1], "-c", 2))
    {
        // Compress
        Matrixf M;
        ReadMatrixFromFile(M, std::string(argv[2]));
        Scale(M, 2.0f, M);
        ZerotreeEncode(std::string(argv[3]), M);
    }
    else if (argc == 4 && !strncmp(argv[1], "-d", 2))
    {
        // Decompress
        Matrixf M;
        ZerotreeDecode(M, std::string(argv[2]));
        WriteMatrixToFile(M, std::string(argv[3]));
    }
    else
    {
        printf("[-c|-d] [input] [output]\n");
    }
    return 0;
}
