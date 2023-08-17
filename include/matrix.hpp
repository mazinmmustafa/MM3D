#ifndef MATRIX_HPP
#define MATRIX_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"

namespace basic_lib{
// Definitions
struct Matrix_Size{
    int rows=0, cols=0;
};

class Matrix{
private:
    Matrix_Size matrix_size;
    cmplx **data=NULL;
    cmplx **lu_data=NULL;
    int *P_data=NULL;
    int is_allocated=FALSE;
    int is_lu_allocated=FALSE;
public:
    Matrix();
    ~Matrix();
    void allocate(const int rows, const int cols);
    void deallocate();
    void zeros();
    void ones();
    void eye();
    void save(const char *filename);
    void load(const char *filename);
    cmplx operator () (const int i, const int j) const;
    cmplx& operator () (const int i, const int j);
    void lup();
    void solve(const Matrix &b, Matrix &x);
    void inv();
    cmplx det();
    Matrix_Size size();
    void copy(Matrix &matrix);
    void copy_lu(Matrix &matrix);
    void normalize();
};

class Data_Matrix: public Matrix{
    
};

// Functions

}

#endif