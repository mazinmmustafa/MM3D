#ifndef VECTOR_HPP
#define VECTOR_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"

namespace basic_lib{
// Definitions
class Vector{
public:
    cmplx x=0.0, y=0.0, z=0.0;
    Vector();
    ~Vector();
    Vector(const cmplx x, const cmplx y, const cmplx z);
    void create_vector(const cmplx x, const cmplx y, const cmplx z);
    double mag();
    Vector unit();
    int is_equal(const Vector vector);
};

// Functions
Vector operator + (const Vector A, const Vector B);
Vector operator - (const Vector A, const Vector B);
cmplx operator * (const Vector A, const Vector B);
Vector operator * (const Vector A, const cmplx scalar);
Vector operator * (const cmplx scalar, const Vector A);
Vector operator / (const Vector A, const cmplx scalar);
Vector operator ^ (const Vector A, const Vector B);
double mag(const Vector A);
Vector unit(const Vector A);
void disp(const Vector A);

}

#endif