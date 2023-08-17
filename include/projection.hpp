#ifndef PROJECTION_HPP
#define PROJECTION_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "quadl.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "shape.hpp"
#include "engine.hpp"

using namespace basic_lib;

// Definitions
struct Projection_3D_Para{
    double d;
    double l_m, l_p;
    double R_m, R_p;
    double R0;
    Vector u;
    Vector P0;
    Vector p0;
};

struct Projection_1D_Para{
    double l_m=0.0, l_p=0.0;
    double P_m=0.0, P_p=0.0;
    double d=0.0;
    Vector P0;
    Vector p0;
};

// Functions
void projection_3D(const Triangle T, const Vector p, Projection_3D_Para &para);
cmplx T1_projection(const cmplx alpha, const cmplx beta,
    Basis *B_m, Basis *B_n, const char s_m, const char s_n);
Vector T2_projection(const cmplx alpha, const cmplx beta,
    Basis *B_m, Basis *B_n, const char s_m, const char s_n);
    
#endif