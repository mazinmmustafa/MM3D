#ifndef ENGINE_HPP
#define ENGINE_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "quadl.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "shape.hpp"
#include "projection.hpp"

using namespace basic_lib;

// Definitions
struct Singular_Para{
    Vector A, B, C;
    double a=0.0, b=0.0, c=0.0;
    double Am=0.0;
    double sign=0.0;
};

// Functions
cmplx compute_Zmn(Basis *B_m, Basis *B_n, const double near_term, QuadL *quadl, int &flag, const int scenario);
cmplx I1_singular(const Singular_Para *para);
cmplx I2_singular(const Singular_Para *para);
cmplx I3_singular(const Singular_Para *para);
cmplx I4_singular(const Singular_Para *para);
void I_radiation(cmplx s1, cmplx s2, cmplx &I1, cmplx &I2);

#endif