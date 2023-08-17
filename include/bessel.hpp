#ifndef BESSEL_HPP
#define BESSEL_HPP

// Libraries
#include "basic_lib.hpp"

namespace basic_lib{
// Definitions

// Functions
cmplx besselj(double n, cmplx z);
cmplx bessely(double n, cmplx z);
cmplx besseli(double n, cmplx z);
cmplx besselk(double n, cmplx z);
cmplx besselh(int k, double n, cmplx z);

}

#endif