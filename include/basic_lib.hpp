#ifndef BASIC_LIB_HPP
#define BASIC_LIB_HPP

// Libraries
#include <iostream>
#include <fstream>
#include <complex>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

// Definitions
#define TRUE 1
#define FALSE 0
#define cmplx std::complex<double>
#define str std::string

// Constants
namespace basic_lib{
const double pi=3.141592653589793;
const double c0=299792458.0;
const double eps0=8.8541878128E-12;
const double mu0=1.25663706212E-6;
const double eta0=376.730313668;
}

#endif