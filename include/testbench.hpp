#ifndef TESTBENCH_HPP
#define TESTBENCH_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "bessel.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "engine.hpp"
#include "mm3d.hpp"

// Definitions

// Functions
void test_utilities();
void test_range();
void test_bessel();
void test_matrix();
void test_vector();
void test_quadl();
void test_shape();
void test_integrals();
void test_engine();
void test_RCS_1();
void test_RCS_2();
void test_RCS_3();
void test_RCS_sheet();
void test_current_sheet();
void test_sheet_near_field();

#endif