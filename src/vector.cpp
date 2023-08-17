//
#include "vector.hpp"

namespace basic_lib{
// Begin library

Vector::Vector(){}

Vector::~Vector(){}

Vector::Vector(const cmplx x, const cmplx y, const cmplx z){
    this->x = x;
    this->y = y;
    this->z = z;
}

void Vector::create_vector(const cmplx x, const cmplx y, const cmplx z){
    this->x = x;
    this->y = y;
    this->z = z;
}

double Vector::mag(){
    return sqrt(pow(abs(this->x), 2.0)+pow(abs(this->y), 2.0)+pow(abs(this->z), 2.0));
}

Vector Vector::unit(){
    double mag_=Vector::mag();
    return Vector(this->x/mag_, this->y/mag_, this->z/mag_);
}

int Vector::is_equal(const Vector vector){
    return ((round_double(real(this->x), 14)==round_double(real(vector.x), 14))&&
            (round_double(real(this->y), 14)==round_double(real(vector.y), 14))&&
            (round_double(real(this->z), 14)==round_double(real(vector.z), 14))&&
            (round_double(imag(this->x), 14)==round_double(imag(vector.x), 14))&&
            (round_double(imag(this->y), 14)==round_double(imag(vector.y), 14))&&
            (round_double(imag(this->z), 14)==round_double(imag(vector.z), 14))) 
            ? TRUE : FALSE;
}

Vector operator + (const Vector A, const Vector B){
    return Vector(A.x+B.x, A.y+B.y, A.z+B.z);
}

Vector operator - (const Vector A, const Vector B){
    return Vector(A.x-B.x, A.y-B.y, A.z-B.z);
}

cmplx operator * (const Vector A, const Vector B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

Vector operator * (const Vector A, const cmplx scalar){
    return  Vector(scalar*A.x, scalar*A.y, scalar*A.z);
}

Vector operator * (const cmplx scalar, const Vector A){
    return  Vector(scalar*A.x, scalar*A.y, scalar*A.z);
}

Vector operator / (const Vector A, const cmplx scalar){
    return  Vector(A.x/scalar, A.y/scalar, A.z/scalar);
}

Vector operator ^ (const Vector A, const Vector B){
    return Vector(A.y*B.z-A.z*B.y, 
                  A.z*B.x-A.x*B.z,
                  A.x*B.y-A.y*B.x);
}

double mag(const Vector A){
    return sqrt(pow(abs(A.x), 2.0)+pow(abs(A.y), 2.0)+pow(abs(A.z), 2.0));
}

Vector unit(const Vector A){
    double mag_=mag(A);
    return Vector(A.x/mag_, A.y/mag_, A.z/mag_);
}

void disp(const Vector V){
    printf("<(%21.14E, %21.14E),\n (%21.14E, %21.14E),\n (%21.14E, %21.14E)>\n", 
        real(V.x), imag(V.x), real(V.y), imag(V.y), real(V.z), imag(V.z));
}

// End of library
}