//
#include "bessel.hpp"

namespace basic_lib{
// Begin library

// Fortran Functions
extern "C" void zbesj_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesy_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesi_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesk_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesh_f77_(int *k, double *fnu, double *zr, double *zi, double *cyr, double *cyi);

cmplx besselj(double n, cmplx z){
	double zr=real(z);
	double zi=imag(z);
	double cyr, cyi;
	cmplx j=cmplx(0.0, 1.0);
	zbesj_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

cmplx bessely(double n, cmplx z){
	double zr=real(z);
	double zi=imag(z);
	double cyr, cyi;
	cmplx j=cmplx(0.0, 1.0);
	zbesy_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

cmplx besseli(double n, cmplx z){
	double zr=real(z);
	double zi=imag(z);
	double cyr, cyi;
	cmplx j=cmplx(0.0, 1.0);
	zbesi_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

cmplx besselk(double n, cmplx z){
	double zr=real(z);
	double zi=imag(z);
	double cyr, cyi;
	cmplx j=cmplx(0.0, 1.0);
	zbesk_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

cmplx besselh(int k, double n, cmplx z){
	assert(k==1||k==2);
	double zr=real(z);
	double zi=imag(z);
	double cyr, cyi;
	cmplx j=cmplx(0.0, 1.0);
	zbesh_f77_(&k, &n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

// End of library
}