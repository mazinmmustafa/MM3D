#ifndef UTILITIES_HPP
#define UTILITIES_HPP

// Libraries
#include "basic_lib.hpp"

namespace basic_lib{
// Definitions
class File{
private:
    FILE *file_ptr=NULL;
    char *filename=NULL;
    char *option=NULL;
    int is_open=FALSE;
public:
    File();
    ~File();
    void open(const char *filename, const char *option);
    void close();
    void write(const char *format, ...);
    int read(const char *format, ...);
};

class Timer{
private:
    int is_set=FALSE;
    #ifdef _WIN64
    time_t start, stop;
    #endif
    #ifdef __linux__
    struct timespec start, stop;
    #endif
    double elapsed=0.0;
public:
    Timer();
    ~Timer();
    void set();
    void unset();
    void unset_silent();
    double get_elapsed();
};

class Random{
private:
    time_t t;
    void set_seed();
public:
    Random();
    ~Random();
    int rand_int(const int a, const int b);
    double rand_double(const double a, const double b);
};

class Range{
private:
    int Ns=0;
    int is_allocated=FALSE;
    double *data=NULL;
public:
    Range();
    ~Range();
    void linspace(const double x_min, const double x_max, const int Ns);
    void logspace(const double x_min, const double x_max, const int Ns);
    double operator() (const int i) const;
    void deallocate();
    int size();
};

// Functions
void check_error(const int condition, const char *error_msg);
void progress_bar(const unsigned long i, const unsigned long N, const char *msg);
/*
#ifdef _WIN64
void usleep(__int64 usec);
#endif
*/
void disp(const cmplx z);
void disp(const double x);
void disp(const int n);
double deg2rad(const double x);
double rad2deg(const double x);
double sinc(double x);
cmplx sinc(cmplx z);
double round_double(const double x, const int n);

}

#endif