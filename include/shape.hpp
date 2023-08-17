#ifndef SHAPE_HPP
#define SHAPE_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "vector.hpp"
#include "list.hpp"

using namespace basic_lib;

// Definitions
struct Triangle{
    Vector v1, v2, v3;
    Vector n;
};

struct Basis{
    Vector r_m, r_p, e_m, e_p;
    Vector n_m, n_p;
    double L=0;
    double A_m=0.0, A_p=0.0;
    Vector L_mm, L_mp, L_pm, L_pp;
};

class Shape{
private:
    Timer T;
    const double mesh_tol=1.0E-10;
    int N_triangles=0;
    int N_bases=0;
    Basis *bases_list=NULL;
    int is_bases_allocated=FALSE;
    double lambda=1.0;
    int is_set=FALSE;
public:
    Shape();
    void set(const double lambda);
    void mesh(const int is_log);
    void load_mesh();
    void log(const char *filename);
    int get_N_basis(){ return this->N_bases; }
    Basis* get_basis_list(){ return this->bases_list; }
    ~Shape();
    void create_sphere(const double radius, const double lc);
    void create_sheet(const double l, const double w, const double lc);
};


// Functions
double get_area(Vector v1, Vector v2, Vector v3);
void get_basis(Basis &B);

#endif