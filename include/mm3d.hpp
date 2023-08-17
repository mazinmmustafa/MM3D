#ifndef MM3D_HPP
#define MM3D_HPP

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
struct RCS{
    double theta=0.0, phi=0.0;
};

struct Scattered_Field{
    Vector E, H;
};

class Engine{
private:
public:
    QuadL quadl;
    int is_quadl_allocated=FALSE;
    int is_shape_allocated=FALSE;
    int is_Z_mn_allocated=FALSE;
    int is_V_m_allocated=FALSE;
    int is_I_n_allocated=FALSE;
    Matrix Z_mn, V_m, I_n;
    double lambda, lc, near_term;
    Shape shape;
    int N_basis=0;
    Basis *basis_list=NULL;
    //
    Engine(const double lambda, const double lc, const double near_term);
    ~Engine();
    void set_quadl(const int N_quadl, const double tol, const int k_max);
    void compute_Z_mn();
    void test_Z_mn(const int m, const int n);
    void save_Z_mn(const char *filename);
    void load_Z_mn(const char *filename);
    void mesh(const int is_log);
    void load_mesh(const int is_log);
    void set_Zmn_element(const int m, const int n, const cmplx Z);
    cmplx get_Zmn_element(const int m, const int n);
    void compute_V_m_plane_wave(const cmplx E_TM, const cmplx E_TE,
        const double theta_i, const double phi_i);
    void save_V_m(const char *filename);
    void load_V_m(const char *filename);
    void compute_I_n();
    RCS compute_rcs(const double theta_s, const double phi_s);
    void export_currents();
    Scattered_Field compute_near_field(Vector r);
    // Shapes
    void create_sphere(const double radius);
    void create_sheet(const double l, const double w);
};

// Functions


#endif