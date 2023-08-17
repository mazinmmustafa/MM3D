//
#include "testbench.hpp"

using namespace basic_lib;

void test_utilities(){

    printf("Testing utilities:\n");

    check_error(1==0, "This is not an error!");

    File file;

    file.open("data/misc/test_file.txt", "w");
    file.write("%21.14E %21.14E %21.14E\n", pi, c0, mu0);
    file.write("%21.14E %21.14E %21.14E\n", eps0, pi, eta0);
    file.close();

    file.open("data/misc/test_file.txt", "a");
    file.write("The next line is appended:\n");
    file.write("%21.14E %21.14E %21.14E\n", pi, c0, mu0);
    file.close();

    file.open("data/misc/test_file.txt", "r");
    double x, y;
    file.read("%lf %lf", &x, &y);
    printf("%21.14E %21.14E\n", x, y);
    file.read("%lf", &x);
    printf("%21.14E\n", x);
    file.read("%lf", &x);
    printf("%21.14E\n", x);
    file.close();

    Timer T;
    T.unset();
    T.set();
    
    int N=600;
    Random random_gen;
    double r_min=-4.0, r_max=+3.0;

    File new_file;
    new_file.open("data/misc/randon.dat", "w");

    T.set();
    for (int i=0; i<N; i++){
        new_file.write("%3d %21.14E\n", i, random_gen.rand_double(r_min, r_max));
        usleep(1000);
        progress_bar(i, N, " dummy loop");
    }
    T.unset();

    new_file.close();

}

void test_range(){
    
    printf("Testing range:\n");

    Range x_range;
    File file;

    x_range.linspace(0.0, 10.0, 60);
    file.open("data/test_basic_lib/test_range/data1.dat", "w");
    for (int i=0; i<x_range.size(); i++){
        double x, y;
        x = x_range(i);
        y = sin(sin(x)+x);
        file.write("%21.14E %21.14E\n", x, y);
    }
    file.close();

    x_range.logspace(1.0E-3, 1.0E+1, 50);
    file.open("data/test_basic_lib/test_range/data2.dat", "w");
    for (int i=0; i<x_range.size(); i++){
        double x, y;
        x = x_range(i);
        y = sin(x)*exp(-x);
        file.write("%21.14E %21.14E\n", x, y);
    }
    file.close();

}

void test_bessel(){

    printf("Testing bessel:\n");

    int Ns=100;
    double R=4.0;

    Range x_range, y_range;
    x_range.linspace(-R, +R, Ns);
    y_range.linspace(-R, +R, Ns);

    File file;
    file.open("data/test_basic_lib/test_bessel/data1.dat", "w");
    double x, y;
    cmplx z;
    for (int i=0; i<x_range.size(); i++){
        x = x_range(i);
        for (int j=0; j<y_range.size(); j++){
            y = y_range(j);
            z = cmplx(x, y);
            file.write("%21.14E %21.14E %21.14E\n", x, y, log10(abs(besselh(2, 1, z))));
        }
        file.write("\n");
    }
    file.close();

}

void test_matrix(){

    printf("Testing matrix:\n");

    Matrix A;
    A.allocate(4, 3);
    A.eye();

    A.save("data/misc/A.dat");

    Matrix B;
    B.allocate(4, 3);
    B.load("data/misc/A.dat");

    disp(A(0, 0));
    B(1, 0) = -1.0;
    disp(B(1, 0));
    A.deallocate();

    A.allocate(3, 3);
    A(0, 0) = +0.0; A(0, 1) = +1.0; A(0, 2) = +2.0; 
    A(1, 0) = -4.0; A(1, 1) = +0.0; A(1, 2) = +0.0; 
    A(2, 0) = +6.0; A(2, 1) = +0.0; A(2, 2) = +1.0; 

    A.lup();
    A.lup();
    B.allocate(3, 3);
    A.copy_lu(B);
    B.save("data/misc/LU.dat");

    Matrix b, x;
    b.allocate(3, 1);
    x.allocate(3, 1);
    b(0, 0) = +1.0; 
    b(1, 0) = +2.0;
    b(2, 0) = +3.0;

    A.solve(b, x);
    x.save("data/misc/x.dat");
    disp(A.det());
    A.inv();
    A.save("data/misc/inv.dat");
    disp(A.det());

}

void test_vector(){

    printf("Testing vector:\n");

    Vector A(1.0, 2.0, 3.0);
    Vector B=A;
    disp(A^B);
    disp(A+B);
    disp(A-B);
    disp(A+4.0*B/3.0-0.5*A.unit());
    disp(A*B);
    disp(mag(A)*mag(A));
    disp(pow(A*unit(A), 2.0));
    
}

struct Func_Args{
    double a;
};

cmplx func_1D(const cmplx x, void *args){
    Func_Args *Args=(Func_Args*)args;
    double a=Args->a;
    return exp(-a*a*x*x)+sin(x);
}

cmplx func_2D(const cmplx x, const cmplx y, void *args){
    Func_Args *Args=(Func_Args*)args;
    double a=Args->a;
    return cos(x+1.0)*log(2.0-y)+exp(-a*a*(pow(x-0.2, 2.0)+pow(y-0.2, 2.0)));
}

cmplx func_2D_tri(const cmplx x, const cmplx y, void *args){
    assert(args!=NULL);
    return cos(x+1.0)*log(2.0-y);
}

cmplx func_4D_tri(const cmplx x, const cmplx y, const cmplx x_, const cmplx y_, void *args){
    assert(args!=NULL);
    return cos(x_+y)*cos(y_)*x;
}

void test_quadl(){ 

    printf("Testing quadl:\n");

    int N_quadl=3;
    double tol=1.0E-12;
    int k_max=15;
    int flag;
    QuadL quadl;
    quadl.set(N_quadl, tol, k_max);
    quadl.set(N_quadl, tol, k_max);

    Func_Args args;
    Timer T;
    
    args.a = 10.0;
    T.set();
    disp(quadl.integral_1D(func_1D, &args, -2.0, +4.0, flag));
    T.unset();
    printf("flag = %d\n", flag);

    args.a = sqrt(80.0);
    T.set();
    disp(quadl.integral_2D(func_2D, &args, -2.0, 1.0, -1.0, 1.0, flag));
    T.unset();
    printf("flag = %d\n", flag);

    N_quadl = 8;
    quadl.set(N_quadl, tol, k_max);

    Triangle_2D T1;
    Vector v1(0.0, 0.0, 0.0);
    Vector v2(1.0, 0.0, 0.0);
    Vector v3(0.0, 1.0, 0.0);
    T1.v1 = v1; T1.v2 = v2; T1.v3 = v3;

    T.set();
    disp(quadl.integral_2D_tri(func_2D_tri, &args, T1, flag));
    T.unset();
    printf("flag = %d\n", flag);

    Triangle_2D T2;
    Vector v1_(0.0, 0.0, 0.0);
    Vector v2_(2.0, 0.0, 0.0);
    Vector v3_(0.0, 2.0, 0.0);
    T2.v1 = v1_; T2.v2 = v2_; T2.v3 = v3_;
    T.set();
    disp(quadl.integral_4D_tri(func_4D_tri, &args, T1, T2, flag));
    T.unset();
    printf("flag = %d\n", flag);

    quadl.unset();
    quadl.unset();

}

void test_shape(){

    // Definitions
    const double lc=0.1;    
    const double lambda=1.0;
    const double radius=0.1*lambda;

    //
    Shape shape;
    shape.set(lambda);
    shape.create_sphere(radius, lc);
    shape.mesh(TRUE);
    shape.load_mesh();
    shape.log("mesh/shape_log.txt");
    
}

void test_integrals(){

    // Definitions
    double near_term=1.0;
    int N_quadl=3;
    double tol=1.0E-4;
    int k_max=15;
    QuadL quadl;
    quadl.set(N_quadl, tol, k_max);
    quadl.set(N_quadl, tol, k_max);
    //

    Basis B_m, B_n;

    Timer T;
    int flag=0;

    // Singular integrals

    Singular_Para para;

    Vector v1=Vector(+0.0, +0.4, -0.1);
    Vector v2=Vector(-0.3, +0.0, +0.0);
    Vector v3=Vector(+0.2, +0.0, +0.2);

    para.A = v3-v2;
    para.B = v3-v1;
    para.C = v2-v1;

    para.a = para.A.mag(); 
    para.b = para.B.mag(); 
    para.c = para.C.mag();

    para.sign = +1.0;
    para.Am = get_area(para.A, para.B, para.C);

    printf("Signular term I:\n");
    T.set();
    disp(I1_singular(&para));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Signular term 2:\n");
    T.set();
    disp(I2_singular(&para));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Signular term 3:\n");
    T.set();
    disp(I3_singular(&para));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Signular term 4:\n");
    T.set();
    disp(I4_singular(&para));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario I
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.2, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.3, +0.0, +0.0);

    B_n.r_m = Vector(+0.0, -0.1, +0.0);
    B_n.e_m = Vector(+0.2, +0.0, +0.0);
    B_n.r_p = Vector(+0.0, +0.2, +0.0);
    B_n.e_p = Vector(-0.3, +0.0, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario I:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario II Case I
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.3, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.2, +0.0, +0.0);

    B_n.r_m = Vector(+0.2, +0.3, +0.0);
    B_n.e_m = Vector(+0.0, +0.2, +0.0);
    B_n.r_p = Vector(-0.2, +0.0, +0.0);
    B_n.e_p = Vector(+0.3, +0.0, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario II Case I:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario II Case III
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.3, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.2, +0.0, +0.0);

    B_n.r_m = Vector(-0.1, +0.2, +0.0);
    B_n.e_m = Vector(-0.2, +0.0, +0.0);
    B_n.r_p = Vector(+0.3, +0.0, +0.0);
    B_n.e_p = Vector(+0.0, +0.2, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario II Case III:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario II Case II
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.3, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.2, +0.0, +0.0);

    B_n.r_m = Vector(-0.2, +0.0, +0.0);
    B_n.e_m = Vector(+0.3, +0.0, +0.0);
    B_n.r_p = Vector(+0.2, +0.4, +0.0);
    B_n.e_p = Vector(+0.0, +0.2, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario II Case II:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario II Case IV
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.3, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.2, +0.0, +0.0);

    B_n.r_m = Vector(+0.3, +0.0, +0.0);
    B_n.e_m = Vector(+0.0, +0.2, +0.0);
    B_n.r_p = Vector(-0.2, +0.1, +0.0);
    B_n.e_p = Vector(-0.2, +0.0, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario II Case IV:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario III Case I
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.1, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.3, +0.0, +0.0);

    B_n.r_m = Vector(+0.3, -0.1, +0.0);
    B_n.e_m = Vector(+0.1, +0.0, +0.0);
    B_n.r_p = Vector(-0.3, +0.0, +0.0);
    B_n.e_p = Vector(+0.0, -0.1, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario III Case I:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario III Case III
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.1, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.3, +0.0, +0.0);

    B_n.r_m = Vector(-0.1, -0.3, +0.0);
    B_n.e_m = Vector(+0.0, -0.1, +0.0);
    B_n.r_p = Vector(+0.1, +0.0, +0.0);
    B_n.e_p = Vector(-0.3, +0.0, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario III Case III:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);
    
    // Scenario III Case II
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.1, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.3, +0.0, +0.0);

    B_n.r_m = Vector(-0.3, +0.0, +0.0);
    B_n.e_m = Vector(+0.0, -0.1, +0.0);
    B_n.r_p = Vector(-0.3, -0.2, +0.0);
    B_n.e_p = Vector(+0.1, +0.0, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario III Case II:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario III Case IV
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.1, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.2, +0.0);
    B_m.e_p = Vector(-0.3, +0.0, +0.0);

    B_n.r_m = Vector(+0.1, +0.0, +0.0);
    B_n.e_m = Vector(-0.3, +0.0, +0.0);
    B_n.r_p = Vector(-0.2, -0.2, +0.0);
    B_n.e_p = Vector(+0.0, -0.1, +0.0);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario III Case IV:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

    // Scenario IV Far
    
    B_m.r_m = Vector(+0.0, -0.1, +0.0);
    B_m.e_m = Vector(+0.0, +0.0, +0.0);
    B_m.r_p = Vector(+0.0, +0.1, +0.0);
    B_m.e_p = Vector(-0.1, +0.0, +0.0);

    B_n.r_m = Vector(+0.0, -0.1, +1.1);
    B_n.e_m = Vector(+0.0, +0.0, +1.1);
    B_n.r_p = Vector(+0.0, +0.1, +1.1);
    B_n.e_p = Vector(-0.1, +0.0, +1.1);

    get_basis(B_m);
    get_basis(B_n);

    printf("\nScenario IV Far:\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, FALSE));
    T.unset();
    printf("flag = %d\n", flag);

    printf("Scenario IV Near (all numerical):\n");
    T.set();
    disp(compute_Zmn(&B_m, &B_n, near_term, &quadl, flag, 6));
    T.unset();
    printf("flag = %d\n", flag);

}

void test_engine(){

    // Definition
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=5.0;
    const double near_term=10.0*lambda;
    const double lc=0.1;

    //
    Timer T;
    const double radius=0.1*lambda;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sphere(radius);
    engine.mesh(FALSE);
    engine.load_mesh(TRUE);
    //

    // int m, n;
    // m = 151;
    // n = 150;
    // T.set();
    // engine.test_Z_mn(m, n);
    // T.unset();

    T.set();
    engine.compute_Z_mn();
    engine.save_Z_mn("data/misc/Z_mn.dat");
    T.unset();

    engine.load_Z_mn("data/misc/Z_mn.dat");
    // disp(engine.get_Zmn_element(100, 100));
    
}

void test_RCS_1(){

    // Definition
    const double GHz=1.0E+9;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=0.25*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double radius=0.5;
    const int Ns=1001;
    double theta_i=deg2rad(+0.0);
    double phi_i=deg2rad(+0.0);
    cmplx E_TM, E_TE;
    const double theta_s_min=deg2rad(-180.0);
    const double theta_s_max=deg2rad(+180.0);
    double phi_s;

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sphere(radius);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    T.set();
    engine.compute_Z_mn();
    T.unset();
    // engine.save_Z_mn("data/misc/Z_mn.dat");
    // engine.load_Z_mn("data/misc/Z_mn.dat");
   
    //
    Range theta_s;
    theta_s.linspace(theta_s_min, theta_s_max, Ns);
    RCS sigma;
    File file;

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents();
    file.open("data/sphere/RCS_1/data1.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
    
    //
    E_TM = +0.0;
    E_TE = +1.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sphere/RCS_1/data2.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
}

void test_RCS_2(){

    // Definition
    const double GHz=1.0E+9;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=0.5*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double radius=0.5;
    const int Ns=1001;
    double theta_i=deg2rad(+0.0);
    double phi_i=deg2rad(+0.0);
    cmplx E_TM, E_TE;
    const double theta_s_min=deg2rad(-180.0);
    const double theta_s_max=deg2rad(+180.0);
    double phi_s;

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sphere(radius);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    T.set();
    engine.compute_Z_mn();
    T.unset();
    // engine.save_Z_mn("data/misc/Z_mn.dat");
    // engine.load_Z_mn("data/misc/Z_mn.dat");
   
    //
    Range theta_s;
    theta_s.linspace(theta_s_min, theta_s_max, Ns);
    RCS sigma;
    File file;

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sphere/RCS_2/data1.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
    
    //
    E_TM = +0.0;
    E_TE = +1.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sphere/RCS_2/data2.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
}

void test_RCS_3(){

    // Definition
    const double GHz=1.0E+9;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=0.75*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double radius=0.5;
    const int Ns=1001;
    double theta_i=deg2rad(+0.0);
    double phi_i=deg2rad(+0.0);
    cmplx E_TM, E_TE;
    const double theta_s_min=deg2rad(-180.0);
    const double theta_s_max=deg2rad(+180.0);
    double phi_s;

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sphere(radius);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    // T.set();
    // engine.compute_Z_mn();
    // T.unset();
    // engine.save_Z_mn("data/misc/Z_mn.dat");
    engine.load_Z_mn("data/misc/Z_mn.dat");
   
    //
    Range theta_s;
    theta_s.linspace(theta_s_min, theta_s_max, Ns);
    RCS sigma;
    File file;

    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    // engine.export_currents(); 

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents(); 
    file.open("data/sphere/RCS_3/data1.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
    
    //
    E_TM = +0.0;
    E_TE = +1.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sphere/RCS_3/data2.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();

}

void test_RCS_sheet(){

    // Definition
    const double GHz=1.0E+9;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=0.5*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double l=0.6;
    const double w=0.4;
    const int Ns=1001;
    double theta_i=deg2rad(+0.0);
    double phi_i=deg2rad(+0.0);
    cmplx E_TM, E_TE;
    const double theta_s_min=deg2rad(-180.0);
    const double theta_s_max=deg2rad(+180.0);
    double phi_s;

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sheet(l, w);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    T.set();
    engine.compute_Z_mn();
    T.unset();
    // engine.save_Z_mn("data/misc/Z_mn.dat");
    // engine.load_Z_mn("data/misc/Z_mn.dat");
   
    //
    Range theta_s;
    theta_s.linspace(theta_s_min, theta_s_max, Ns);
    RCS sigma;
    File file;

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sheet/RCS_1/data1.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
    
    //
    E_TM = +0.0;
    E_TE = +1.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    file.open("data/sheet/RCS_1/data2.dat", "w");
    phi_s = deg2rad(+0.0);
    for (int i=0; i<Ns; i++){
        sigma = engine.compute_rcs(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta_s(i)), 
            10.0*log10(sigma.theta*pow(lambda, 2.0)), 10.0*log10(sigma.phi*pow(lambda, 2.0)));
        progress_bar(i, Ns, "computing RCS...");
    }
    file.close();
}

void test_current_sheet(){

    // Definition
    const double GHz=1.0E+9;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=1.0*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double l=0.4;
    const double w=0.3;
    double theta_i=deg2rad(+0.0);
    double phi_i=deg2rad(+0.0);
    cmplx E_TM, E_TE;

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sheet(l, w);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    T.set();
    engine.compute_Z_mn();
    T.unset();
    engine.save_Z_mn("data/misc/Z_mn.dat");
    // engine.load_Z_mn("data/misc/Z_mn.dat");

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents();
    
}

void test_sheet_near_field(){

    // Definition
    const double GHz=1.0E+9;
    const double cm=1.0E-2;
    const int N_quadl=3;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lc=0.1;
    //
    const double freq=1.0*GHz;
    const double lambda=c0/freq;
    const double near_term=5.0*lambda;
    const double l=40.0*cm;
    const double w=25.0*cm;
    double theta_i=deg2rad(+30.0);
    double phi_i=deg2rad(+60.0);
    cmplx E_TM, E_TE;
    const int Ns=201;
    const double x_min=-100.0*cm;
    const double x_max=+100.0*cm;
    const double y=0.0*cm;
    const double z=+50.0*cm;

    Range x;
    x.linspace(x_min, x_max, Ns);

    //
    Timer T;
    Engine engine(lambda, lc, near_term);
    engine.set_quadl(N_quadl, tol, k_max);
    engine.create_sheet(l, w);
    engine.mesh(FALSE);
    engine.load_mesh(FALSE);
    //

    T.set();
    engine.compute_Z_mn();
    T.unset();
    // engine.save_Z_mn("data/misc/Z_mn.dat");
    // engine.load_Z_mn("data/misc/Z_mn.dat");

    //
    E_TM = +1.0;
    E_TE = +0.0;
    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);
    engine.compute_I_n();
    engine.export_currents();

    //
    Scattered_Field field;
    File file1, file2;
    file1.open("data/sheet/near_field/data1.dat", "w");
    file2.open("data/sheet/near_field/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing near field...");
        Vector r=Vector(x(i), y, z);
        field = engine.compute_near_field(r);
        file1.write("%21.14E %21.14E %21.14E %21.14E\n", 
            x(i)/cm, abs(field.E.x), abs(field.E.y), abs(field.E.z));
        file2.write("%21.14E %21.14E %21.14E %21.14E\n", 
            x(i)/cm, abs(field.H.x)/1.0E-3, abs(field.H.y)/1.0E-3, abs(field.H.z)/1.0E-3);
    }
    file1.close();
    file2.close();
}