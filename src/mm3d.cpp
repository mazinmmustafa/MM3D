//
#include "mm3d.hpp"

#define DEBUG TRUE

// Engine class

Engine::Engine(const double lambda, const double lc, const double near_term){
    check_error(lambda<=0.0, "invalid value of lambda!");
    check_error(lc<=0.0, "invalid value for lc!");
    check_error(near_term<0.0, "invalid value for near_term!");
    this->lambda = lambda;
    this->lc = lc;
    this->near_term = near_term;
}

Engine::~Engine(){
    this->is_quadl_allocated = FALSE;
    this->is_shape_allocated = FALSE;
    this->is_Z_mn_allocated = FALSE;
    this->is_V_m_allocated = FALSE;
    this->is_I_n_allocated = FALSE;
}

void Engine::set_quadl(const int N_quadl, const double tol, const int k_max){
    this->quadl.set(N_quadl, tol, k_max);
    this->is_quadl_allocated = TRUE;
}

void Engine::mesh(const int is_log){
    this->shape.mesh(is_log);
}

void Engine::load_mesh(const int is_log){
    this->shape.set(this->lambda);
    this->shape.load_mesh();
    this->N_basis = this->shape.get_N_basis();
    this->basis_list = this->shape.get_basis_list();
    assert(this->basis_list!=NULL);
    if (is_log){
        this->shape.log("mesh/shape_log.txt");
    }
    printf("%d basis functions were loaded!\n", this->N_basis);
    this->is_shape_allocated = TRUE;
}

void Engine::test_Z_mn(const int m, const int n){
    check_error(!this->is_quadl_allocated, "no quadrature is set!");
    Basis B_m, B_n;
    int flag;
    B_m = this->basis_list[m];
    B_n = this->basis_list[n];
    flag = 0;
    cmplx ans = compute_Zmn(&B_m, &B_n, this->near_term, &this->quadl, flag, FALSE);
    disp(ans);
    if (flag){
        printf("Warning: no convergence at (%d, %d)\n", m, n);
    }
}


void Engine::compute_Z_mn(){
    check_error(!this->is_quadl_allocated, "no quadrature is set!");
    if (this->is_Z_mn_allocated){
        this->Z_mn.deallocate();
    }
    int N=this->N_basis;
    Basis B_m, B_n;
    int flag;
    this->Z_mn.allocate(N, N);
    unsigned long counter=0;
    const int MAX=100;
    char *msg=(char*)calloc(MAX, sizeof(char));
    //
    Timer T;
    T.set();
    printf("estimating Z_mn computation time...\n");
    fflush(stdout);
    for (int n=0; n<1; n++){
        B_n = this->basis_list[n];
        for (int m=n; m<N; m++){
            B_m = this->basis_list[m];
            flag = 0;
            this->Z_mn(m, n) = compute_Zmn(&B_m, &B_n, this->near_term, &this->quadl, flag, FALSE);
            if (flag){
                printf("Warning: no convergence!\n");
            }
            counter++;
        }
        for (int k=n+1; k<N; k++){
            this->Z_mn(n, k) = this->Z_mn(k, n);
            counter++;
        }
    }
    T.unset_silent();
    double time_est=T.get_elapsed();
    time_est*=(double)(N*N)/(2.0*N-1.0);
    printf("estimated time is %0.2f hours\n", time_est/(3600.0)); 
    for (int n=1; n<N; n++){
        B_n = this->basis_list[n];
        for (int m=n; m<N; m++){
            sprintf(msg, "computing Z_mn term (%d, %d)...", m, n);
            progress_bar(counter, (unsigned long)N*(unsigned long)N, msg);
            B_m = this->basis_list[m];
            flag = 0;
            this->Z_mn(m, n) = compute_Zmn(&B_m, &B_n, this->near_term, &this->quadl, flag, FALSE);
            if (flag){
                printf("Warning: no convergence!\n");
            }
            counter++;
        }
        for (int k=n+1; k<N; k++){
            this->Z_mn(n, k) = this->Z_mn(k, n);
            counter++;
        }
    }
    free(msg);
    this->quadl.set(4, 1.0E-3, 15);
    this->is_Z_mn_allocated = TRUE;
}

void Engine::save_Z_mn(const char *filename){
    this->Z_mn.save(filename);
}

void Engine::load_Z_mn(const char *filename){
    this->Z_mn.allocate(this->N_basis, this->N_basis);
    printf("loading Z_mn matrix...");
    this->Z_mn.load(filename);
    printf(", done!\n");
    this->is_Z_mn_allocated = TRUE;
}

void Engine::save_V_m(const char *filename){
    this->V_m.save(filename);
}

void Engine::load_V_m(const char *filename){
    this->V_m.allocate(this->N_basis, this->N_basis);
    printf("loading V_m matrix...");
    this->V_m.load(filename);
    printf(", done!\n");
    this->is_V_m_allocated = TRUE;
}

void Engine::set_Zmn_element(const int m, const int n, const cmplx Z){
    this->Z_mn(m, n) = Z;
}

cmplx Engine::get_Zmn_element(const int m, const int n){
    return this->Z_mn(m, n);
}

void Engine::create_sphere(const double radius){
    check_error(radius<=0.0, "invalid sphere radius!");
    Engine::shape.set(this->lambda);
    Engine::shape.create_sphere(radius, this->lc);
}

void Engine::create_sheet(const double l, const double w){
    check_error(l<=0.0||w<=0.0, "invalid sheet dimensions!");
    Engine::shape.set(this->lambda);
    Engine::shape.create_sheet(l, w, this->lc);
}

void Engine::compute_V_m_plane_wave(const cmplx E_TM, const cmplx E_TE,
    const double theta_i, const double phi_i){
    check_error(!this->is_shape_allocated, "there is no shape defined yet!");
    if (this->is_V_m_allocated){
        this->V_m.deallocate();
    }
    const int N=this->N_basis;
    this->V_m.allocate(N, 1);
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    Vector ki=Vector(k*sin(theta_i)*cos(phi_i), 
                     k*sin(theta_i)*sin(phi_i),
                     k*cos(theta_i));
    Vector Theta_i=Vector(+cos(theta_i)*cos(phi_i),
                          +cos(theta_i)*sin(phi_i),
                          -sin(theta_i));
    Vector Phi_i=Vector(-sin(phi_i), +cos(phi_i), 0.0);
    cmplx I1, I2;
    cmplx theta_pp, theta_pm, theta_mp, theta_mm;
    double Lm;
    Vector x_m, x_p;
    Vector L_pp, L_pm, L_mp, L_mm;
    Vector r_m, r_p;
    Basis Bm;
    cmplx factor_p, factor_m;
    fflush(stdout);
    printf("computing V_m matrix...");
    for (int m=0; m<N; m++){
        Bm = this->basis_list[m];
        r_p = Bm.r_p;
        r_m = Bm.r_m;
        L_pp = Bm.L_pp;
        L_pm = Bm.L_pm;
        L_mp = Bm.L_mp;
        L_mm = Bm.L_mm;
        Lm = Bm.L;
        theta_pp = L_pp*ki;
        theta_pm = L_pm*ki;
        theta_mp = L_mp*ki;
        theta_mm = L_mm*ki;
        factor_p = Lm*exp(+j*ki*r_p);
        factor_m = Lm*exp(+j*ki*r_m);
        //
        I_radiation(-theta_pp, -theta_mp, I1, I2);
        x_p = factor_p*(I1*L_pp+I2*L_mp);
        I_radiation(+theta_mm, +theta_pm, I1, I2);
        x_m = factor_m*(I1*L_mm+I2*L_pm);
        this->V_m(m, 0) = E_TM*(x_p+x_m)*Theta_i+
                          E_TE*(x_p+x_m)*Phi_i;
    }
    printf(", done!\n");
    this->is_V_m_allocated = TRUE;
}

void Engine::compute_I_n(){
    check_error(!this->is_shape_allocated, "there is no shape defined yet!");
    check_error(!this->is_Z_mn_allocated, "Z_mn was not calculated yet!");
    check_error(!this->is_V_m_allocated, "V_m was not calculated yet!");
    if (this->is_I_n_allocated){
        this->I_n.deallocate();
    }
    this->I_n.allocate(this->N_basis, 1);
    printf("solving I_n matrix...");
    this->Z_mn.solve(this->V_m, this->I_n);
    printf(", done!\n");
    this->is_I_n_allocated = TRUE;
}

RCS Engine::compute_rcs(const double theta_s, const double phi_s){
    check_error(!this->is_I_n_allocated, "I_n was not calculated yet!");
    const int N=this->N_basis;
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    const double eta=eta0;
    Vector ks=Vector(k*sin(theta_s)*cos(phi_s), 
                     k*sin(theta_s)*sin(phi_s),
                     k*cos(theta_s));
    Vector Theta_s=Vector(+cos(theta_s)*cos(phi_s),
                          +cos(theta_s)*sin(phi_s),
                          -sin(theta_s));
    Vector Phi_s=Vector(-sin(phi_s), +cos(phi_s), 0.0);
    cmplx I1, I2;
    cmplx theta_pp, theta_pm, theta_mp, theta_mm;
    double Lm;
    Vector x_m, x_p;
    Vector L_pp, L_pm, L_mp, L_mm;
    Vector r_m, r_p;
    Basis Bn;
    cmplx factor_p, factor_m;
    cmplx sum_theta = 0.0;
    cmplx sum_phi = 0.0;
    for (int n=0; n<N; n++){
        Bn = this->basis_list[n];
        r_p = Bn.r_p;
        r_m = Bn.r_m;
        L_pp = Bn.L_pp;
        L_pm = Bn.L_pm;
        L_mp = Bn.L_mp;
        L_mm = Bn.L_mm;
        Lm = Bn.L;
        theta_pp = L_pp*ks;
        theta_pm = L_pm*ks;
        theta_mp = L_mp*ks;
        theta_mm = L_mm*ks;
        factor_p = Lm*exp(+j*ks*r_p);
        factor_m = Lm*exp(+j*ks*r_m);
        //
        I_radiation(-theta_pp, -theta_mp, I1, I2);
        x_p = factor_p*(I1*L_pp+I2*L_mp);
        I_radiation(+theta_mm, +theta_pm, I1, I2);
        x_m = factor_m*(I1*L_mm+I2*L_pm);
        //
        sum_theta+=(eta/2.0)*(x_p+x_m)*Theta_s*this->I_n(n, 0);
        sum_phi+=(eta/2.0)*(x_p+x_m)*Phi_s*this->I_n(n, 0);
    }
    RCS sigma;
    sigma.theta = 4.0*pi*pow(abs(sum_theta), 2.0);
    sigma.phi = 4.0*pi*pow(abs(sum_phi), 2.0);
    return sigma;
}

struct Current_Info{
    Vector v;
    cmplx I=0.0;
    char s='-';
    double L=0.0, A=0.0;
};

void sub_divide_triangles(Vector v1, Vector v2, Vector v3,
    Current_Info info[3], int k, int k_max, Vector rp){
    const double alpha=0.0, beta=1.0/2.0;
    k++;
    if (k<k_max){
        Vector vc=v1+alpha*(v2-v1)+beta*(v3-v1);
        Vector v1_1=v2;
        Vector v2_1=vc;
        Vector v3_1=v1;
        Vector v1_2=v3;
        Vector v2_2=vc;
        Vector v3_2=v2;
        sub_divide_triangles(v1_1, v2_1, v3_1, info, k, k_max, v1_1);
        sub_divide_triangles(v1_2, v2_2, v3_2, info, k, k_max, v1_2);
    }
    if (k==k_max){
        Vector rp_new = rp+(1.0/3.0)*(v2-v1)+(1.0/3.0)*(v3-v1);
        File file;
        file.open("mesh/current.pos", "a");
        Vector J;
        for (int i=0; i<3; i++){
            if (info[i].s=='-'){
                J = J+info[i].I*(rp_new-info[i].v)*info[i].L/(2.0*info[i].A);
            }
            if (info[i].s=='+'){
                J = J+info[i].I*(info[i].v-rp_new)*info[i].L/(2.0*info[i].A);
            }
        }
        double c=J.mag();
        file.write("ST(%11.4E, %11.4E, %11.4E, %11.4E, %11.4E, %11.4E, %11.4E, %11.4E, %11.4E)",
            real(v1.x), real(v1.y), real(v1.z), 
            real(v2.x), real(v2.y), real(v2.z), 
            real(v3.x), real(v3.y), real(v3.z));
        file.write("{%11.4E, %11.4E, %11.4E};\n", c, c, c);
        file.close();
    }
}

void Engine::export_currents(){
    check_error(!this->is_I_n_allocated, "I_n was not calculated yet!");
    File file;
    file.open("mesh/current.pos", "w");
    file.write("View \"background mesh\" {\n");
    file.close();
    int N=this->N_basis;
    Basis B_m, B_n;
    Current_Info info[3];
    const int k_max=6;
    fflush(stdout);
    printf("computing currents...");
    List<Triangle> list;
    for (int m=0; m<N; m++){
        B_m = this->basis_list[m];
        int counter=0;
        info[0].I = 0.0; info[1].I = 0.0; info[2].I = 0.0;
        info[counter].I = this->I_n(m, 0);
        info[counter].v = B_m.r_p;
        info[counter].L = B_m.L;
        info[counter].A = B_m.A_p;
        info[counter].s = '+';
        counter++;
        for (int n=0; n<N; n++){
            B_n = this->basis_list[n];
            // Scenario II case I
            if ((B_m.r_p.is_equal(B_n.e_m))&&
                (B_m.e_p.is_equal(B_n.r_p))&&
                (B_m.e_m.is_equal(B_n.e_p))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_p;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_p;
                info[counter].s = '+';
                counter++;
            }else
            // Scenario II case III
            if ((B_m.r_p.is_equal(B_n.e_p))&&
                (B_m.e_p.is_equal(B_n.e_m))&&
                (B_m.e_m.is_equal(B_n.r_p))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_p;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_p;
                info[counter].s = '+';
                counter++;
            }else
            // Scenario II case II
            if ((B_m.r_p.is_equal(B_n.e_p))&&
                (B_m.e_p.is_equal(B_n.r_m))&&
                (B_m.e_m.is_equal(B_n.e_m))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_m;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_m;
                info[counter].s = '-';
                counter++;
            }else
            // Scenario II case IV
            if ((B_m.r_p.is_equal(B_n.e_m))&&
                (B_m.e_p.is_equal(B_n.e_p))&&
                (B_m.e_m.is_equal(B_n.r_m))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_m;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_m;
                info[counter].s = '-';
                counter++;
            }
        }
        assert(counter<=3);
        int flag=FALSE;
        for (int i=0; i<list.len(); i++){
            if (B_m.r_p.is_equal(list(i).v1)&&
                B_m.e_p.is_equal(list(i).v2)&&
                B_m.e_m.is_equal(list(i).v3)){
                flag = TRUE;
            }
            if (B_m.r_p.is_equal(list(i).v3)&&
                B_m.e_p.is_equal(list(i).v1)&&
                B_m.e_m.is_equal(list(i).v2)){
                flag = TRUE;
            }
            if (B_m.r_p.is_equal(list(i).v2)&&
                B_m.e_p.is_equal(list(i).v3)&&
                B_m.e_m.is_equal(list(i).v1)){
                flag = TRUE;
            }
        }
        if (!flag){
            sub_divide_triangles(B_m.r_p, B_m.e_p, B_m.e_m, info, 0, k_max, B_m.r_p);
        }
        Triangle T;
        T.v1=B_m.r_p; T.v2=B_m.e_p; T.v3=B_m.e_m;
        list.append(T);
    }
    for (int m=0; m<N; m++){
        B_m = this->basis_list[m];
        int counter=0;
        info[0].I = 0.0; info[1].I = 0.0; info[2].I = 0.0;
        info[counter].I = this->I_n(m, 0);
        info[counter].v = B_m.r_m;
        info[counter].L = B_m.L;
        info[counter].A = B_m.A_m;
        info[counter].s = '-';
        counter++;
        for (int n=0; n<N; n++){
            B_n = this->basis_list[n];
            // Scenario III case I
            if ((B_m.r_m.is_equal(B_n.e_p))&&
                (B_m.e_m.is_equal(B_n.e_m))&&
                (B_m.e_p.is_equal(B_n.r_p))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_p;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_p;
                info[counter].s = '+';
                counter++;
            }else
            // Scenario III case III
            if ((B_m.r_m.is_equal(B_n.e_m))&&
                (B_m.e_m.is_equal(B_n.r_p))&&
                (B_m.e_p.is_equal(B_n.e_p))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_p;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_p;
                info[counter].s = '+';
                counter++;
            }else
            // Scenario III case II
            if ((B_m.r_m.is_equal(B_n.e_m))&&
                (B_m.e_m.is_equal(B_n.e_p))&&
                (B_m.e_p.is_equal(B_n.r_m))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_m;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_m;
                info[counter].s = '-';
                counter++;
            }else
            // Scenario III case IV
            if ((B_m.r_m.is_equal(B_n.e_p))&&
                (B_m.e_m.is_equal(B_n.r_m))&&
                (B_m.e_p.is_equal(B_n.e_m))){
                info[counter].I = this->I_n(n, 0);
                info[counter].v = B_n.r_m;
                info[counter].L = B_n.L;
                info[counter].A = B_n.A_m;
                info[counter].s = '-';
                counter++;
            }
        }
        assert(counter<=3);
        int flag=FALSE;
        for (int i=0; i<list.len(); i++){
            if (B_m.r_m.is_equal(list(i).v1)&&
                B_m.e_m.is_equal(list(i).v2)&&
                B_m.e_p.is_equal(list(i).v3)){
                flag = TRUE;
            }
            if (B_m.r_m.is_equal(list(i).v3)&&
                B_m.e_m.is_equal(list(i).v1)&&
                B_m.e_p.is_equal(list(i).v2)){
                flag = TRUE;
            }
            if (B_m.r_m.is_equal(list(i).v2)&&
                B_m.e_m.is_equal(list(i).v3)&&
                B_m.e_p.is_equal(list(i).v1)){
                flag = TRUE;
            }
        }
        if (!flag){
            sub_divide_triangles(B_m.r_m, B_m.e_m, B_m.e_p, info, 0, k_max, B_m.r_m);
        }
        Triangle T;
        T.v1=B_m.r_m; T.v2=B_m.e_m; T.v3=B_m.e_p;
        list.append(T);
    }
    file.open("mesh/current.pos", "a");
    file.write("};\n");
    file.close();
    printf(", done!\n");
}


struct Integrand_Args{
    Engine *engine=NULL;
    double k=2.0*pi;
    Vector unit_vector;
    Vector r;
};

cmplx near_field_E_integrand(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    const cmplx j=cmplx(0.0, 1.0);
    const double k=args->k;
    const double eta=eta0;
    Engine *engine=args->engine;
    Basis B_m;
    Vector rho_m, rho_p;
    int N=engine->N_basis;
    Vector sum_E;
    Vector R_m, R_p;
    Vector r=args->r;
    double Lm;
    cmplx g_m, g_p;
    cmplx I;
    Vector A, B, C, D;
    for (int m=0; m<N; m++){
        B_m = engine->basis_list[m];
        rho_m = alpha*B_m.L_mm+beta*B_m.L_pm;
        rho_p = alpha*B_m.L_pp+beta*B_m.L_mp;
        R_m = r-(B_m.r_m+rho_m);
        R_p = r-(B_m.r_p-rho_p);
        Lm = B_m.L;
        g_m = exp(-j*k*R_m.mag())/(4.0*pi*R_m.mag());
        g_p = exp(-j*k*R_p.mag())/(4.0*pi*R_p.mag());
        I = engine->I_n(m, 0);
        A = (rho_m*g_m);
        B = (rho_p*g_p);
        C = 2.0*((1.0+j*k*R_m.mag())/(R_m.mag()))*g_m*R_m.unit();
        D = 2.0*((1.0+j*k*R_p.mag())/(R_p.mag()))*g_p*R_p.unit();
        sum_E = sum_E+(-j*k*eta*(A+B)+(j*eta/k)*(C-D))*I*Lm;
    }
    return sum_E*args->unit_vector;
}

cmplx near_field_H_integrand(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    const cmplx j=cmplx(0.0, 1.0);
    const double k=args->k;
    Engine *engine=args->engine;
    Basis B_m;
    Vector rho_m, rho_p;
    int N=engine->N_basis;
    Vector sum_H;
    Vector R_m, R_p;
    Vector r=args->r;
    double Lm;
    cmplx g_m, g_p;
    cmplx I;
    Vector A, B;
    for (int m=0; m<N; m++){
        B_m = engine->basis_list[m];
        rho_m = alpha*B_m.L_mm+beta*B_m.L_pm;
        rho_p = alpha*B_m.L_pp+beta*B_m.L_mp;
        R_m = r-(B_m.r_m+rho_m);
        R_p = r-(B_m.r_p-rho_p);
        Lm = B_m.L;
        g_m = exp(-j*k*R_m.mag())/(4.0*pi*R_m.mag());
        g_p = exp(-j*k*R_p.mag())/(4.0*pi*R_p.mag());
        I = engine->I_n(m, 0);
        A = ((1.0+j*k*R_m.mag())/(R_m.mag()))*g_m*(rho_m^R_m.unit());
        B = ((1.0+j*k*R_p.mag())/(R_p.mag()))*g_p*(rho_p^R_p.unit());
        sum_H = sum_H+(A+B)*I*Lm;
    }
    return sum_H*args->unit_vector;
}

Scattered_Field Engine::compute_near_field(Vector position){
    check_error(!this->is_I_n_allocated, "I_n was not calculated yet!");
    Vector r=position/this->lambda;
    Integrand_Args args;
    args.k = 2.0*pi;
    args.r = r;
    args.engine = this;
    int flag;
    Scattered_Field field;
    Triangle_2D T;
    T.v1 = Vector(+0.0, +0.0, +0.0);
    T.v2 = Vector(+1.0, +0.0, +0.0);
    T.v3 = Vector(+0.0, +1.0, +0.0);
    //
    args.unit_vector = Vector(1.0, 0.0, 0.0);
    field.E.x = this->quadl.integral_2D_tri(near_field_E_integrand, &args, T, flag);
    field.H.x = this->quadl.integral_2D_tri(near_field_H_integrand, &args, T, flag);
    if (flag){ printf("warning: no convergence!\n"); }
    args.unit_vector = Vector(0.0, 1.0, 0.0);
    field.E.y = this->quadl.integral_2D_tri(near_field_E_integrand, &args, T, flag);
    field.H.y = this->quadl.integral_2D_tri(near_field_H_integrand, &args, T, flag);
    if (flag){ printf("warning: no convergence!\n"); }
    args.unit_vector = Vector(0.0, 0.0, 1.0);
    field.E.z = this->quadl.integral_2D_tri(near_field_E_integrand, &args, T, flag);
    field.H.z = this->quadl.integral_2D_tri(near_field_H_integrand, &args, T, flag);
    if (flag){ printf("warning: no convergence!\n"); }
    //
    return field;
}