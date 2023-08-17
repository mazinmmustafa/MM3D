//
#include "engine.hpp"

#define DEBUG FALSE

struct Integrand_Args{
    Basis *B_m=NULL;
    Basis *B_n=NULL;
    double k=2.0*pi;
    char s_m, s_n;
    Singular_Para para_I_p;
    Singular_Para para_I_m;
    Singular_Para para_II;
    cmplx s1, s2;
};

double R_mn(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, const Basis *B_m,
    const Basis *B_n, const char s_m, const char s_n){
    int flag=0;
    Vector R;
    if (s_m=='+'&&s_n=='+'){
        R = (B_m->r_p)-(B_n->r_p)-alpha*(B_m->L_pp)+alpha_*(B_n->L_pp)-beta*(B_m->L_mp)+beta_*(B_n->L_mp);
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        R = (B_m->r_p)-(B_n->r_m)-alpha*(B_m->L_pp)-alpha_*(B_n->L_mm)-beta*(B_m->L_mp)-beta_*(B_n->L_pm);
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        R = (B_m->r_m)-(B_n->r_p)+alpha*(B_m->L_mm)+alpha_*(B_n->L_pp)+beta*(B_m->L_pm)+beta_*(B_n->L_mp);
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        R = (B_m->r_m)-(B_n->r_m)+alpha*(B_m->L_mm)-alpha_*(B_n->L_mm)+beta*(B_m->L_pm)-beta_*(B_n->L_pm);
        flag++;
    }
    assert(flag==1);
    return R.mag();
}

cmplx g_mn(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, const Basis *B_m,
    const Basis *B_n, const double k, const char s_m, const char s_n){
    const cmplx j=cmplx(0.0, 1.0);
    double R=R_mn(alpha, beta, alpha_, beta_, B_m, B_n, s_m, s_n);
    return exp(-j*k*R)/(4.0*pi*R);
}

cmplx psi_IV_far_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    const Basis *B_m=args->B_m;
    const Basis *B_n=args->B_n;
    const double k=args->k;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx g;
    Vector rho_m, rho_n;
    if (s_m=='+'&&s_n=='+'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '+', '+');
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        rho_n = alpha_*B_n->L_pp+beta_*B_n->L_mp;
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '+', '-');
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        rho_n = alpha_*B_n->L_mm+beta_*B_n->L_pm;
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '-', '+');
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        rho_n = alpha_*B_n->L_pp+beta_*B_n->L_mp;
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '-', '-');
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        rho_n = alpha_*B_n->L_mm+beta_*B_n->L_pm;
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return L_m*L_n*(rho_m*rho_n)*g;
}

cmplx phi_IV_far_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    const Basis *B_m=args->B_m;
    const Basis *B_n=args->B_n;
    const double k=args->k;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx g;
    if (s_m=='+'&&s_n=='+'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '+', '+');
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '+', '-');
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '-', '+');
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        g = g_mn(alpha, beta, alpha_, beta_, B_m, B_n, k, '-', '-');
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return 4.0*L_m*L_n*g;
}

cmplx phi_IV_near_integrand_I(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const Basis *B_m=args->B_m;
    const Basis *B_n=args->B_n;
    const double k=args->k;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    double R;
    if (s_m=='+'&&s_n=='+'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '+', '+');
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '+', '-');
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '-', '+');
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '-', '-');
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return (-j*k/pi)*L_m*L_n*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx phi_IV_near_integrand_II(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    Basis *B_m=args->B_m;
    Basis *B_n=args->B_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx term;
    double A_n=0.0;
    if (s_m=='+'&&s_n=='+'){
        term = T1_projection(alpha, beta, B_m, B_n, '+', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        term = T1_projection(alpha, beta, B_m, B_n, '+', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        term = T1_projection(alpha, beta, B_m, B_n, '-', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        term = T1_projection(alpha, beta, B_m, B_n, '-', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return ((L_m*L_n)/(2.0*pi*A_n))*term;
}

cmplx psi_IV_near_integrand_I(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const Basis *B_m=args->B_m;
    const Basis *B_n=args->B_n;
    const double k=args->k;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx factor=0.0;
    Vector rho_m, rho_n;
    double R;
    if (s_m=='+'&&s_n=='+'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '+', '+');
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        rho_n = alpha_*B_n->L_pp+beta_*B_n->L_mp;
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '+', '-');
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        rho_n = alpha_*B_n->L_mm+beta_*B_n->L_pm;
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '-', '+');
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        rho_n = alpha_*B_n->L_pp+beta_*B_n->L_mp;
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        R = R_mn(alpha, beta, alpha_, beta_, B_m, B_n, '-', '-');
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        rho_n = alpha_*B_n->L_mm+beta_*B_n->L_pm;
        flag++;
    }
    assert(flag==1);
    factor = rho_m*rho_n;
    double L_m=B_m->L;
    double L_n=B_n->L;
    return (-j*k/(4.0*pi))*L_m*L_n*factor*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx psi_IV_near_integrand_II(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    Basis *B_m=args->B_m;
    Basis *B_n=args->B_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx term;
    double A_n=0.0;
    cmplx factor=0.0;
    Vector rho_m, rho_n;
    Projection_3D_Para para;
    Vector p;
    Triangle T;
    if (s_m=='+'&&s_n=='+'){
        term = T1_projection(alpha, beta, B_m, B_n, '+', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        T.v1 = B_n->r_p; T.v2 = B_n->e_p; T.v3 = B_n->e_m;
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        p = B_m->r_p-rho_m;
        projection_3D(T, p, para);
        rho_n = B_n->r_p-para.p0;
        factor = +1.0*(rho_m*rho_n);
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        term = T1_projection(alpha, beta, B_m, B_n, '+', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        T.v1 = B_n->r_m; T.v2 = B_n->e_m; T.v3 = B_n->e_p;
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        p = B_m->r_p-rho_m;
        projection_3D(T, p, para);
        rho_n = B_n->r_m-para.p0;
        factor = -1.0*(rho_m*rho_n);
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        term = T1_projection(alpha, beta, B_m, B_n, '-', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        T.v1 = B_n->r_p; T.v2 = B_n->e_p; T.v3 = B_n->e_m;
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        p = B_m->r_m+rho_m;
        projection_3D(T, p, para);
        rho_n = B_n->r_p-para.p0;
        factor = +1.0*(rho_m*rho_n);
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        term = T1_projection(alpha, beta, B_m, B_n, '-', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        T.v1 = B_n->r_m; T.v2 = B_n->e_m; T.v3 = B_n->e_p;
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        p = B_m->r_m+rho_m;
        projection_3D(T, p, para);
        rho_n = B_n->r_m-para.p0;
        factor = -1.0*(rho_m*rho_n);
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return ((L_m*L_n)/(8.0*pi*A_n))*factor*term;
}

cmplx psi_IV_near_integrand_III(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    Basis *B_m=args->B_m;
    Basis *B_n=args->B_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    int flag=0;
    cmplx term;
    double A_n=0.0;
    Vector rho_m;
    if (s_m=='+'&&s_n=='+'){
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        term = -1.0*rho_m*T2_projection(alpha, beta, B_m, B_n, '+', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        flag++;
    }
    if (s_m=='+'&&s_n=='-'){
        rho_m = alpha*B_m->L_pp+beta*B_m->L_mp;
        term = +1.0*rho_m*T2_projection(alpha, beta, B_m, B_n, '+', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        flag++;
    }
    if (s_m=='-'&&s_n=='+'){
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        term = -1.0*rho_m*T2_projection(alpha, beta, B_m, B_n, '-', '+'); 
        A_n = get_area(B_n->r_p, B_n->e_p, B_n->e_m);
        flag++;
    }
    if (s_m=='-'&&s_n=='-'){
        rho_m = alpha*B_m->L_mm+beta*B_m->L_pm;
        term = +1.0*rho_m*T2_projection(alpha, beta, B_m, B_n, '-', '-'); 
        A_n = get_area(B_n->r_m, B_n->e_m, B_n->e_p);
        flag++;
    }
    assert(flag==1);
    double L_m=B_m->L;
    double L_n=B_n->L;
    return ((L_m*L_n)/(8.0*pi*A_n))*term;
}

cmplx I1_singular_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const double k=args->k;
    //
    Vector B, C;
    double a;
    int flag=0;
    if (args->s_m=='+'){
        B = args->para_I_p.B;
        C = args->para_I_p.C;
        a = args->para_I_p.a;
        flag++;
    }
    if (args->s_m=='-'){
        B = args->para_I_m.B;
        C = args->para_I_m.C;
        a = args->para_I_m.a;
        flag++;
    }
    assert(flag==1);
    //
    double R;
    R = ((alpha-alpha_)*B+(beta-beta_)*C).mag();
    return (-j*k*a*a/(pi))*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I1_singular(const Singular_Para *para){
    double p;
    //
    double a=para->a;
    double b=para->b;
    double c=para->c;
    //
    p = 0.5*(a+b+c);
    return -(a*a/(3.0*pi))*(log(1.0-a/p)/a+
                            log(1.0-b/p)/b+
                            log(1.0-c/p)/c);
}

cmplx I2_singular_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const double k=args->k;
    Vector B, C;
    double a, b, c;
    int flag=0;
    if (args->s_m=='+'){
        B = args->para_I_p.B;
        C = args->para_I_p.C;
        a = args->para_I_p.a;
        b = args->para_I_p.b;
        c = args->para_I_p.c;
        flag++;
    }
    if (args->s_m=='-'){
        B = args->para_I_m.B;
        C = args->para_I_m.C;
        a = args->para_I_m.a;
        b = args->para_I_m.b;
        c = args->para_I_m.c;
        flag++;
    }
    assert(flag==1);
    double R;
    R = ((alpha-alpha_)*B+(beta-beta_)*C).mag();
    cmplx factor=alpha*alpha_*b*b+beta*beta_*c*c+(alpha*beta_+beta*alpha_)*(B*C);
    return (-j*k*a*a/(4.0*pi))*factor*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I2_singular(const Singular_Para *para){
    double p;
    //
    double Am=para->Am;
    double a=para->a;
    double b=para->b;
    double c=para->c;
    //
    p = 0.5*(a+b+c);
    double term=(10.0+3.0*(c*c-a*a)/(b*b)-3.0*(a*a-b*b)/(c*c))*a
               -(5.0-3.0*(a*a-b*b)/(c*c)-2.0*(b*b-c*c)/(a*a))*b
               -(5.0+3.0*(c*c-a*a)/(b*b)+2.0*(b*b-c*c)/(a*a))*c
               +(a*a-3.0*b*b-3.0*c*c-8.0*Am*Am/(a*a))*(2.0/a)*log(1.0-a/p)
               +(a*a-2.0*b*b-4.0*c*c+6.0*Am*Am/(b*b))*(4.0/b)*log(1.0-b/p)
               +(a*a-4.0*b*b-2.0*c*c+6.0*Am*Am/(c*c))*(4.0/c)*log(1.0-c/p);
    return +(a*a/(480.0*pi))*term;
}

cmplx I3_singular(const Singular_Para *para){
    double p;
    //
    double a=para->a;
    double b=para->b;
    double c=para->c;
    //
    p = 0.5*(a+b+c);
    return -(b*c/(3.0*pi))*(log(1.0-a/p)/a+
                            log(1.0-b/p)/b+
                            log(1.0-c/p)/c);
}

cmplx I3_singular_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const double k=args->k;
    //
    Vector A=args->para_II.A;
    Vector B=args->para_II.B;
    Vector C=args->para_II.C;
    double b=args->para_II.b;
    double c=args->para_II.c;
    //
    double R;
    R = ((alpha+alpha_-1.0)*A+beta*B-beta_*C).mag();
    return (-j*k*b*c/(pi))*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I4_singular_integrand(const cmplx alpha, const cmplx beta,
    const cmplx alpha_, const cmplx beta_, void *args_in){
    const cmplx j=cmplx(0.0, 1.0);
    Integrand_Args *args=(Integrand_Args*)args_in;
    const double k=args->k;
    //
    Vector A=args->para_II.A;
    Vector B=args->para_II.B;
    Vector C=args->para_II.C;
    double a=args->para_II.a;
    double b=args->para_II.b;
    double c=args->para_II.c;
    double sign=args->para_II.sign;
    //
    double R;
    R = ((alpha+alpha_-1.0)*A+beta*B-beta_*C).mag();
    cmplx factor=beta*beta_*(B*C)-alpha*alpha_*a*a+A*(alpha*beta_*C-beta*alpha_*B);
    return sign*(-j*k*b*c/(4.0*pi))*factor*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I4_singular(const Singular_Para *para){
    double p;
    //
    Vector A=para->A;
    Vector B=para->B;
    Vector C=para->C;
    double Am=para->Am;
    double a=para->a;
    double b=para->b;
    double c=para->c;
    double sign=para->sign;
    //
    p = 0.5*(a+b+c);
    double term=(-10.0+(c*c-a*a)/(b*b)-(a*a-b*b)/(c*c))*a
               +(5.0+(a*a-b*b)/(c*c)-6.0*(b*b-c*c)/(a*a))*b
               +(5.0-(c*c-a*a)/(b*b)+6.0*(b*b-c*c)/(a*a))*c
               +(2.0*a*a-b*b-c*c+4.0*Am*Am/(a*a))*(12.0/a)*log(1.0-a/p)
               +(9.0*a*a-3.0*b*b-c*c+4.0*Am*Am/(b*b))*(2.0/b)*log(1.0-b/p)
               +(9.0*a*a-b*b-3.0*c*c+4.0*Am*Am/(c*c))*(2.0/c)*log(1.0-c/p);
    return +sign*(b*c/(960.0*pi))*term;
}

int is_near_term_check(const Triangle T_m, const Triangle T_n, const double near_term){
    double dist[9];
    dist[0] = (T_m.v1-T_n.v1).mag();
    dist[1] = (T_m.v1-T_n.v2).mag();
    dist[2] = (T_m.v1-T_n.v3).mag();
    dist[3] = (T_m.v2-T_n.v1).mag();
    dist[4] = (T_m.v2-T_n.v2).mag();
    dist[5] = (T_m.v2-T_n.v3).mag();
    dist[6] = (T_m.v3-T_n.v1).mag();
    dist[7] = (T_m.v3-T_n.v2).mag();
    dist[8] = (T_m.v3-T_n.v3).mag();
    double max_dist=dist[0];
    for (int i=1; i<9; i++){
        if (dist[i]>max_dist){
            max_dist = dist[i];
        }
    }
    return max_dist<=near_term ? TRUE : FALSE;
}

int is_near_term(Basis *B_m, Basis *B_n, const double near_term){
    int flag=0;
    Triangle T_m_p; T_m_p.v1 = B_m->r_p; T_m_p.v2 = B_m->e_p; T_m_p.v3 = B_m->e_m;
    Triangle T_m_m; T_m_m.v1 = B_m->r_m; T_m_m.v2 = B_m->e_m; T_m_m.v3 = B_m->e_p;
    Triangle T_n_p; T_n_p.v1 = B_n->r_p; T_n_p.v2 = B_n->e_p; T_n_p.v3 = B_n->e_m;
    Triangle T_n_m; T_n_m.v1 = B_n->r_m; T_n_m.v2 = B_n->e_m; T_n_m.v3 = B_n->e_p;
    flag = is_near_term_check(T_m_p, T_n_p, near_term); if (flag){return TRUE;}
    flag = is_near_term_check(T_m_p, T_n_m, near_term); if (flag){return TRUE;}
    flag = is_near_term_check(T_m_m, T_n_p, near_term); if (flag){return TRUE;}
    flag = is_near_term_check(T_m_m, T_n_m, near_term); if (flag){return TRUE;}
    return FALSE;
}

cmplx compute_Zmn(Basis *B_m, Basis *B_n, const double near_term, QuadL *quadl, int &flag, 
    const int scenario){
    assert(near_term>=0.0);
    //
    const double k=2.0*pi;
    const double eta=eta0;
    const cmplx j=cmplx(0.0, 1.0); 
    //
    Integrand_Args args;
    args.B_m = B_m;
    args.B_n = B_n;
    args.k = k;
    //
    Triangle_2D T_m;
    T_m.v1 = Vector(+0.0, +0.0, +0.0);
    T_m.v2 = Vector(+1.0, +0.0, +0.0);
    T_m.v3 = Vector(+0.0, +1.0, +0.0);
    Triangle_2D T_n;
    T_n.v1 = Vector(+0.0, +0.0, +0.0);
    T_n.v2 = Vector(+1.0, +0.0, +0.0);
    T_n.v3 = Vector(+0.0, +1.0, +0.0);
    //
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    cmplx A=0.0, B=0.0;
    //
    check_error(scenario<0||scenario>7, "invalid scenario!");
    int scenario_local=0;
    Singular_Para para_I_p, para_I_m, para_II;
    if (scenario==0){
        // Senario I
        if ((B_m->r_p.is_equal(B_n->r_p))&&
            (B_m->e_m.is_equal(B_n->e_m))&&
            (B_m->r_m.is_equal(B_n->r_m))&&
            (B_m->e_p.is_equal(B_n->e_p))){
            //
            para_I_p.A = B_m->e_m-B_m->e_p;
            para_I_p.B = B_m->e_m-B_m->r_p;
            para_I_p.C = B_m->e_p-B_m->r_p;
            para_I_p.a = para_I_p.A.mag();
            para_I_p.b = para_I_p.B.mag();
            para_I_p.c = para_I_p.C.mag();
            para_I_p.Am = B_m->A_p;
            //
            para_I_m.A = B_m->e_p-B_m->e_m;
            para_I_m.B = B_m->e_p-B_m->r_m;
            para_I_m.C = B_m->e_m-B_m->r_m;
            para_I_m.a = para_I_m.A.mag();
            para_I_m.b = para_I_m.B.mag();
            para_I_m.c = para_I_m.C.mag();
            para_I_m.Am = B_m->A_m;
            //
            scenario_local = 1;
            if (DEBUG){
                printf("Scenario I:\n");
            }
            goto label;
        }
        // Scenario II case I
        if ((B_m->r_p.is_equal(B_n->e_m))&&
            (B_m->e_p.is_equal(B_n->r_p))&&
            (B_m->e_m.is_equal(B_n->e_p))){
            para_II.A = B_m->e_p-B_m->r_p;
            para_II.B = B_m->e_m-B_m->r_p;
            para_II.C = B_m->e_m-B_m->e_p;
            para_II.Am = B_m->A_p;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = +1.0;
            scenario_local = 2;
            if (DEBUG){
                printf("Scenario II case I:\n");
            }
            goto label;
        }
        // Scenario II case III
        if ((B_m->r_p.is_equal(B_n->e_p))&&
            (B_m->e_p.is_equal(B_n->e_m))&&
            (B_m->e_m.is_equal(B_n->r_p))){
            para_II.A = B_m->e_m-B_m->r_p;
            para_II.B = B_m->e_p-B_m->r_p;
            para_II.C = B_m->e_p-B_m->e_m;
            para_II.Am = B_m->A_p;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = +1.0;
            scenario_local = 2;
            if (DEBUG){
                printf("Scenario II case III:\n");
            }
            goto label;
        }
        // Scenario II case II
        if ((B_m->r_p.is_equal(B_n->e_p))&&
            (B_m->e_p.is_equal(B_n->r_m))&&
            (B_m->e_m.is_equal(B_n->e_m))){
            para_II.A = B_m->e_p-B_m->r_p;
            para_II.B = B_m->e_m-B_m->r_p;
            para_II.C = B_m->e_m-B_m->e_p;
            para_II.Am = B_m->A_p;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = -1.0;
            scenario_local = 3;
            if (DEBUG){
                printf("Scenario II case II:\n");
            }
            goto label;
        }
        // Scenario II case IV
        if ((B_m->r_p.is_equal(B_n->e_m))&&
            (B_m->e_p.is_equal(B_n->e_p))&&
            (B_m->e_m.is_equal(B_n->r_m))){
            para_II.A = B_m->e_m-B_m->r_p;
            para_II.B = B_m->e_p-B_m->r_p;
            para_II.C = B_m->e_p-B_m->e_m;
            para_II.Am = B_m->A_p;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = -1.0;
            scenario_local =3;
            if (DEBUG){
                printf("Scenario II case IV:\n");
            }
            goto label;
        }
        // Scenario III case I
        if ((B_m->r_m.is_equal(B_n->e_p))&&
            (B_m->e_m.is_equal(B_n->e_m))&&
            (B_m->e_p.is_equal(B_n->r_p))){
            para_II.A = B_m->e_p-B_m->r_m;
            para_II.B = B_m->e_m-B_m->r_m;
            para_II.C = B_m->e_m-B_m->e_p;
            para_II.Am = B_m->A_m;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = -1.0;
            scenario_local = 4;
            if (DEBUG){
                printf("Scenario III case I:\n");
            }
            goto label;
        }
        // Scenario III case III
        if ((B_m->r_m.is_equal(B_n->e_m))&&
            (B_m->e_m.is_equal(B_n->r_p))&&
            (B_m->e_p.is_equal(B_n->e_p))){
            para_II.A = B_m->e_m-B_m->r_m;
            para_II.B = B_m->e_p-B_m->r_m;
            para_II.C = B_m->e_p-B_m->e_m;
            para_II.Am = B_m->A_m;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = -1.0;
            scenario_local = 4;
            if (DEBUG){
                printf("Scenario III case III:\n");
            }
            goto label;
        }
        // Scenario III case II
        if ((B_m->r_m.is_equal(B_n->e_m))&&
            (B_m->e_m.is_equal(B_n->e_p))&&
            (B_m->e_p.is_equal(B_n->r_m))){
            para_II.A = B_m->e_p-B_m->r_m;
            para_II.B = B_m->e_m-B_m->r_m;
            para_II.C = B_m->e_m-B_m->e_p;
            para_II.Am = B_m->A_m;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = +1.0;
            scenario_local = 5;
            if (DEBUG){
                printf("Scenario III case II:\n");
            }
            goto label;
        }
        // Scenario III case IV
        if ((B_m->r_m.is_equal(B_n->e_p))&&
            (B_m->e_m.is_equal(B_n->r_m))&&
            (B_m->e_p.is_equal(B_n->e_m))){
            para_II.A = B_m->e_m-B_m->r_m;
            para_II.B = B_m->e_p-B_m->r_m;
            para_II.C = B_m->e_p-B_m->e_m;
            para_II.Am = B_m->A_m;
            para_II.a = para_II.A.mag();
            para_II.b = para_II.B.mag();
            para_II.c = para_II.C.mag();
            para_II.sign = +1.0;
            scenario_local = 5;
            if (DEBUG){
                printf("Scenario III case IV:\n");
            }
            goto label;
        }
        // scenario IV near
        if (is_near_term(B_m, B_n, near_term)){
            scenario_local = 6;
            if (DEBUG){
                printf("Scenario IV near:\n");
            }
            goto label;
        // scenario IV far    
        }else{
            scenario_local = 7;
            if (DEBUG){
                printf("Scenario IV far:\n");
            }
            goto label;
        }
    }else{
        scenario_local = scenario;
    }
    label: 
    assert(scenario_local>0&&scenario_local<8);
    args.para_I_p = para_I_p;
    args.para_I_m = para_I_m;
    args.para_II = para_II;
    // Senario I
    switch (scenario_local){
        case 1:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(I2_singular_integrand, &args, T_m, T_n, flag)+
                     I2_singular(&para_I_p);
            phi_pp = quadl->integral_4D_tri(I1_singular_integrand, &args, T_m, T_n, flag)+
                     I1_singular(&para_I_p);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(I2_singular_integrand, &args, T_m, T_n, flag)+
                     I2_singular(&para_I_m);
            phi_mm = quadl->integral_4D_tri(I1_singular_integrand, &args, T_m, T_n, flag)+
                     I1_singular(&para_I_m);
            break;
        case 2:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(I4_singular_integrand, &args, T_m, T_n, flag)+
                     I4_singular(&para_II);
            phi_pp = quadl->integral_4D_tri(I3_singular_integrand, &args, T_m, T_n, flag)+
                     I3_singular(&para_II);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            break;
        case 3:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(I4_singular_integrand, &args, T_m, T_n, flag)+
                     I4_singular(&para_II);
            phi_pm = quadl->integral_4D_tri(I3_singular_integrand, &args, T_m, T_n, flag)+
                     I3_singular(&para_II);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            break;
        case 4:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(I4_singular_integrand, &args, T_m, T_n, flag)+
                     I4_singular(&para_II);
            phi_mp = quadl->integral_4D_tri(I3_singular_integrand, &args, T_m, T_n, flag)+
                     I3_singular(&para_II);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            break;
        case 5:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(I4_singular_integrand, &args, T_m, T_n, flag)+
                     I4_singular(&para_II);
            phi_mm = quadl->integral_4D_tri(I3_singular_integrand, &args, T_m, T_n, flag)+
                     I3_singular(&para_II);
            break;
        case 6:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(psi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_II, &args, T_m, flag)+
                     quadl->integral_2D_tri(psi_IV_near_integrand_III, &args, T_m, flag);
            phi_mm = quadl->integral_4D_tri(phi_IV_near_integrand_I, &args, T_m, T_n, flag)+
                     quadl->integral_2D_tri(phi_IV_near_integrand_II, &args, T_m, flag);
            break;
        case 7:
            args.s_m = '+'; args.s_n = '+';
            psi_pp = quadl->integral_4D_tri(psi_IV_far_integrand, &args, T_m, T_n, flag);
            phi_pp = quadl->integral_4D_tri(phi_IV_far_integrand, &args, T_m, T_n, flag);
            args.s_m = '+'; args.s_n = '-';
            psi_pm = quadl->integral_4D_tri(psi_IV_far_integrand, &args, T_m, T_n, flag);
            phi_pm = quadl->integral_4D_tri(phi_IV_far_integrand, &args, T_m, T_n, flag);
            args.s_m = '-'; args.s_n = '+';
            psi_mp = quadl->integral_4D_tri(psi_IV_far_integrand, &args, T_m, T_n, flag);
            phi_mp = quadl->integral_4D_tri(phi_IV_far_integrand, &args, T_m, T_n, flag);
            args.s_m = '-'; args.s_n = '-';
            psi_mm = quadl->integral_4D_tri(psi_IV_far_integrand, &args, T_m, T_n, flag);
            phi_mm = quadl->integral_4D_tri(phi_IV_far_integrand, &args, T_m, T_n, flag);
            break;
    }
    A = psi_pp+psi_pm+psi_mp+psi_mm;
    B = phi_pp-phi_pm-phi_mp+phi_mm;
    return +j*k*eta*A-j*(eta/k)*B;
}

cmplx I_radiation_1_integrand(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    cmplx s1=args->s1;
    cmplx s2=args->s2;
    cmplx j=cmplx(0.0, 1.0);
    return alpha*exp(+j*s1*alpha)*exp(+j*s2*beta);
}

cmplx I_radiation_2_integrand(const cmplx alpha, const cmplx beta, void *args_in){
    Integrand_Args *args=(Integrand_Args*)args_in;
    cmplx s1=args->s1;
    cmplx s2=args->s2;
    cmplx j=cmplx(0.0, 1.0);
    return beta*exp(+j*s1*alpha)*exp(+j*s2*beta);
}

void I_radiation(cmplx s1, cmplx s2, cmplx &I1, cmplx &I2){
    const int N_quadl=4;
    const double tol=1.0E-4;
    const int k_max=10;
    Triangle_2D T;
    T.v1 = Vector(+0.0, +0.0, +0.0);
    T.v2 = Vector(+1.0, +0.0, +0.0);
    T.v3 = Vector(+0.0, +1.0, +0.0);
    Integrand_Args args;
    args.s1 = s1;
    args.s2 = s2;
    //
    QuadL quadl;
    quadl.set(N_quadl, tol, k_max);
    int flag;
    I1 = quadl.integral_2D_tri(I_radiation_1_integrand, &args, T, flag); assert(!flag);
    I2 = quadl.integral_2D_tri(I_radiation_2_integrand, &args, T, flag); assert(!flag);
}
