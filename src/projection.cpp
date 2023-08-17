//
#include "projection.hpp"

//  

void projection_1D(Vector v1, Vector v2, Vector p, Projection_1D_Para &para){
    check_error(v1.is_equal(v2), "projection is undefined!");
    Vector v21=v2-v1;
    double alpha=real(v21*(p-v1)/pow(v21.mag(), 2.0));
    Vector p0=v1+alpha*v21;
    para.d = (p0-p).mag();
    para.P_m = (p-v1).mag();
    para.P_p = (p-v2).mag();
    para.p0 = p0;
    para.P0 = p0-p;
    if (alpha<0.0){
        para.l_m = +(p0-v1).mag();
        para.l_p = +(p0-v2).mag();
    }else
    if (alpha>1.0){
        para.l_m = -(p0-v1).mag();
        para.l_p = -(p0-v2).mag();
    }else{
        para.l_m = -(p0-v1).mag();
        para.l_p = +(p0-v2).mag();
    }
}

void projection_3D(const Triangle T, const Vector p, Projection_3D_Para &para){
    Vector v1=T.v1;
    Vector v2=T.v2;
    Vector v3=T.v3;
    Vector v21=v2-v1;
    Vector v31=v3-v1;
    //
    cmplx a1=v21*v21; cmplx a2=v21*v31; cmplx a3=v21*(v1-p);
    cmplx b1=v31*v31; cmplx b2=v31*v21; cmplx b3=v31*(v1-p);
    cmplx alpha=(a2*b3-b1*a3)/(a1*b1-a2*b2);
    cmplx beta =(b2*a3-a1*b3)/(a1*b1-a2*b2);
    Vector p0=v1+alpha*v21+beta*v31;
    //
    para.p0 = p0;
    para.d = (p0-p).mag();
    para.R_m = (p-v2).mag();
    para.R_p = (p-v3).mag();
    Projection_1D_Para para_1D;
    projection_1D(v2, v3, p0, para_1D);
    para.l_m = para_1D.l_m;
    para.l_p = para_1D.l_p;
    para.P0 = para_1D.P0;
    para.R0 = (para_1D.p0-p).mag();
    projection_1D(v2, v3, v1, para_1D);
    para.u = para_1D.P0.unit();
}

cmplx T1_projection_terms(Projection_3D_Para para){
    cmplx factor = (para.P0.unit()*para.u);
    double P0=para.P0.mag();
    double R0=para.R0;
    double R_p=para.R_p;
    double R_m=para.R_m;
    double l_p=para.l_p;
    double l_m=para.l_m;
    double d=para.d;
    double A=P0*log((R_p+l_p)/(R_m+l_m));
    double B=atan2(P0*l_p, R0*R0+d*R_p);
    double C=atan2(P0*l_m, R0*R0+d*R_m);
    return factor*(A-d*(B-C));
}

cmplx T1_projection(const cmplx alpha, const cmplx beta,
    Basis *B_m, Basis *B_n, const char s_m, const char s_n){
    Vector p;
    Projection_3D_Para para;
    Triangle T;
    cmplx sum=0.0;
    int flag=0;    
    cmplx factor;
    if (s_m=='+'){
        p = B_m->r_p-alpha*B_m->L_pp-beta*B_m->L_mp;
        flag++;
    }
    if (s_m=='-'){
        p = B_m->r_m+alpha*B_m->L_mm+beta*B_m->L_pm;
        flag++;
    }
    if (s_n=='+'){
        // Side 1
        T.v1 = B_n->r_p; T.v2 = B_n->e_p; T.v3 = B_n->e_m;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        // Side 2
        T.v1 = B_n->e_p; T.v2 = B_n->e_m; T.v3 = B_n->r_p;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        // Side 3
        T.v1 = B_n->e_m; T.v2 = B_n->r_p; T.v3 = B_n->e_p;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        flag++;
    }
    if (s_n=='-'){
        // Side 1
        T.v1 = B_n->r_m; T.v2 = B_n->e_m; T.v3 = B_n->e_p;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        // Side 2
        T.v1 = B_n->e_m; T.v2 = B_n->e_p; T.v3 = B_n->r_m;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        // Side 3
        T.v1 = B_n->e_p; T.v2 = B_n->r_m; T.v3 = B_n->e_m;
        projection_3D(T, p, para);
        sum+=T1_projection_terms(para);
        flag++;
    }
    assert(flag==2);
    return sum;
}

Vector T2_projection_terms(Projection_3D_Para para){
    double R0=para.R0;
    double R_p=para.R_p;
    double R_m=para.R_m;
    double l_p=para.l_p;
    double l_m=para.l_m;
    double terms=R0*R0*log((R_p+l_p)/(R_m+l_m))+R_p*l_p-R_m*l_m;
    return 0.5*para.u*terms;
}

Vector T2_projection(const cmplx alpha, const cmplx beta,
    Basis *B_m, Basis *B_n, const char s_m, const char s_n){
    Vector p;
    Projection_3D_Para para;
    Triangle T;
    Vector sum;
    int flag=0;    
    cmplx factor;
    if (s_m=='+'){
        p = B_m->r_p-alpha*B_m->L_pp-beta*B_m->L_mp;
        flag++;
    }
    if (s_m=='-'){
        p = B_m->r_m+alpha*B_m->L_mm+beta*B_m->L_pm;
        flag++;
    }
    if (s_n=='+'){
        // Side 1
        T.v1 = B_n->r_p; T.v2 = B_n->e_p; T.v3 = B_n->e_m;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        // Side 2
        T.v1 = B_n->e_p; T.v2 = B_n->e_m; T.v3 = B_n->r_p;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        // Side 3
        T.v1 = B_n->e_m; T.v2 = B_n->r_p; T.v3 = B_n->e_p;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        flag++;
    }
    if (s_n=='-'){
        // Side 1
        T.v1 = B_n->r_m; T.v2 = B_n->e_m; T.v3 = B_n->e_p;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        // Side 2
        T.v1 = B_n->e_m; T.v2 = B_n->e_p; T.v3 = B_n->r_m;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        // Side 3
        T.v1 = B_n->e_p; T.v2 = B_n->r_m; T.v3 = B_n->e_m;
        projection_3D(T, p, para);
        sum = sum + T2_projection_terms(para);
        flag++;
    }
    assert(flag==2);
    return sum;
}
