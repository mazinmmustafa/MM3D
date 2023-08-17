//
#include "quadl.hpp"

namespace basic_lib{
// Begin library
extern "C" void cgqf_f77_(int *rule, int *order, double *x, double *w);

QuadL::QuadL(){}

QuadL::~QuadL(){
    QuadL::unset();
}

void QuadL::set(const int N_quadl, const double tol, const int k_max){
    check_error(N_quadl<1, "invalid quadrature order!");
    check_error(tol<=0.0, "invalid quadrature tolerance!");
    check_error(k_max<0, "invalid recursion depth!");
    if (this->is_set){
        QuadL::unset();
    }
    this->N_quadl = N_quadl;
    this->tol = tol;
    this->k_max = k_max;
    this->x = (double*)calloc(N_quadl, sizeof(double));
    assert(this->x!=NULL);
    this->w = (double*)calloc(N_quadl, sizeof(double));
    assert(this->w!=NULL);
    int rule=1;
    cgqf_f77_(&rule, &(this->N_quadl), this->x, this->w);
    this->is_set = TRUE;
}

void QuadL::unset(){
    if (this->is_set){
        free(this->x);
        free(this->w);
    }
    this->is_set = FALSE;
}

cmplx QuadL::integral_1D_quadl(cmplx func(cmplx, void*), void *args, 
    const double a, const double b){
    double h_p=(b+a)/2.0;
    double h_m=(b-a)/2.0;
    cmplx sum=0.0;
    double w_n;
    double x_n;
    for (int n=0; n<this->N_quadl; n++){
        w_n = this->w[n];
        x_n = this->x[n];
        sum+=h_m*w_n*func(h_m*x_n+h_p, args);
    }
    return sum;
}

cmplx QuadL::integral_1D_recursive(cmplx func(cmplx, void*), void *args, 
    const double a, const double b, int &flag, int k, cmplx I_p){
    double m=(a+b)/2.0;
    cmplx I1=QuadL::integral_1D_quadl(func, args, a, m);
    cmplx I2=QuadL::integral_1D_quadl(func, args, m, b);
    cmplx I_n=I1+I2;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<this->k_max&&error>0.0){
        k++;
        I1 = QuadL::integral_1D_recursive(func, args, a, m, flag, k, I1);
        I2 = QuadL::integral_1D_recursive(func, args, m, b, flag, k, I2);
        I_n = I1+I2;
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_1D(cmplx func(cmplx, void*), void *args, 
    const double a, const double b, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_1D_recursive(func, args, a, b, flag, 0, 0.0);
}

cmplx QuadL::integral_2D_quadl(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y){
    double h_x_p=(b_x+a_x)/2.0;
    double h_x_m=(b_x-a_x)/2.0;
    double h_y_p=(b_y+a_y)/2.0;
    double h_y_m=(b_y-a_y)/2.0;
    cmplx sum=0.0;
    double w_m;
    double x_m;
    double w_n;
    double x_n;
    for (int m=0; m<this->N_quadl; m++){
        w_m = this->w[m];
        x_m = this->x[m];
        for (int n=0; n<this->N_quadl; n++){
            w_n = this->w[n];
            x_n = this->x[n];
            sum+=h_x_m*h_y_m*w_m*w_n*func(h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, args);
        }
    }
    return sum;
}

cmplx QuadL::integral_2D_recursive(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y, int &flag, int k, cmplx I_p){
    double m_x=(a_x+b_x)/2.0;
    double m_y=(a_y+b_y)/2.0;
    cmplx I1=integral_2D_quadl(func, args, a_x, m_x, a_y, m_y);
    cmplx I2=integral_2D_quadl(func, args, a_x, m_x, m_y, b_y);
    cmplx I3=integral_2D_quadl(func, args, m_x, b_x, a_y, m_y);
    cmplx I4=integral_2D_quadl(func, args, m_x, b_x, m_y, b_y);
    cmplx I_n=I1+I2+I3+I4;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<k_max&&error>0.0){
        k++;
        I1 = integral_2D_recursive(func, args, a_x, m_x, a_y, m_y, flag, k, I1);
        I2 = integral_2D_recursive(func, args, a_x, m_x, m_y, b_y, flag, k, I2);
        I3 = integral_2D_recursive(func, args, m_x, b_x, a_y, m_y, flag, k, I3);
        I4 = integral_2D_recursive(func, args, m_x, b_x, m_y, b_y, flag, k, I4);
        I_n = I1+I2+I3+I4; 
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_2D(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_2D_recursive(func, args, a_x, b_x, a_y, b_y, flag, 0, 0.0);
}

static double get_area(const Triangle_2D T){
    double a=(T.v1-T.v2).mag();
    double b=(T.v1-T.v3).mag();
    double c=(T.v2-T.v3).mag();
    double p=(a+b+c)/2.0;
    return sqrt(p*(p-a)*(p-b)*(p-c));
}

cmplx QuadL::integral_2D_tri_quadl(cmplx func(cmplx, cmplx, void*), void *args, 
    const Triangle_2D T){
    double A=get_area(T);
    cmplx sum=0.0;
    double alpha_mn;
    double beta_mn;
    double w_mn;    
    double x_m, x_n;
    double w_m, w_n;
    cmplx alpha, beta;
    for (int m=0; m<this->N_quadl; m++){
        x_m = this->x[m];
        w_m = this->w[m];
        for (int n=0; n<this->N_quadl; n++){
            x_n = this->x[n];
            w_n = this->w[n];
            alpha_mn = (1.0/8.0)*(3.0+3.0*x_m-x_n-x_m*x_n);
            beta_mn  = (1.0/8.0)*(3.0+3.0*x_n-x_m-x_m*x_n);
            w_mn = (1.0/16.0)*(2.0-x_m-x_n)*w_m*w_n;
            alpha = (T.v2.x-T.v1.x)*alpha_mn+(T.v3.x-T.v1.x)*beta_mn+T.v1.x;
            beta  = (T.v2.y-T.v1.y)*alpha_mn+(T.v3.y-T.v1.y)*beta_mn+T.v1.y;
            sum+=w_mn*func(alpha, beta, args);
        }
    }
    return 2.0*A*sum;
}

cmplx QuadL::integral_2D_tri_recursive(cmplx func(cmplx, cmplx, void*), void *args, 
    const Triangle_2D T, int &flag, int k, cmplx I_p){
    Triangle_2D T1, T2;
    Vector m=0.5*(T.v2+T.v3);
    T1.v1 = m; T1.v2 = T.v1; T1.v3 = T.v2;
    T2.v1 = m; T2.v2 = T.v3; T2.v3 = T.v1;
    cmplx I1=integral_2D_tri_quadl(func, args, T1);
    cmplx I2=integral_2D_tri_quadl(func, args, T2);
    cmplx I_n=I1+I2;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<k_max&&error>0.0){
        k++;
        I1 = integral_2D_tri_recursive(func, args, T1, flag, k, I1);
        I2 = integral_2D_tri_recursive(func, args, T2, flag, k, I2);
        I_n = I1+I2; 
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_2D_tri(cmplx func(cmplx, cmplx, void*), void *args, 
    const Triangle_2D T, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_2D_tri_recursive(func, args, T, flag, 0, 0.0);
}

cmplx QuadL::integral_4D_tri_quadl(cmplx func(cmplx, cmplx, cmplx, cmplx, void*), void *args, 
    const Triangle_2D Ta, const Triangle_2D Tb){
    double A=get_area(Ta);
    double A_=get_area(Tb);
    cmplx sum=0.0;
    double alpha_mn;
    double beta_mn;
    double alpha_mn_;
    double beta_mn_;
    double w_mn, w_mn_;    
    double x_m, x_n, x_m_, x_n_;
    double w_m, w_n, w_m_, w_n_;
    cmplx alpha, beta;
    cmplx alpha_, beta_;
    for (int m=0; m<this->N_quadl; m++){
        x_m = this->x[m];
        w_m = this->w[m];
        for (int n=0; n<this->N_quadl; n++){
            x_n = this->x[n];
            w_n = this->w[n];
            for (int m_=0; m_<this->N_quadl; m_++){
                x_m_ = this->x[m_];
                w_m_ = this->w[m_];
                for (int n_=0; n_<this->N_quadl; n_++){
                    x_n_ = this->x[n_];
                    w_n_ = this->w[n_];
                    alpha_mn = (1.0/8.0)*(3.0+3.0*x_m-x_n-x_m*x_n);
                    beta_mn  = (1.0/8.0)*(3.0+3.0*x_n-x_m-x_m*x_n);
                    w_mn = (1.0/16.0)*(2.0-x_m-x_n)*w_m*w_n;
                    alpha = (Ta.v2.x-Ta.v1.x)*alpha_mn+(Ta.v3.x-Ta.v1.x)*beta_mn+Ta.v1.x;
                    beta  = (Ta.v2.y-Ta.v1.y)*alpha_mn+(Ta.v3.y-Ta.v1.y)*beta_mn+Ta.v1.y;
                    alpha_mn_ = (1.0/8.0)*(3.0+3.0*x_m_-x_n_-x_m_*x_n_);
                    beta_mn_  = (1.0/8.0)*(3.0+3.0*x_n_-x_m_-x_m_*x_n_);
                    w_mn_ = (1.0/16.0)*(2.0-x_m_-x_n_)*w_m_*w_n_;
                    alpha_ = (Tb.v2.x-Tb.v1.x)*alpha_mn_+(Tb.v3.x-Tb.v1.x)*beta_mn_+Tb.v1.x;
                    beta_  = (Tb.v2.y-Tb.v1.y)*alpha_mn_+(Tb.v3.y-Tb.v1.y)*beta_mn_+Tb.v1.y;
                    sum+=w_mn*w_mn_*func(alpha, beta, alpha_, beta_, args);
                }
            }
        }
    }
    return 4.0*A*A_*sum;
}

cmplx QuadL::integral_4D_tri_recursive(cmplx func(cmplx, cmplx, cmplx, cmplx, void*), void *args, 
    const Triangle_2D Ta, const Triangle_2D Tb, int &flag, int k, cmplx I_p){
    Triangle_2D T1, T2, T1_, T2_;
    Vector m=0.5*(Ta.v2+Ta.v3);
    Vector m_=0.5*(Tb.v2+Tb.v3);
    T1.v1 = m; T1.v2 = Ta.v1; T1.v3 = Ta.v2;
    T2.v1 = m; T2.v2 = Ta.v3; T2.v3 = Ta.v1;
    T1_.v1 = m_; T1_.v2 = Tb.v1; T1_.v3 = Tb.v2;
    T2_.v1 = m_; T2_.v2 = Tb.v3; T2_.v3 = Tb.v1;
    cmplx I1=integral_4D_tri_quadl(func, args, T1, T1_);
    cmplx I2=integral_4D_tri_quadl(func, args, T1, T2_);
    cmplx I3=integral_4D_tri_quadl(func, args, T2, T1_);
    cmplx I4=integral_4D_tri_quadl(func, args, T2, T2_);
    cmplx I_n=I1+I2+I3+I4;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<k_max&&error>0.0){
        k++;
        I1 = integral_4D_tri_recursive(func, args, T1, T1_, flag, k, I1);
        I2 = integral_4D_tri_recursive(func, args, T1, T2_, flag, k, I2);
        I3 = integral_4D_tri_recursive(func, args, T2, T1_, flag, k, I3);
        I4 = integral_4D_tri_recursive(func, args, T2, T2_, flag, k, I4);
        I_n = I1+I2+I3+I4; 
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_4D_tri(cmplx func(cmplx, cmplx, cmplx, cmplx, void*), void *args, 
    const Triangle_2D Ta, const Triangle_2D Tb, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_4D_tri_recursive(func, args, Ta, Tb, flag, 0, 0.0);
}

// End of library
}