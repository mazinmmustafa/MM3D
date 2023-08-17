//
#include "matrix.hpp"

namespace basic_lib{
// Begin library

// Matrix class
Matrix::Matrix(){}

Matrix::~Matrix(){
    Matrix::deallocate();
}

void Matrix::allocate(const int rows, const int cols){
    if (this->is_allocated){
        Matrix::deallocate();
    }
    check_error(rows<1||cols<1, "incorrect matrix size!");
    this->matrix_size.rows = rows;
    this->matrix_size.cols = cols;
    this->data = (cmplx**)calloc(rows, sizeof(cmplx*));
    assert(this->data!=NULL);
    for (int i=0; i<rows; i++){
        this->data[i] = (cmplx*)calloc(cols, sizeof(cmplx));
        assert(this->data[i]!=NULL);
    }
    this->is_allocated = TRUE;
    Matrix::zeros();
}

void Matrix::deallocate(){
    if (this->is_allocated){
        for (int i=0; i<this->matrix_size.rows; i++){
            free(this->data[i]);
        }
        free(this->data);
    }
    if (this->is_lu_allocated){
        for (int i=0; i<this->matrix_size.rows; i++){
            free(this->lu_data[i]);
        }
        free(this->lu_data);
        free(this->P_data);
    }
    this->is_allocated = FALSE;
    this->is_lu_allocated = FALSE;
}

void Matrix::zeros(){
    check_error(!this->is_allocated, "matrix is not allocated yet! 1");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            this->data[i][j] = 0.0;
        }
    }
}

void Matrix::ones(){
    check_error(!this->is_allocated, "matrix is not allocated yet! 2");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            this->data[i][j] = 1.0;
        }
    }
}

void Matrix::eye(){
    check_error(!this->is_allocated, "matrix is not allocated yet! 3");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            i==j ? this->data[i][j] = 1.0 : this->data[i][j] = 0.0;
        }
    }
}

void Matrix::copy(Matrix &matrix){
    check_error(!this->is_allocated, "matrix is not allocated yet! 4");
    check_error(!matrix.is_allocated, "matrix is not allocated yet! 5");
    check_error((this->matrix_size.rows!=matrix.matrix_size.rows)||
        this->matrix_size.cols!=matrix.matrix_size.cols, "incompatible matrices!");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            matrix(i, j) = this->data[i][j];
        }
    }
}

void Matrix::copy_lu(Matrix &matrix){
    check_error(!this->is_lu_allocated, "lu matrix is not allocated yet! 6");
    check_error(!matrix.is_allocated, "matrix is not allocated yet! 7");
    check_error((this->matrix_size.rows!=matrix.matrix_size.rows)||
        this->matrix_size.cols!=matrix.matrix_size.cols, "incompatible matrices!");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            matrix(i, j) = this->lu_data[i][j];
        }
    }
}

void Matrix::save(const char *filename){
    check_error(!this->is_allocated, "matrix is not allocated yet! 8");
    File file;
    file.open(filename, "w");
    file.write("Matrix size: %d, %d\n\n", this->matrix_size.rows, this->matrix_size.cols);
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            file.write("Element %d, %d: ", i, j);
            file.write("(%21.14E, %21.14E)\n", real(this->data[i][j]), imag(this->data[i][j]));
        }
        file.write("\n");
    }
    file.close();
}

void Matrix::load(const char *filename){
    check_error(!this->is_allocated, "matrix is not allocated yet! 9");
    File file;
    file.open(filename, "r");
    int rows=0, cols=0;
    file.read("Matrix size: %d, %d\n\n", &rows, &cols);
    check_error(this->matrix_size.rows!=rows, "incompatible matrix size!");
    check_error(this->matrix_size.cols!=cols, "incompatible matrix size!");
    for (int i=0; i<this->matrix_size.rows; i++){
        for (int j=0; j<this->matrix_size.cols; j++){
            int ii, jj;
            file.read("Element %d, %d: ", &ii, &jj);
            assert(i==ii&&j==jj);
            double x, y;
            file.read("(%lf, %lf)\n", &x, &y);
            this->data[i][j] = cmplx(x, y);
        }
        file.read("\n");
    }
    file.close();
}

cmplx Matrix::operator () (const int i, const int j) const{
    check_error(i<0||j<0||i>=this->matrix_size.rows||j>=this->matrix_size.cols, 
        "matrix indices out of range!");
    return this->data[i][j];
}

cmplx& Matrix::operator () (const int i, const int j){
    check_error(i<0||j<0||i>=this->matrix_size.rows||j>=this->matrix_size.cols, 
        "matrix indices out of range!");
    return this->data[i][j];
}

Matrix_Size Matrix::size(){
    return this->matrix_size;
}

void Matrix::lup(){
    check_error(!this->is_allocated, "matrix is not yet located!");
    if (this->is_lu_allocated){
        return;
    }
    int M=this->matrix_size.rows;
    int N=this->matrix_size.cols;
    check_error(M!=N, "matrix is not square!");
    this->lu_data = (cmplx**)calloc(N, sizeof(cmplx*));
    assert(this->lu_data!=NULL);
    for (int i=0; i<N; i++){
        this->lu_data[i] = (cmplx*)calloc(N, sizeof(cmplx));
        assert(this->lu_data[i]!=NULL);
    }
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            this->lu_data[i][j] = this->data[i][j];
        }
    }
    this->is_lu_allocated = TRUE;
    cmplx *ptr;
    int i_max, i_temp; 
    double max_element, abs_element;
    this->P_data=(int*)calloc((N+1), sizeof(int));
    assert(this->P_data!=NULL);
    for (int i=0; i<=N; i++){
        P_data[i]=i;
    }
    for (int i=0; i<N; i++){
        max_element = 0.0;
        i_max = i;
        for (int k=i; k<N; k++){
            if ((abs_element= abs(this->lu_data[k][i]))>max_element){ 
                max_element = abs_element;
                i_max = k;
            }
        }
        if (max_element==0.0){
            printf("Warning: matrix is singular!\n");
        }
        if (i_max!=i) {
            i_temp = this->P_data[i];
            this->P_data[i] = this->P_data[i_max];
            this->P_data[i_max] = i_temp;
            ptr = this->lu_data[i];
            this->lu_data[i] = this->lu_data[i_max];
            this->lu_data[i_max] = ptr;
            P_data[N]++;
        }
        for (int j=i+1; j<N; j++){
            this->lu_data[j][i]/=this->lu_data[i][i];
            for (int k=i+1; k<N; k++){
                this->lu_data[j][k]-=this->lu_data[j][i]*this->lu_data[i][k];
            }    
        }
    }
    this->is_lu_allocated = TRUE;
}

void Matrix::solve(const Matrix &b, Matrix &x){
    check_error(!this->is_allocated, "matrix is not yet located! 10");
    check_error(!b.is_allocated, "matrix is not yet located! 11");
    check_error(!x.is_allocated, "matrix is not yet located! 12");
    check_error(this->matrix_size.rows!=b.matrix_size.rows, "incompatible matrices!");
    check_error(this->matrix_size.rows!=x.matrix_size.rows, "incompatible matrices!");
    check_error(b.matrix_size.cols!=1, "incompatible matrices!");
    check_error(x.matrix_size.cols!=1, "incompatible matrices!");
    if (!this->is_lu_allocated){
        Matrix::lup();
    }
    int N=this->matrix_size.rows;
    for (int i=0; i<N; i++){
        x(i, 0) = b(this->P_data[i], 0);
        for (int k=0; k<i; k++){
            x(i, 0)-=this->lu_data[i][k]*x(k, 0);
        }
    }
    for (int i=N-1; i>=0; i--){
        for (int k=i+1; k<N; k++){
            x(i, 0)-=this->lu_data[i][k]*x(k, 0);
        }
        x(i, 0)/=this->lu_data[i][i];
    }
}

void Matrix::inv(){
    check_error(!this->is_allocated, "matrix is not allocated yet! 13");
    if (!this->is_lu_allocated){
        Matrix::lup();
    }
	int N=this->matrix_size.rows;
    for (int j=0; j<N; j++) {
        for (int i=0; i<N; i++) {
            this->data[i][j] = this->P_data[i] == j ? 1.0 : 0.0;
            for (int k=0; k<i; k++){
                this->data[i][j]-=this->lu_data[i][k]*this->data[k][j];
            }
        }
        for (int i=N-1; i>=0; i--) {
            for (int k=i+1; k<N; k++){
                this->data[i][j]-=this->lu_data[i][k]*this->data[k][j];
            }
            this->data[i][j]/=this->lu_data[i][i];
        }
    }
}

cmplx Matrix::det(){
	if (!this->is_lu_allocated){
        Matrix::lup();
    }
	int N=this->matrix_size.rows;
    cmplx det=this->lu_data[0][0];
    for (int i=1; i<N; i++){
        det*=this->lu_data[i][i];
    }
    return (this->P_data[N]-N)%2 == 0 ? det : -det;
}

void Matrix::normalize(){
    check_error(!this->is_allocated, "matrix is not allocated yet! 14");
	int rows=this->matrix_size.rows;
    int cols=this->matrix_size.cols;
    double max=0.0;
    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            if (max<abs(this->data[i][j])){
                max = abs(this->data[i][j]);
            }
        }
    }
    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            this->data[i][j]/=max; 
        }
    }
}

// End of library
}