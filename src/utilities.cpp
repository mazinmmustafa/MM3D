//
#include "utilities.hpp"

namespace basic_lib{
// Begin library

// Check error function
void check_error(const int condition, const char *error_msg){
    if (condition){
        printf("Error: %s\n", error_msg);
        exit(1);
    }
}

// Progress bar function
void progress_bar(const unsigned long i, const unsigned long N, const char *msg){
    check_error(!(i<N), "incorrect counter for progress bar!");
    const int N_bar=20;
    unsigned long n=i;
    n++;
    printf("\rProgress: [");
    for (unsigned long i=0; i<=n*N_bar/N; i++){
        printf("#");
    }
    for (unsigned long i=n*N_bar/N; i<N_bar; i++){
        printf("-");
    }
    printf("] %3ld%%", (int)100*n/N);
    printf(" %s", msg);
    if (n==N){
        printf(", done!\n");
    }else
    if (n>N){
        printf("\n");
    }
    fflush(stdout);
}

/*
// usleep function
#ifdef _WIN64
void usleep(__int64 usec){ 
    HANDLE timer; 
    LARGE_INTEGER ft; 
    ft.QuadPart = -(10*usec); 
    timer = CreateWaitableTimer(NULL, TRUE, NULL); 
    SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0); 
    WaitForSingleObject(timer, INFINITE); 
    CloseHandle(timer); 
}
#endif
*/

// Misc functions
void disp(const cmplx z){
    printf("(%21.14E, %21.14E)\n", real(z), imag(z));
}

void disp(const double x){
    printf("%21.14E\n", x);
}

void disp(const int n){
    printf("%d\n", n);
}

double deg2rad(const double x){
    return x*pi/180.0;
}

double rad2deg(const double x){
    return x*180.0/pi;
}

double sinc(double x){
    return abs(x)==0.0 ? 1.0 : sin(x)/x;
}

cmplx sinc(cmplx z){
    return abs(z)==0.0 ? 1.0 : sin(z)/z;
}

double round_double(const double x, const int n){
    double r=pow(10.0, n);
    double y=round(x*r);
    return y/r;
}

// File class

File::File(){}

File::~File(){}

void File::open(const char *filename, const char *option){
    check_error(this->is_open, "file is already open!");
    check_error(filename==NULL, "invalid filename!");
    int filename_len=(int)strlen(filename);
    check_error(!(!strcmp(option, "w")||!strcmp(option, "a")||!strcmp(option, "r")), 
        "invalid file option !");
    this->filename = (char*)calloc(filename_len+1, sizeof(char));
    assert(this->filename!=NULL);
    this->option = (char*)calloc(2, sizeof(char));
    assert(this->option!=NULL);
    assert(strcpy(this->filename, filename));
    assert(strcpy(this->option, option));
    this->file_ptr=fopen(this->filename, this->option);
    check_error(this->file_ptr==NULL, "unable to open file!");
    this->is_open = TRUE;
}

void File::close(){
    check_error(!this->is_open, "file is already closed!");
    free(this->filename);
    free(this->option);
    fclose(this->file_ptr);
    this->is_open = FALSE;
}

void File::write(const char *format, ...){
    check_error(!this->is_open, "can't write to unopened file!");
    va_list args;
    va_start(args, format);
    vfprintf(this->file_ptr, format, args);
    va_end(args);
}

int File::read(const char *format, ...){
    check_error(!this->is_open, "can't read from unopened file!");
    va_list args;
    va_start(args, format);
    int n=vfscanf(this->file_ptr, format, args);
    assert(n>-1);
    va_end(args);
    return feof(this->file_ptr);
}

// Timer class

Timer::Timer(){}

Timer::~Timer(){}

void Timer::set(){
    #ifdef _WIN64
    this->start = clock();
    #endif
    #ifdef __linux__
    timespec_get(&(this->start), TIME_UTC);
    #endif
    this->is_set = TRUE;
}

void Timer::unset(){
    if (!this->is_set){
        Timer::set();
    }
    #ifdef _WIN64
        this->stop = clock();
        this->elapsed = (double)(this->stop-this->start)/CLOCKS_PER_SEC;
    #endif
    #ifdef __linux__
        timespec_get(&(this->stop), TIME_UTC);
        this->elapsed = (double)(this->stop.tv_sec-this->start.tv_sec)+
            ((double)(this->stop.tv_nsec-this->start.tv_nsec)/1000000000L);
    #endif    
    std::cout << "Elapsed time is ";
    if (this->elapsed<60.0){
        std::cout << this->elapsed << " seconds" << std::endl;
    }else
    if (this->elapsed<3600.0){
        std::cout << this->elapsed/60.0 << " minutes" << std::endl;
    }else{
        std::cout << this->elapsed/3600.0 << " hours" << std::endl;
    }
    this->is_set = FALSE;
}

void Timer::unset_silent(){
    #ifdef _WIN64
        this->stop = clock();
        this->elapsed = (double)(this->stop-this->start)/CLOCKS_PER_SEC;
    #endif
    #ifdef __linux__
        timespec_get(&(this->stop), TIME_UTC);
        this->elapsed = (double)(this->stop.tv_sec-this->start.tv_sec)+
            ((double)(this->stop.tv_nsec-this->start.tv_nsec)/1000000000L);
    #endif   
}

double Timer::get_elapsed(){
    return this->elapsed;
}

// Random class
void Random::set_seed(){
    srand(time(&this->t));
}

Random::Random(){
    Random::set_seed();
}

Random::~Random(){}

int Random::rand_int(const int a, const int b){
    return b>a ? a+rand()%(b-a+1) : 0;
}

double Random::rand_double(const double a, const double b){
    return b>a ? a+rand()/(RAND_MAX/(b-a)) : 0.0;
}

// Range class

Range::Range(){}

Range::~Range(){
    Range::deallocate();
}

void Range::linspace(const double x_min, const double x_max, const int Ns){
    Range::deallocate();
    check_error(x_max<=x_min, "invalid range!");
    check_error(Ns<1, "invalid number of points!");
    this->Ns = Ns;
    this->data = (double*)calloc(Ns, sizeof(double));
    assert(this->data!=NULL);
    double dx=(x_max-x_min)/(Ns-1.0);
    for (int i=0; i<Ns; i++){
        this->data[i] = x_min+i*dx;
    }
    this->is_allocated = TRUE;
}

void Range::logspace(const double x_min, const double x_max, const int Ns){
    Range::deallocate();
    check_error(x_max<=x_min, "invalid range!");
    check_error(Ns<1, "invalid number of points!");
    this->Ns = Ns;
    this->data = (double*)calloc(Ns, sizeof(double));
    assert(this->data!=NULL);
    double X_min=log10(x_min);
    double X_max=log10(x_max);
    double dX=(X_max-X_min)/(Ns-1.0);
    for (int i=0; i<Ns; i++){
        this->data[i] = pow(10.0, X_min+i*dX);
    }
    this->is_allocated = TRUE;
}

double Range::operator() (const int i) const{
    check_error(!this->is_allocated, "range is not allocated yet!");
    check_error(!(i>=0&&i<this->Ns), "index is out of range!");
    return this->data[i];
}

void Range::deallocate(){
    if (this->is_allocated){
        free(this->data);
        this->data = NULL;
        this->is_allocated = FALSE;
    }
}

int Range::size(){
    return this->Ns;
}

// End of library
}