//
#include "shape.hpp"

// Shape class

Shape::Shape(){
}

void Shape::set(const double lambda){
    check_error(lambda<=0.0, "invalid value for wavelength!");
    this->lambda = lambda;
    this->is_set = TRUE;
}

Shape::~Shape(){
    free(this->bases_list);
    this->is_set = FALSE;
}

void Shape::mesh(const int is_log){
    check_error(!this->is_set, "shape is not set yet!");
    this->T.set();
    assert(!system("gmsh mesh/shape.geo -2 -format stl -save_all -o mesh/shape.stl"));
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import mesh\n");
    if (is_log){
        file.write("mesh.tryMesh(logOption=True)\n");
    }else{
        file.write("mesh.tryMesh(logOption=False)\n");
    }
    file.close();
    #ifdef _WIN64
        assert(!system("python mesh/shape.py"));
    #endif
    #ifdef __linux__
        assert(!system("python3 mesh/shape.py"));
    #endif
    this->T.unset();
}

double get_area(Vector v1, Vector v2, Vector v3){
    double a=(v1-v2).mag();
    double b=(v1-v3).mag();
    double c=(v2-v3).mag();
    double p=(a+b+c)/2.0;
    return sqrt(p*(p-a)*(p-b)*(p-c));
}

void get_basis(Basis &B){
    B.L = (B.e_p-B.e_m).mag();
        B.A_m = get_area(B.r_m, B.e_m, B.e_p);
        B.A_p = get_area(B.r_p, B.e_p, B.e_m);
        B.L_mm = +1.0*(B.e_m-B.r_m);
        B.L_pm = +1.0*(B.e_p-B.r_m);
        B.L_mp = -1.0*(B.e_m-B.r_p);
        B.L_pp = -1.0*(B.e_p-B.r_p);
        B.n_m = (B.L_mm^B.L_pm).unit();
        B.n_p = (B.L_pp^B.L_mp).unit();
}

void Shape::load_mesh(){
    check_error(!this->is_set, "shape is not set yet!");
    File file;
    file.open("mesh/mesh_info.txt", "r");
    file.read("%d", &this->N_bases);
    file.close();
    file.open("mesh/basis.dat", "r");
    this->bases_list = (Basis*)calloc(this->N_bases, sizeof(Basis));
    double x, y, z;
    printf("reading mesh data...");
    for (int i=0; i<this->N_bases; i++){
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].r_m.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].r_m.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].r_m.z = round_double(z, -int(log10(this->mesh_tol)));
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].e_m.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].e_m.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].e_m.z = round_double(z, -int(log10(this->mesh_tol)));
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].r_p.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].r_p.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].r_p.z = round_double(z, -int(log10(this->mesh_tol)));
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].e_p.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].e_p.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].e_p.z = round_double(z, -int(log10(this->mesh_tol)));
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].n_m.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].n_m.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].n_m.z = round_double(z, -int(log10(this->mesh_tol)));
        file.read("%lf", &x); file.read("%lf", &y); file.read("%lf", &z);
        this->bases_list[i].n_p.x = round_double(x, -int(log10(this->mesh_tol)));
        this->bases_list[i].n_p.y = round_double(y, -int(log10(this->mesh_tol)));
        this->bases_list[i].n_p.z = round_double(z, -int(log10(this->mesh_tol)));
        get_basis(this->bases_list[i]);
    }
    printf(", done!\n");
    file.close();
}

void Shape::log(const char *filename){
    check_error(!this->is_set, "shape is not set yet!");
    File file;
    file.open(filename, "w");
    file.write("number of bases functions %d\n\n", this->N_bases);
    for (int i=0; i<this->N_bases; i++){
        file.write("basis %d:\n", i);
        file.write("r_m: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_m.x),
            real(this->bases_list[i].r_m.y),
            real(this->bases_list[i].r_m.z));
        file.write("e_m: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].e_m.x),
            real(this->bases_list[i].e_m.y),
            real(this->bases_list[i].e_m.z));
        file.write("r_p: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_p.x),
            real(this->bases_list[i].r_p.y),
            real(this->bases_list[i].r_p.z));
        file.write("e_p: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].e_p.x),
            real(this->bases_list[i].e_p.y),
            real(this->bases_list[i].e_p.z));
        file.write("n_m: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].n_m.x),
            real(this->bases_list[i].n_m.y),
            real(this->bases_list[i].n_m.z));
        file.write("n_p: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].n_p.x),
            real(this->bases_list[i].n_p.y),
            real(this->bases_list[i].n_p.z));
        file.write("L_mm: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_mm.x),
            real(this->bases_list[i].L_mm.y),
            real(this->bases_list[i].L_mm.z));
        file.write("L_mp: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_mp.x),
            real(this->bases_list[i].L_mp.y),
            real(this->bases_list[i].L_mp.z));
        file.write("L_pm: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_pm.x),
            real(this->bases_list[i].L_pm.y),
            real(this->bases_list[i].L_pm.z));
        file.write("L_pp: (%21.14E, %21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_pp.x),
            real(this->bases_list[i].L_pp.y),
            real(this->bases_list[i].L_pp.z));
        file.write("L: %21.14E\n", this->bases_list[i].L);
        file.write("A_m: %21.14E\n", this->bases_list[i].A_m);
        file.write("A_p: %21.14E\n", this->bases_list[i].A_p);
        file.write("\n");
    }
    file.close();
}

void Shape::create_sphere(const double radius, const double lc){
    check_error(lc<=0.0, "invalid value for lc!");
    check_error(radius<=0.0, "invalid radius for a sphere!");
    //
    File file;
    file.open("mesh/shape.geo", "w");
    file.write("SetFactory(\"OpenCASCADE\");\n");
    file.write("Mesh.Algorithm = 6;\n");
    file.write("Mesh.Smoothing = 2;\n");
    file.write("Mesh.RefineSteps = 10;\n");
    file.write("Mesh.SmoothNormals = 1;\n");
    file.write("Mesh.SmoothRatio = 1.2;\n\n");
    file.write("// Mesh Size\n");
    file.write("lc = {%21.14E};\n", lc);
    file.write("xc = {%21.14E};\n", 0.0/this->lambda);
    file.write("yc = {%21.14E};\n", 0.0/this->lambda);
    file.write("zc = {%21.14E};\n", 0.0/this->lambda);
    file.write("r = {%21.14E};\n", radius/this->lambda);
    file.write("\n");
    file.write("// Shape\n");
    file.write("Sphere(1) = {xc, yc, zc, r};\n");
    file.write("\n");
    file.write("// Loops\n");
    file.write("Physical Surface(\"sphere\", 1) = {1};\n");
    file.write("MeshSize{ PointsOf{ Surface{:}; } } = lc;\n");
    file.close();
}

void Shape::create_sheet(const double l, const double w, const double lc){
    check_error(lc<=0.0, "invalid value for lc!");
    check_error(l<=0.0||w<=0.0, "invalid sheet dimensions!");
    //
    File file;
    file.open("mesh/shape.geo", "w");
    file.write("SetFactory(\"OpenCASCADE\");\n");
    file.write("Mesh.Algorithm = 6;\n");
    file.write("Mesh.Smoothing = 2;\n");
    file.write("Mesh.RefineSteps = 10;\n");
    file.write("Mesh.SmoothNormals = 1;\n");
    file.write("Mesh.SmoothRatio = 1.2;\n\n");
    file.write("// Mesh Size\n");
    file.write("lc = {%21.14E};\n", lc);
    file.write("Point(1) = {%21.14E, %21.14E,  %21.14E,  %21.14E};\n", -(l/2.0)/this->lambda, -(w/2.0)/this->lambda, 0.0, lc);
    file.write("Point(2) = {%21.14E, %21.14E,  %21.14E,  %21.14E};\n", +(l/2.0)/this->lambda, -(w/2.0)/this->lambda, 0.0, lc);
    file.write("Point(3) = {%21.14E, %21.14E,  %21.14E,  %21.14E};\n", +(l/2.0)/this->lambda, +(w/2.0)/this->lambda, 0.0, lc);
    file.write("Point(4) = {%21.14E, %21.14E,  %21.14E,  %21.14E};\n", -(l/2.0)/this->lambda, +(w/2.0)/this->lambda, 0.0, lc);
    file.write("// Sheet\n");
    file.write("// Lines\n");
    file.write("Line(1) = {1, 2};\n");
    file.write("Line(2) = {2, 3};\n");
    file.write("Line(3) = {3, 4};\n");
    file.write("Line(4) = {4, 1};\n");
    file.write("// Loops\n");
    file.write("Curve Loop(1) = {1, 2, 3, 4};\n");
    file.write("Plane Surface(1) = {1};\n");
    file.write("Physical Surface(\"sheet\", 1) = {1};\n");
    file.write("MeshSize{ PointsOf{ Surface{:}; } } = lc;\n");
    file.close();
}
