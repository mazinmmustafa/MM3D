import mesh
import sys
import os

try:
    lc = float(sys.argv[1])
    isLog = int(float(sys.argv[2]))
    assert(isLog==0 or isLog==1)
except:
    print("ERROR: Missing lc value!")
    exit(1)

# Dimensions
xc = 0.0
yc = 0.0
zc = 0.0
r = 3.2

# Write the mesh file
file = open("mesh/mesh.geo", "w")
file.write("SetFactory(\"OpenCASCADE\");\n")
file.write("Mesh.Algorithm = 5;\n")
file.write("Mesh.Smoothing = 3;\n")
file.write("// Mesh Size\n")
file.write("lc = {:21.14E};\n".format(lc))
file.write("xc = {:21.14E};\n".format(xc))
file.write("yc = {:21.14E};\n".format(yc))
file.write("zc = {:21.14E};\n".format(zc))
file.write("r = {:21.14E};\n".format(r))
file.write("\n")
file.write("// Sheet\n")
file.write("Sphere(1) = {xc, yc, zc, r};\n")
file.write("\n")
file.write("// Loops\n")
file.write("Physical Surface(1) = {1};\n")
file.write("MeshSize{ PointsOf{ Surface{:}; } } = lc;\n")
file.close()

if isLog==1:
    logFlag = True
else:
    logFlag = False

mesh.tryMesh(logFlag)
