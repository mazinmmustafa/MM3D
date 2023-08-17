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
L = 0.4
W = 0.6

# Write the mesh file
file = open("mesh/mesh.geo", "w")
file.write("Mesh.Algorithm = 5;\n")
file.write("Mesh.Smoothing = 3;\n")
file.write("// Mesh Size\n")
file.write("lc = {:21.14E};\n".format(lc))
file.write("\n")
file.write("// Points\n")
file.write("// Sheet\n")
file.write("Point(1) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, -W/2.0, 0.0, lc))
file.write("Point(2) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, -W/2.0, 0.0, lc))
file.write("Point(3) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, +W/2.0, 0.0, lc))
file.write("Point(4) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, +W/2.0, 0.0, lc))
file.write("\n")
file.write("// Lines\n")
file.write("Line(1) = {1, 2};\n")
file.write("Line(2) = {2, 3};\n")
file.write("Line(3) = {3, 4};\n")
file.write("Line(4) = {4, 1};\n")
file.write("\n")
file.write("// Loops\n")
file.write("Curve Loop(1) = {1, 2, 3, 4};\n")
file.write("Plane Surface(1) = {1};\n")
file.close()

if isLog==1:
    logFlag = True
else:
    logFlag = False

mesh.tryMesh("mesh.geo", logFlag)
