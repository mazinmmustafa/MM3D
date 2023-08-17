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
L = 1.0
W = 1.0
H = 1.0

# Write the mesh file
file = open("mesh/mesh.geo", "w")
file.write("Mesh.Algorithm = 5;\n")
file.write("Mesh.Smoothing = 3;\n")
file.write("// Mesh Size\n")
file.write("lc = {:21.14E};\n".format(lc))
file.write("\n")
file.write("// Points\n")
file.write("// Sheet\n")
file.write("Point(1) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, -W/2.0, -H/2.0, lc))
file.write("Point(2) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, -W/2.0, -H/2.0, lc))
file.write("Point(3) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, +W/2.0, -H/2.0, lc))
file.write("Point(4) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, +W/2.0, -H/2.0, lc))
file.write("Point(5) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, +W/2.0, +H/2.0, lc))
file.write("Point(6) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(+L/2.0, -W/2.0, +H/2.0, lc))
file.write("Point(7) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, +W/2.0, +H/2.0, lc))
file.write("Point(8) = {{{:21.14E}, {:21.14E}, {:21.14E}, {:21.14E}}};\n".format(-L/2.0, -W/2.0, +H/2.0, lc))
file.write("\n")
file.write("// Lines\n")
file.write("Line(1) = {3, 7};\n")
file.write("Line(2) = {7, 5};\n")
file.write("Line(3) = {5, 4};\n")
file.write("Line(4) = {4, 3};\n")
file.write("Line(5) = {3, 1};\n")
file.write("Line(6) = {2, 4};\n")
file.write("Line(7) = {2, 6};\n")
file.write("Line(8) = {6, 8};\n")
file.write("Line(9) = {8, 1};\n")
file.write("Line(10) = {1, 2};\n")
file.write("Line(11) = {8, 7};\n")
file.write("Line(12) = {6, 5};\n")
file.write("\n")
file.write("// Loops\n")
file.write("Line Loop(13) = {7, 8, 9, 10};\n")
file.write("Plane Surface(14) = {13};\n")
file.write("Line Loop(15) = {-6, -10, -5, -4};\n")
file.write("Plane Surface(16) = {15};\n")
file.write("Line Loop(17) = {3, 4, 1, 2};\n")
file.write("Plane Surface(18) = {17};\n")
file.write("Line Loop(19) = {12, -2, -11, -8};\n")
file.write("Plane Surface(20) = {19};\n")
file.write("Line Loop(21) = {-7, 6, -3, -12};\n")
file.write("Plane Surface(22) = {21};\n")
file.write("Line Loop(23) = {-9, 11, -1, 5};\n")
file.write("Plane Surface(24) = {23};\n")
file.write("Surface Loop(25) = {14, 22, 20, 18, 16, 24};\n")
file.close()

if isLog==1:
    logFlag = True
else:
    logFlag = False

mesh.tryMesh("mesh.geo", logFlag)
