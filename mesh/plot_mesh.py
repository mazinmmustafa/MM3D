from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

fileName = "mesh/mesh_data.dat"
fileBasis = "mesh/basis.dat"

MeshData = []
with open(fileName, "r") as file:
    for newLine in file:
        MeshData.append(newLine.split())
        continue

BasisData = []
with open(fileBasis, "r") as file:
    for newLine in file:
        BasisData.append(newLine.split())
        continue

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

for tri in MeshData:
    x = [float(tri[0]), float(tri[3]), float(tri[6]), float(tri[0])]
    y = [float(tri[1]), float(tri[4]), float(tri[7]), float(tri[1])]
    z = [float(tri[2]), float(tri[5]), float(tri[8]), float(tri[2])]
    vertices = [list(zip(x,y,z))]
    poly = Poly3DCollection(vertices, alpha=0.0)
    ax.add_collection3d(poly)

    xL, yL, zL = [x[0], x[1]], [y[0], y[1]], [z[0], z[1]]
    # ax.scatter(xL, yL, zL, c='red', s=2)
    ax.plot(xL, yL, zL, color='black')
    xL, yL, zL = [x[1], x[2]], [y[1], y[2]], [z[1], z[2]]
    # ax.scatter(xL, yL, zL, c='red', s=2)
    ax.plot(xL, yL, zL, color='black')
    xL, yL, zL = [x[2], x[0]], [y[2], y[0]], [z[2], z[0]]
    # ax.scatter(xL, yL, zL, c='red', s=2)
    ax.plot(xL, yL, zL, color='black')

for basis in BasisData:
    x = [float(basis[0]), float(basis[3]), float(basis[6]), float(basis[9])]
    y = [float(basis[1]), float(basis[4]), float(basis[7]), float(basis[10])]
    z = [float(basis[2]), float(basis[5]), float(basis[8]), float(basis[11])]

    # Red->Blue
    # xL, yL, zL = [x[0], 0.5*(x[1]+x[3])], [y[0], 0.5*(y[1]+y[3])], [z[0], 0.5*(z[1]+z[3])]
    # ax.plot(xL, yL, zL, color='red', linewidth=1)
    # xL, yL, zL = [x[2], 0.5*(x[1]+x[3])], [y[2], 0.5*(y[1]+y[3])], [z[2], 0.5*(z[1]+z[3])]
    # ax.plot(xL, yL, zL, color='blue', linewidth=1)

    # Mesh
    xL, yL, zL = [x[0], x[1], x[2], x[3]], [y[0], y[1], y[2], y[3]], [z[0], z[1], z[2], z[3]]
    ax.plot(xL, yL, zL, color='black', linewidth=0.2)
    
ax.set_aspect('equal', adjustable='box', anchor='C')
ax.set_xlim(-1, +1)
ax.set_ylim(-1, +1)
ax.set_zlim(-1, +1)

plt.show()