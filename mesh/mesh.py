def ProcessMesh(meshFileName, trianglesFileName):
    # Call the raw mesh file
    print("Processing \"{}\"...".format(meshFileName))

    # Class Vertex
    class Vertex():
        def __init__(self, x=0, y=0, z=0):
            self.x = x
            self.y = y
            self.z = z
            pass

    # Class Triangle
    class Triangle():
        def __init__(self, v=[], n=0):
            self.v = v
            self.n = n
            pass

    # Obtain triangles from mesh file
    triangles = []
    n = Vertex()
    T = Triangle()
    N = 0
    startFlag = False
    with open(meshFileName, "r") as file:
        for line in file:
            parsedLine = line.split()
            if parsedLine[0]=="facet" and startFlag==False:
                n.x = float(parsedLine[2])
                n.y = float(parsedLine[3])
                n.z = float(parsedLine[4])
                startFlag = True
                count = 0
                T.n = n
                T.v = []
                pass
            if parsedLine[0]=="vertex" and startFlag==True and count<3:
                v = Vertex(float(parsedLine[1]),
                           float(parsedLine[2]),
                           float(parsedLine[3]))
                T.v.append(v)
                count+=1
                if count==3:
                    newVertecies = [T.v[0], T.v[1], T.v[2]]
                    newNormal = Vertex(T.n.x, T.n.y, T.n.z)
                    newTriangle = Triangle(newVertecies, newNormal)
                    triangles.append(newTriangle)
                    startFlag = False
                    N+=1
                    pass
                pass
            continue
        pass

    # Output number of triangles
    print("{} triangles have been extracted".format(N))

    # Write triagnles file
    trianglesFile = open(trianglesFileName, "w")
    for i in range(0, N):
        trianglesFile.write("{:22.14E} {:22.14E} {:22.14E}".format(triangles[i].v[0].x, triangles[i].v[0].y, triangles[i].v[0].z))
        trianglesFile.write("{:22.14E} {:22.14E} {:22.14E}".format(triangles[i].v[1].x, triangles[i].v[1].y, triangles[i].v[1].z))
        trianglesFile.write("{:22.14E} {:22.14E} {:22.14E}".format(triangles[i].v[2].x, triangles[i].v[2].y, triangles[i].v[2].z))
        trianglesFile.write("{:22.14E} {:22.14E} {:22.14E}".format(triangles[i].n.x, triangles[i].n.y, triangles[i].n.z))
        trianglesFile.write("\n")
        pass

    trianglesFile.close()
    pass
    
def creatBasis(meshFileName, basisFileName, logOption):
    # Call the raw mesh file
    logFileName = "mesh/log.txt"
    print("Creating basis functions...")
    # Class Logging
    class Logging:
        def __init__(self, fileName=""):
            self.fileName = fileName
            self.file = []
            pass
        def openLog(self):
            if self.file==[]:
                self.file = open(self.fileName, "w")
                pass
            pass
        def closeLog(self):
            if self.file!=[]:
                self.file.close()
                pass
            pass
        def addLog(self, newLine):
            if self.file!=[]:
                self.file.write(newLine)
                pass
            pass

    # Class Vertex
    class Vertex():
        def __init__(self, x=0, y=0, z=0):
            self.x = x
            self.y = y
            self.z = z
            pass

    # Class Triangle
    class Triangle():
        def __init__(self, v=[], n=0):
            self.v = v
            self.n = n
            self.adjacent = []
            self.BasisIn = 0
            self.BasisOut = 0
            self.index = 0
            self.filled = []
            pass
        def addAdjacent(self, adjacent):
            self.adjacents.append(adjacent.index)
            pass
        def addBasisIn(self):
            self.BasisIn+=1
            pass
        def addBasisOut(self):
            self.BasisOut+=1
            pass

    # Class Basis
    class Basis:
        def __init__(self):
            pass
        def addBasis(self, file, triS, triD, basisCount):
            if file!=[]:
                flag = 0
                # Scenario 1
                if isVertexEqual(triS.v[1], triD.v[2]) \
                and isVertexEqual(triS.v[2], triD.v[1]):
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 2
                if isVertexEqual(triS.v[1], triD.v[0]) \
                and isVertexEqual(triS.v[2], triD.v[2]):
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 3
                if isVertexEqual(triS.v[1], triD.v[1]) \
                and isVertexEqual(triS.v[2], triD.v[0]):
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 4
                if isVertexEqual(triS.v[2], triD.v[2]) \
                and isVertexEqual(triS.v[0], triD.v[1]):
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 5
                if isVertexEqual(triS.v[2], triD.v[0]) \
                and isVertexEqual(triS.v[0], triD.v[2]):
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 6
                if isVertexEqual(triS.v[2], triD.v[1]) \
                and isVertexEqual(triS.v[0], triD.v[0]):
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 7
                if isVertexEqual(triS.v[0], triD.v[2]) \
                and isVertexEqual(triS.v[1], triD.v[1]):
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 8
                if isVertexEqual(triS.v[0], triD.v[0]) \
                and isVertexEqual(triS.v[1], triD.v[2]):
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 1
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Scenario 9
                if isVertexEqual(triS.v[0], triD.v[1]) \
                and isVertexEqual(triS.v[1], triD.v[0]):
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.v[n].x, triS.v[n].y, triS.v[n].z)
                    file.write(newLine)
                    n = 2
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    n = 0
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.v[n].x, triD.v[n].y, triD.v[n].z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triS.n.x, triS.n.y, triS.n.z)
                    file.write(newLine)
                    newLine = "{:22.14E} {:22.14E} {:22.14E}".format(triD.n.x, triD.n.y, triD.n.z)
                    file.write(newLine)
                    file.write("\n")
                    flag+=1
                    pass
                # Check results
                if flag<1:
                    print("ERROR: Basis {} to {} was not found!".format(triS.index, triD.index))
                    exit(1)
                    pass
                elif flag>1:
                    print("ERROR: More than 1 basis assigned!")
                    exit(1)
                else:
                    return basisCount+1
            pass

    # Check equal vertices
    def isVertexEqual(v1, v2):
        if v1.x==v2.x and v1.y==v2.y and v1.z==v2.z:
            return True
        else:
            return False
        pass

    # Check equal triangles
    def isTriangeEqual(T1, T2):
        if isVertexEqual(T1.v[0], T2.v[0]) and \
           isVertexEqual(T1.v[1], T2.v[1]) and \
           isVertexEqual(T1.v[2], T2.v[2]):
            return True
        else:
            return False
        pass

    # Check adjacent triangles
    def isTriangleAdjacent(triS, triD):
        count = 0
        for vertexS in triS.v:
            for vertexD in triD.v:
                if isVertexEqual(vertexS, vertexD):
                    count+=1
                continue
            continue
        if count==2:
            return True
        else:
            return False

    # Start log file
    if logOption:
        newLogFile = Logging(logFileName)
        newLogFile.openLog()

    # Read triangles file
    triangles = []
    count = 1
    with open(meshFileName, "r") as meshFile:
        for line in meshFile:
            parsedLine = line.split()
            newVertex1 = Vertex(float(parsedLine[0]),
                                float(parsedLine[1]),
                                float(parsedLine[2]))
            newVertex2 = Vertex(float(parsedLine[3]),
                                float(parsedLine[4]),
                                float(parsedLine[5]))
            newVertex3 = Vertex(float(parsedLine[6]),
                                float(parsedLine[7]),
                                float(parsedLine[8]))
            newVertecies = [newVertex1, newVertex2, newVertex3]
            newNormal = Vertex(float(parsedLine[9]),
                               float(parsedLine[10]),
                               float(parsedLine[11]))
            newTriangle = Triangle(newVertecies, newNormal)
            newTriangle.index = count
            triangles.append(newTriangle)
            count+=1
            pass
        pass

    newLine = "{} triangles were found\n".format(len(triangles))
    if logOption:
        newLogFile.addLog(newLine)

    # Find adjacent triangles
    print("Searching for adjacent triangles...")
    from progress.bar import Bar
    with Bar('Progress:', fill='#', suffix='%(percent).1f%% - %(eta)ds', max = len(triangles)) as bar:
        for triS in triangles:
            for triD in triangles:
                if isTriangleAdjacent(triS, triD) and not (isTriangeEqual(triS, triD)):
                    triS.adjacent.append(triD.index)
                    if len(triS.adjacent)>3 or len(triD.adjacent)>3:
                        print("ERROR: More than 3 adjacents were found!")
                        exit(1)
                    pass
                continue
            bar.next()
            continue

    # Log adjacent triangles
    if logOption:
        for tri in triangles:
            for n in tri.adjacent:
                newLine = "Triangle {} is adjacent to {}\n".format(tri.index, n)
                newLogFile.addLog(newLine)
                continue
            continue
    
    # Find a starting basis
    start = triangles[0].index
        
    # Obtain basis functions
    print("Forming basis functions...")
    basisFile = open(basisFileName, "w")
    basisCount = 0
    for i in range(start-1, len(triangles)):
        triS = triangles[i]
        # Check adjacents
        for n in triS.adjacent:
            triD = triangles[n-1]
            if (triD.BasisIn+triD.BasisOut)<len(triD.adjacent) \
                and triS.index not in triD.filled:
                triS.BasisOut+=1
                triD.BasisIn+=1
                triS.filled.append(triD.index)
                triD.filled.append(triS.index)
                newBasis = Basis()
                basisCount = newBasis.addBasis(basisFile, triS, triD, basisCount)
                if logOption:
                    newLine = "{}->{}\n".format(triS.index, triD.index)
                    newLogFile.addLog(newLine)
            continue
    basisFile.close()

    # Final Check
    for tri in triangles:
        if not (tri.BasisIn+tri.BasisOut == len(tri.adjacent) and \
                len(tri.filled) == len(tri.adjacent)):
            print("Basis checking failed!\n")
            return 1
        continue

    print("Basis functions processing was successful!")
    print("{} basis functions have been created".format(basisCount))

    with open("mesh/mesh_info.txt", "w") as file:
        file.write("{}".format(basisCount))
    if logOption:
        newLogFile.closeLog()
    return 0
    
def tryMesh(logOption=False):
    meshFileName = "mesh/shape.stl"
    meshDataFileName = "mesh/mesh_data.dat"
    basisFileName = "mesh/basis.dat"    
    ProcessMesh(meshFileName, meshDataFileName)
    n = creatBasis(meshDataFileName, basisFileName, logOption)
    return n
