# Import required packages 
import math
import numpy as np
import subprocess as sp
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import scipy.io as sio


# Functions 

EPSILON = 1e-5
SKIP_THRESHOLD = 1e2

def execute_command(command):
    """Executes a shell command and returns stdout and stderr."""
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def gettingFacets(filename, Tracer):
    """Function to get facets from a given filename and Tracer value."""
    exe = ["./getFacet1", filename] if Tracer == 1 else ["./getFacet2", filename]
    _, stderr = execute_command(exe)
    temp2 = stderr.split("\n")
    temp2 = [line for line in temp2 if line.strip()]  # Remove empty lines
    segs = []

    if len(temp2) > SKIP_THRESHOLD:
        for n1 in range(0, len(temp2) - 1, 2):  # Process every other line
            temp3 = temp2[n1].split()
            temp4 = temp2[n1+1].split()
            r1, z1 = float(temp3[1]), float(temp3[0])
            r2, z2 = float(temp4[1]), float(temp4[0])
            segs.append(((r1, z1),(r2, z2)))
            segs.append(((-r1, z1),(-r2, z2)))

    return segs

def getSegs(place,a):
    segs = gettingFacets(place,a)
    if (len(segs) == 0):
        print("Problem in the available file %s" % place)
    x_coords = np.array([seg[0][0] for seg in segs])
    y_coords = np.array([seg[0][1] for seg in segs])
    data = np.stack((x_coords, y_coords), axis=1)
    dataPositive = data[data[:,0]>0]
    R = np.array([[0,-1],[1,0]])
    data = np.dot(R,dataPositive.T).T
    x = np.array(data[:,0])
    y = np.array(data[:,1])
    sorted_indices = np.argsort(x)
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]
    return x_sorted,y_sorted 
def getPlot(x,y):
    data = np.stack((x, y), axis=1)
    adict = {}
    adict['data'] = data
    sio.savemat('/media/vatsal/Spreading/TestPrecursorFilm/RegimeMap simulations/3044/3304.mat',adict)
    # with open('saveFacet3044.txt', 'w') as f:
    #     for row in data:
    #         row_data = '\t'.join([str(i) for i in row]) # tab separated values
    #         f.write(row_data + '\n')
existingPath = "intermediate"
pathToFolder = os.path.join(os.getcwd(), existingPath)
maxFiles = len([f for f in os.listdir(pathToFolder) if os.path.isfile(os.path.join(pathToFolder, f))])
dm = 5.0e-2
m = int(3.5/dm)
time = np.empty(m)
n = 2000
xData = np.empty((m,n))
yData = np.empty((m,n))
print(len(xData))
for ti in range(m):
    t = ti*(dm)
    time[ti] = t
    print(f"Time is {t}")
    place = f"intermediate/snapshot-{t:5.4f}"
    x, y = getSegs(place,1)
    x = np.pad(x,(0,n-len(x)))
    y = np.pad(y,(0,n-len(y)))
    xData[ti] = x
    yData[ti] = y
getPlot(xData,yData)




