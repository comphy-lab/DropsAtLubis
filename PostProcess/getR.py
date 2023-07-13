'''
Title: getR.py 
Function: Finds the minimum radius at the interface. i iterates each value in Facet
until lowest value is found for which both sides increase aka minimum. 

'''

# Import required packages 
import math
import numpy as np
import subprocess as sp
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

# Set plot parameters to LaTeX
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['font.size'] = 12

# Functions 
def gettingFacets(filename):
    # Execute the "getFacet" program and capture its output
    exe = ["./getFacet", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2,z2)))
                    segs.append(((-r1, z1),(-r2,z2)))
                    skip = True
    return segs

def getSegs(place):
    segs = gettingFacets(place)
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

def gettingR(x,y,t):
    n = 10
    idx = y > 0.01
    y = y[idx]
    x = x[idx]
    idx2 = x < 0.01
    y = y[idx2]
    x = x[idx2]
    xIndex = np.zeros(len(x))
    yIndex = np.zeros(len(y))
    yMinimum = np.zeros(len(y))
    xMinimum = np.zeros(len(x))
    for i in range(30,len(y)-15):
        h = 0
        distance = abs(np.sqrt((x[i]-x[i+1])**2 + (y[i]-y[i+1])**2))
        if distance > 0.4:
            break
        for k in range(n):
            if y[i-k] > y[i] and y[i+k] > y[i]:
                h = h + 1 
            elif y[i-k] > y[i] and y[i+k] < y[i]:
                break
            elif y[i-k] < y[i] and y[i+k] > y[i]:
                break
            elif y[i-k] < y[i] and y[i+k] < y[i]:
                break
            if h == 4:
                xIndex[i] = 1
                yIndex[i] = 1
                xMinimum[i] = x[i]
                yMinimum[i] = y[i]
                break
    if len(yMinimum) == 0: 
        print("No Minimum found at time %f"%t)
        xMin = 0
        yMin = 0
    else:
        boolArray = xIndex.astype(bool)
        if all(not val for val in boolArray):
            print("No minimum found")
            xMin = 0
            yMin = 0 
            return xMin, yMin
        nonzeroYMinimum = yMinimum[boolArray]
        nonzeroxMinimum = xMinimum[boolArray]
        globalMinimumIdx = nonzeroYMinimum.argmin()
        xMin = nonzeroxMinimum[globalMinimumIdx]
        yMin = nonzeroYMinimum[globalMinimumIdx]
        delta = abs(xMin - x[np.where(x == xMin)[0][0]+15])  
        print(xMin,x[np.where(x == xMin)[0][0]+15])
        if delta > 0.15:
            print("No minimum found")
            xMin = 0
            yMin = 0
            return xMin, yMin
        print("Minimum found at x = %f and y = %f"%(xMin,yMin))
    return xMin, yMin
    
def getPlot(x,y,xMin,yMin,time):
    # fig, ax = plt.subplots()
    # ax.scatter(x,y)
    # ax.plot(xMin,yMin,'ro')
    # plt.show()
    # fig, ax = plt.subplots()
    # time = [time[i] for i in range(len(yMin)) if yMin[i] !=0]
    # R = [yMin[i] for i in range(len(yMin)) if yMin[i] !=0]
    # ax.spines["top"].set_color("none")
    # ax.spines["right"].set_color("none")
    # ax.set_title('Minumum distance vs. t0',fontsize = 15)
    # ax.set_ylabel('Radius',fontsize = 15)
    # ax.set_xlabel('$t/t_{\gamma}$',fontsize = 15)
    # ax.scatter(time,R,s=3,marker='o',color='black')
    # plt.show()
    data = np.stack((xMin, yMin), axis=1)
    Rot= np.array([[0,1],[-1,0]])
    data = np.dot(Rot,data.T).T
    xRrot = np.array(data[:,0])
    Rrot = np.array(data[:,1])
    # fig, ax = plt.subplots()
    # ax.spines["top"].set_color("none")
    # ax.spines["right"].set_color("none")
    # ax.set_title('Minimum distance vs. x position',fontsize = 15)
    # ax.set_ylabel('Radius',fontsize = 15)
    # ax.set_xlabel('x position',fontsize = 15)
    # ax.scatter(xRrot,Rrot,s=3,marker='o',color='black')
    # plt.show()
    dataa = np.stack((time, R,xRrot,Rrot), axis=1)
    with open('Rdata.txt', 'w') as f:
        for row in dataa:
            row_data = '\t'.join([str(i) for i in row]) # tab separated values
            f.write(row_data + '\n')

existingPath = "intermediate"
pathToFolder = os.path.join(os.getcwd(), existingPath)
maxFiles = len([f for f in os.listdir(pathToFolder) if os.path.isfile(os.path.join(pathToFolder, f))])
dm = 5.0e-3
m = int(maxFiles/dm*5.0e-4)
R = np.empty(m)
xR = np.empty(m)
time = np.empty(m)

for ti in range(m):
    t = ti*(dm)
    time[ti] = t
    print(f"Time is {t}")
    place = f"intermediate/snapshot-{t:5.4f}"
    x, y = getSegs(place)
    xMin, yMin = gettingR(x,y,t)
    R[ti] = yMin
    xR[ti] = xMin
getPlot(x,y,xR,R,time)
