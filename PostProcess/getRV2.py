import numpy as np
import subprocess as sp
import os
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

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

def gettingR(x,y):
    yFlip = np.flip(y)
    n = 10 # Precision in the minimum 
    h = 1
    k = 0
    cont = 0
    for i in range(10,len(y)-10):
        h = 1
        for k in range(1,n):
            if y[i-k] > y[i] and y[i+k] > y[i]:
                h = h + 1 
            elif y[i-k] > y[i] and y[i+k] < y[i]:
                break
            elif y[i-k] < y[i] and y[i+k] > y[i]:
                break
            elif y[i-k] < y[i] and y[i+k] < y[i]:
                break
            if h == 7:
                print("Local minimum found")
                xMin = x[i]
                yMin = y[i]
                print(xMin, yMin)
                cont = 1
                break    
        if cont == 1:
            break
    if xMin is not None:
        xMinima = xMin
        yMinima = yMin
    else:
        print("No minimum value found.")
        xMinima = None
        xMinima = None
    return xMinima, yMinima

def getPlot(x,y,xMin,yMin):
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.plot(xMin,yMin,'ro')
    plt.show()
    fig, ax = plt.subplots()
    ax.spines["top"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.set_title('Minumum distance vs. t0',fontsize = 15)
    ax.set_ylabel('Radius',fontsize = 15)
    ax.set_xlabel('$t/t_{\gamma}$',fontsize = 15)
    ax.scatter(time,R,s=3,marker='o',color='black')
    plt.show()

    data = np.stack((xR, R), axis=1)
    Rot= np.array([[0,1],[-1,0]])
    data = np.dot(Rot,data.T).T
    xRrot = np.array(data[:,0])
    Rrot = np.array(data[:,1])
    fig, ax = plt.subplots()
    ax.spines["top"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.set_title('Minimum distance vs. x position',fontsize = 15)
    ax.set_ylabel('Radius',fontsize = 15)
    ax.set_xlabel('x position',fontsize = 15)
    ax.scatter(xRrot,Rrot,s=3,marker='o',color='black')
    plt.show()

    dataa = np.stack((xR, R,xRrot,Rrot), axis=1)
    with open('Rdata.txt', 'w') as f:
        for row in dataa:
            row_data = '\t'.join([str(i) for i in row]) # tab separated values
            f.write(row_data + '\n')

existing_path = "intermediate"
path_to_folder = os.path.join(os.getcwd(), existing_path)
maxFiles = len([f for f in os.listdir(path_to_folder) if os.path.isfile(os.path.join(path_to_folder, f))])
dm = 5.0e-4
m = int(maxFiles/dm*5.0e-4)
print(m)
R = np.empty(m)
xR = np.empty(m)
time = np.empty(m)

for ti in range(m):
    t = ti*(dm)
    time[ti] = t
    print(f"Time is {t}")
    place = f"intermediate/snapshot-{t:5.4f}"
    x, y = getSegs(place)
    xMin, yMin = gettingR(x,y)
    R[ti] = yMin
    xR[ti] = xMin
getPlot(x,y,xMin,yMin)