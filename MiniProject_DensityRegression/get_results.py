from network import NetworkHandler
from fingerprints import fingerprints
from castep_density import Castep_density
import matplotlib.pyplot as plt
import numpy as np
import pickle

def transform2d(x,y,z):
    xlen = max(x)
    ylen = max(y)

    xint = x[31]-x[0]
    yint = y[31]-y[0]

    xsize = int((xlen/xint)+0.0001)
    ysize = int((ylen/yint)+0.0001)

    X = np.linspace(0.0,xlen,xsize+1)
    Y = np.linspace(0.0,ylen,ysize+1)

    Z = np.zeros((len(X),len(Y)))

    indx = (x*xsize/xlen + 0.0001).astype(int)
    indy = (y*ysize/ylen + 0.0001).astype(int)

    Z[indx,indy] = z

    #print(Z)

    return X,Y,Z

def results_along_a_line(netName,cellName,pltline,ensembleplot=True,wDensity=True):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.0)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName)

    grid = C.supercell.grid

    line_idx = np.where(pltline[0]==grid[:,1])[0]
    line_idx = np.where(pltline[1]==grid[line_idx,2])[0]

    line_FPs = C.supercell.FP[line_idx,:]

    x = C.supercell.grid[line_idx,0]

    if (wDensity):
        #density stored with fp data
        f=open("./FP_data/{}".format(cellName),'rb')
        dict = pickle.load(f)
        f.close()
        density = dict["density"][line_idx]

    if (ensembleplot):
        ensemble_mean,ensemble_std = C.ensemble_predict(line_FPs)
        plt.errorbar(x,ensemble_mean,yerr=ensemble_std,fmt='o')
        if (wDensity):
            plt.plot(x,density)
        plt.show()
    else:
        mean,std = C.NetHandler.predict(line_FPs,ensemble=False)
        nmeans=mean.shape[1]
        for i in range(nmeans):
            plt.errorbar(x,mean[:,i],yerr=std[:,i])
        plt.show()
    return

def std_2d(netName,cellName,plane):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.0)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName)

    grid = C.supercell.grid

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    _, z = C.ensemble_predict(plane_FPs)

    #X,Y,Z = transform2d(x,y,z)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=01.0,alpha=1)
    plt.colorbar(sc)
    plt.show()

def ErrorVsR(netName,FP_dir,Cell_dir):
    if(not os.path.isdir(FP_dir) or not os.path.isdir(Cell_dir)):
        print("invalid fingerrpint or cell directory")
        return
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.0)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    files = os.path.listdir(FP_dir)
    for file in files:
        C.get_cell_data(file)
        mean,std = C.ensemble_predict()
        H =

#results_along_a_line("DensityEnsemble","Hcell1.41.pckl",pltline=np.asarray([0.1,0.0]))
#std_2d("DensityEnsemble","Hcell1.41.pckl",plane=0.1)
