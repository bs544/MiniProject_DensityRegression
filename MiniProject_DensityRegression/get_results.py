from network import NetworkHandler
from fingerprints import fingerprints
from castep_density import Castep_density
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

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

def shiftxy(x,y,shift=[5.0,5.0]):
    xmax=max(x)
    ymax=max(y)
    x_wrap = x-xmax+shift[0]
    y_wrap = y-ymax+shift[1]
    x = x+shift[0]
    y = y+shift[1]
    x = np.where(x>xmax,x_wrap,x)
    y = np.where(y>ymax,y_wrap,y)
    return x,y

def results_along_a_line(netName,cellName,pltline,ensembleplot=True,wDensity=True):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName,include_density=True)
    print(np.mean(C.supercell.train_density))

    grid = C.supercell.grid

    if (pltline is None):
        pltline=np.zeros(2)
        pltline[0]=grid[1,0]
        pltline[1]=grid[2,0]#assumes spherical cell

    line_idx = np.where(pltline[0]==grid[:,1])[0]
    line_idx = np.where(pltline[1]==grid[line_idx,2])[0]

    line_FPs = C.supercell.FP[line_idx,:]

    x = C.supercell.grid[line_idx,0]

    if (wDensity):
        #density stored with fp data
        density = C.supercell.train_density[line_idx]

    if (ensembleplot):
        ensemble_mean,ensemble_std = C.ensemble_predict(line_FPs)
        plt.errorbar(x,ensemble_mean,yerr=ensemble_std,fmt='o')
        if (wDensity):
            plt.plot(x[:-1],density[:-1],'r')
            plt.legend(["Density Correction","Ensemble Prediction"])
        plt.show()
    else:
        mean,std = C.NetHandler.predict(line_FPs,ensemble=False)
        nmeans=mean.shape[1]
        for i in range(nmeans):
            plt.errorbar(x[:-1],mean[:-1,i],yerr=std[:-1,i])
        plt.show()
    plt.close()
    return

def std_2d(netName,cellName,plane,rel_std=True,rel_cap=2.0):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName)

    grid = C.supercell.grid

    if (plane is None):
        plane=grid[1,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    _, z = C.ensemble_predict(plane_FPs)
    if (rel_std):
        z =np.abs(z)/np.abs(_)
        if(rel_cap is not None):
            z = np.where(z>rel_cap,rel_cap,z)


    #X,Y,Z = transform2d(x,y,z)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker='.',s=5.0,alpha=1)
    plt.colorbar(sc)
    plt.show()

    return

def FPmap(netName,cellName,plane=None):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName)

    grid = C.supercell.grid

    if (plane is None):
        plane=grid[1,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z = np.linalg.norm(plane_FPs,axis=1)

    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker='.',s=5.0,alpha=1)
    plt.colorbar(sc)
    plt.show()

def densityMap(netName,cellName,plane=None):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    # N = NetworkHandler(fp,name=netName)
    # N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp)
    #check to see the network was properly loaded
    #C.setupNetwork()

    C.get_cell_data(cellName,include_density=True)

    grid = C.supercell.grid

    if (plane is None):
        plane=grid[1,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    if (True):
        filename = "{}/{}".format("FP_data",cellName)
        print(filename)
        f = open(filename,'rb')
        dict = pickle.load(f)
        f.close()
        density=dict["density"]
        z = density[plane_idx]
    else:
        z = C.supercell.fin_density
        z = z[plane_idx]

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker='.',s=5.0,alpha=1)
    plt.colorbar(sc)
    plt.show()
    plt.close()
    return

def NetOutputMap(netName,cellName,plane=None):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellName)

    grid = C.supercell.grid

    if (plane is None):
        plane=grid[1,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]
    allmean,_ = C.ensemble_predict(C.supercell.FP)

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z, _ = C.ensemble_predict(plane_FPs)


    #X,Y,Z = transform2d(x,y,z)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker='.',s=5.0,alpha=1)
    plt.colorbar(sc)
    plt.show()
    return

def VarWithR(netName,FP_dir,Cell_dir):
    if(not os.path.isdir(FP_dir) or not os.path.isdir(Cell_dir)):
        print("invalid fingerrpint or cell directory")
        return
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    files = os.listdir(FP_dir)
    H = []
    rmse = []
    r=[]
    for file in files:
        C.get_cell_data(file,include_density=True)
        mean,std = C.ensemble_predict()
        density = C.supercell.train_density
        r.append(np.linalg.norm(C.supercell.cart_coords[0,:]-C.supercell.cart_coords[1,:]))
        if (len(density)==len(mean)):
            rmse_ = np.sqrt(np.mean(np.square(density-mean)))
            rmse.append(rmse_)
        else:
            print("Incompatible mean and density")
        H.append(C.getH(std))

    plt.plot(r,H,'b.')
    plt.plot(r,rmse,'r.')
    plt.legend(["H","RMSE"])
    plt.show()
    plt.close()

    plt.plot(H,rmse)
    plt.show()
    plt.close()

    return

#results_along_a_line("smallInnerEnsemble","Hcell1.5.pckl",pltline=np.asarray([0.1,0.0]))
#std_2d("DensityEnsemble","Hcell1.41.pckl",plane=0.1)
