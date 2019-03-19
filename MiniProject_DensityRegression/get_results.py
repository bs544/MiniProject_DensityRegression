from MiniProject_DensityRegression.network import NetworkHandler
from MiniProject_DensityRegression.fingerprints import fingerprints
from MiniProject_DensityRegression.castep_density import Castep_density
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

s_ = 15
colour = ['b','m','k','c','g','r','y']

trainsets = {"InnerEnsemble":[0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5],"debug_decay":[0.5,0.75,1.0,1.25],"AlternatingEnsemble":[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5],"AlternatingEnsemblie_big":[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5]}

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
    x_diff = x[1]-x[0]
    y_diff = y[1]-y[0]
    #print(x_diff,y_diff)
    x_wrap = x-xmax+shift[0]-x_diff
    y_wrap = y-ymax+shift[1]-y_diff
    x = x+shift[0]
    y = y+shift[1]
    x = np.where(x>xmax,x_wrap,x)
    y = np.where(y>ymax,y_wrap,y)
    return x,y

def results_along_a_line(netName,cellName,pltline=None,ensembleplot=True,wDensity=True,save=False):
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
            plt.plot(x[:-1],density[:-1],'r--')
            plt.legend(["Density Correction","Ensemble Prediction"])
        plt.xlabel('x/$\AA$')
        plt.ylabel("Density Correction e/$\AA ^3$")
        plt.show()
    else:
        mean,std = C.NetHandler.predict(line_FPs,ensemble=False,standard=False)
        print(mean.shape)
        # mean = np.mean(mean,axis=1)
        # plt.plot(x[:-1],mean[:-1])
        nmeans=mean.shape[1]
        leg = []
        for i in range(nmeans):
            plt.plot(x[:-1],mean[:-1,i],colour[i])
            leg.append("Net_{}".format(i))
        if (wDensity):
            plt.plot(x[:-1],density[:-1],'r--')
            leg.append("Density Correction")
        plt.legend(leg)
        plt.xlabel('x/$\AA$')
        plt.ylabel("Density Correction e/$\AA ^3$")

        plt.show()
    plt.close()
    return

def std_2d(netName,cellName,plane=None,rel_std=True,rel_cap=2.0,save=False):
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
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

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

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1)
    plt.colorbar(sc)
    plt.show()

    return

def FPmap(netName,cellName,plane=None,save=False):
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
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]
    # np.save("FPs.npy",plane_FPs)

    plane_FPs -= np.mean(plane_FPs,axis=0)

    std = np.std(plane_FPs,axis=0)
    plane_FPs = np.where(std[None,:]>1e-5,plane_FPs/10*std[None,:],plane_FPs)
    #print(np.argmax(plane_FPs,axis=0))

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z = np.linalg.norm(plane_FPs,axis=1)

    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1)#,vmax=0.1)
    plt.colorbar(sc)
    plt.show()

def densityMap(netName,cellName,plane=None,save=False):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=2.5)

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
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

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

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1)
    plt.colorbar(sc)
    plt.show()
    plt.close()
    return

def NetInputMap(netName,cellName,plane=None,save=False):
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
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]

    mean = N.datamean
    std = N.datastd

    plane_FPs -= mean#np.mean(plane_FPs,axis=0)

    #std = np.std(plane_FPs,axis=0)
    plane_FPs = np.where(std[None,:]>1e-5,plane_FPs/10*std[None,:],plane_FPs)
    #print(np.argmax(plane_FPs,axis=0))

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z = np.linalg.norm(plane_FPs,axis=1)

    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1)#,vmax=0.1)
    plt.colorbar(sc)
    plt.show()
    return

def NetErrorMap(netname,cellname,plane=None,relative=True,save=False):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=3,nmax=3,r_c=5.5)

    #setup and load network
    N = NetworkHandler(fp,name=netname)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellname)

    grid = C.supercell.grid

    if (plane is None):
        plane=grid[1,2]
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]
    # allmean,std = C.ensemble_predict(C.supercell.FP)

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z, std_ = C.ensemble_predict(plane_FPs)

    filename = "{}/{}".format("FP_data",cellname)
    f = open(filename,'rb')
    dict = pickle.load(f)
    f.close()
    density=dict["density"]
    z_ = density[plane_idx]

    if (relative):
        z = np.abs((z-z_)/std_)
        #z = np.abs((z-z_)/z_)
        maxerr = np.max(z)
        vmax = min([3,maxerr])
    else:
        z = np.abs(z-z_)
        vmax = np.max(z)

    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1,vmax=vmax)
    plt.colorbar(sc,label='|Prediction Error/Uncertainty|')
    plt.ylabel("y ($\AA$)")
    plt.xlabel("x ($\AA$)")

    plt.show()
    return

def ErrorStdPlot(netname,cellname,save=False):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=3,nmax=3,r_c=5.5)

    #setup and load network
    N = NetworkHandler(fp,name=netname)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    C.get_cell_data(cellname)

    grid = C.supercell.grid

    FPs = C.supercell.FP
    # allmean,std = C.ensemble_predict(C.supercell.FP)


    z, std_ = C.ensemble_predict(FPs)

    filename = "{}/{}".format("FP_data",cellname)
    f = open(filename,'rb')
    dict = pickle.load(f)
    f.close()
    density=dict["density"]
    z_ = density
    print(np.max(z_))

    err = np.abs(z-z_)

    plt.plot(std_,err,'b.')
    plt.xlabel("Ensemble Predicted Variance")
    plt.ylabel("Absolue Ensemble Error")
    plt.show()
    return

def NetOutputMap(netName,cellName,plane=None,save=False):
    #using precalculated fingerprints and cell data for everything

    #define descriptor (optional)
    fp = fingerprints(lmax=3,nmax=3,r_c=5.5)

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
    else:
        idx = np.argmin(plane-grid[:,2])
        plane = grid[idx,2]

    plane_idx = np.where(plane==grid[:,2])[0]

    plane_FPs = C.supercell.FP[plane_idx,:]
    allmean,_ = C.ensemble_predict(C.supercell.FP)

    x = C.supercell.grid[plane_idx,0]
    y = C.supercell.grid[plane_idx,1]

    x,y = shiftxy(x,y)

    z, _ = C.ensemble_predict(plane_FPs)


    #X,Y,Z = transform2d(x,y,z)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x,y,c=z,cmap=cm,marker=',',s=s_,alpha=1)
    plt.colorbar(sc)
    if (save):
        plt.savefig("{}-{}_Out.pdf".format(netName,cellName))
    plt.show()
    return

def VarWithR(netName,FP_dir,Cell_dir):
    if(not os.path.isdir(FP_dir) or not os.path.isdir(Cell_dir)):
        print("invalid fingerrpint or cell directory")
        return
    #using precalculated fingerprints and cell data for everything


    if (not os.path.isdir("{}_data".format(netName))):
        os.mkdir("{}_data".format(netName))
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
        mae=[]
        for file in files:
            C.get_cell_data(file,include_density=True)
            mean,std = C.ensemble_predict()
            density = C.supercell.train_density
            r.append(np.linalg.norm(C.supercell.cart_coords[0,:]-C.supercell.cart_coords[1,:]))
            if (len(density)==len(mean)):
                rmse_ = np.sqrt(np.mean(np.square(density-mean)))
                rmse.append(rmse_)
                mae_ = np.mean(np.abs(density-mean))
                mae.append(mae_)
            else:
                print("Incompatible mean and density")
            H.append(C.getH(std))

        H = np.asarray(H)
        rmse = np.asarray(rmse)
        r = np.asarray(r)
        mae = np.asarray(mae)

        np.save("{}_data/r.npy".format(netName),r)
        np.save("{}_data/rmse.npy".format(netName),rmse)
        np.save("{}_data/H.npy".format(netName),H)
        np.save("{}_data/mae.npy".format(netName),mae)

    else:
        files = os.listdir(FP_dir)
        train_idx = []
        test_idx = []
        r = np.load("{}_data/r.npy".format(netName))
        rmse = np.load("{}_data/rmse.npy".format(netName))
        H = np.load("{}_data/H.npy".format(netName))
        mae = np.load("{}_data/mae.npy".format(netName))
        for i, file in enumerate(files):
            if (any((abs(trainsets[netName]-r[i]))<1e-5)):
                train_idx.append(i)
            else:
                test_idx.append(i)


    test_H = H[test_idx]
    train_H = H[train_idx]
    test_r = r[test_idx]
    train_r = r[train_idx]
    train_rmse = rmse[train_idx]
    test_rmse = rmse[test_idx]
    train_mae = mae[train_idx]
    test_mae = mae[test_idx]


    plt.plot(train_r,train_H,'bx')
    plt.plot(test_r,test_H,'rx')
    plt.legend(["Train set","Test set"])
    plt.xlabel("Interatomic Distance ($\AA$)")
    plt.ylabel("H")
    plt.show()
    plt.close()
    plt.plot(train_r,train_mae,'bx')
    plt.plot(test_r,test_mae,'rx')
    plt.legend(["Train set","Test set"])
    plt.title("Ensemble 2")
    plt.xlabel("Interatomic Distance ($\AA$)")
    plt.ylabel("Mean Absolue Error (e/$\AA^3$)")
    plt.show()
    plt.close()

    plt.plot(train_mae,train_H,'bx')
    plt.plot(test_mae,test_H,'rx')
    plt.legend(["Train set","Test set"])
    plt.ylabel("H")
    plt.xlabel("Mean Absolute Error (e/$\AA^3$)")
    plt.show()
    plt.close()

    return

def TotErrAndH(netNames):
    #loads all of the numpy files with the MAE and H values for the networks specified and plots them all
    mae = []
    H = []
    leg_convert = {"InnerEnsemble":"Ensemble 2","debug_decay":"Ensemble 3","AlternatingEnsemble":"Ensemble 1"}
    leg = []
    for name in netNames:
        dir = "{}_data".format(name)
        if (os.path.isdir(dir)):
            errarray = "{}/mae.npy".format(dir)
            harray = "{}/H.npy".format(dir)
            H_ = np.load(harray)
            mae_ = np.load(errarray) * 1000
        else:
            print("network data not present")
        plt.plot(mae_,H_,'x')
        en_name = leg_convert[name]
        leg.append(en_name)
    plt.legend(leg)
    plt.xlabel("Mean Absolute Error of Cell (meV/$\AA^3$)")
    plt.ylabel("Global Uncertainty, H")
    plt.show()
    return

def dSCFplots(netNames,SCFs):
    H = []
    r = []
    dSCF = []
    leg = []
    min = 0
    max = 0
    linepoint = -6.1
    leg_convert = {"InnerEnsemble":"Ensemble 2","debug_decay":"Ensemble 3","AlternatingEnsemble":"Ensemble 1"}
    for name in netNames:
        dSCF_ = SCFs[name]-SCFs["Unchanged"]
        dir = "{}_data".format(name)
        if (not os.path.isdir(dir)):
            print("No directory found")
            return
        harray = "{}/H.npy".format(dir)
        rarray = "{}/r.npy".format(dir)
        H_ = np.load(harray)
        r_ = np.load(rarray)
        idx = np.argsort(r_)
        r_ = r_[idx]
        H_ = H_[idx]
        plt.plot(H_[:len(dSCF_)],dSCF_,'x')
        H.append(H_[:len(dSCF_)])
        r.append(r_[:len(dSCF_)])
        leg.append(leg_convert[name])
        dSCF.append(dSCF_)
        max_ = np.max(dSCF_)
        min_ = np.min(dSCF_)
        if (max_>max):
            max = max_
        if (min_<min):
            min = min_

    y = np.linspace(min-1,max+1,2)
    x = np.asarray([linepoint,linepoint])
    plt.plot(x,y,'r--')
    leg.append("$H_{cut}$")
    plt.legend(leg,loc=4)
    plt.xlabel("Global Uncertainty, $H$")
    plt.ylabel("dSCF")
    plt.ylim([min-1,max+1])
    plt.show()
    return


def NetDensityWrite(netName,FP_dir,Cell_dir,files=None):
    if(not os.path.isdir(FP_dir) or not os.path.isdir(Cell_dir)):
        print("invalid fingerrpint or cell directory")
        return
    #using precalculated fingerprints and cell data for everything


    if (not os.path.isdir("{}_data".format(netName))):
        os.mkdir("{}_data".format(netName))
    #define descriptor (optional)
    fp = fingerprints(lmax=4,nmax=5,r_c=4.5)

    #setup and load network
    N = NetworkHandler(fp,name=netName)
    N.load()

    #setup the castep interface (don't really need it here, but it has everything conveniently stored)
    C = Castep_density(fp,N)
    #check to see the network was properly loaded
    C.setupNetwork()

    if (files is None):
        files = os.listdir(FP_dir)
    H = []
    rmse = []
    r=[]
    mae=[]
    for i, file in enumerate(files):
        if ("0.5" in file):
            den_file = "{}_data/{}{}".format(netName,file[:-4],"initial_den")
            print("Writing density to: ", den_file)
            C.get_cell_data(file,include_density=True)
            # print(np.max(C.supercell.fin_density))
            # mean,_ = C.ensemble_predict()
            # print(np.mean(mean))
            C.setCellDensities(den_file,taper=False)

    return

def NetLossPlot(netName,save=False):
    fp = fingerprints()
    N = NetworkHandler(fp,name=netName)
    N.load()
    x = np.linspace(0.0,N.nEpochs,len(N.loss[0]))
    for i, loss in enumerate(N.loss):
        plt.plot(x,loss)
    #plt.show()
    plt.close()
    x_ = np.linspace(0.0,N.nEpochs,len(N.test_rmse[0]))
    leg = []

    for i, rmse in enumerate(N.test_rmse):
        print(len(rmse))
        plt.plot(x_,rmse,colour[i])
        leg.append("Net_{}".format(i))
    plt.legend(leg)
    if (save):
        plt.savefig(NetRMSE.pdf)
    plt.show()
    plt.close()
    return



#results_along_a_line("smallInnerEnsemble","Hcell1.5.pckl",pltline=np.asarray([0.1,0.0]))
#std_2d("DensityEnsemble","Hcell1.41.pckl",plane=0.1)
