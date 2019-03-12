from parsers.structure_class import supercell
from fingerprints import fingerprints
from network import NetworkHandler
import numpy as np
from scipy import stats
import pickle
from parsers.dft.parser_castep import wrap_inhouse
import os
from operator import mul
from functools import reduce
import time

class Castep_density():
    def __init__(self,descriptor=None,network_handler=None,calc_FP=False,trainNet=False,data_dir='./data/',train_dir='../CastepCalculations/DenNCells/'):
        self.traindir = train_dir
        self.datadir = data_dir
        self.trainNet = trainNet
        self.calc_FP = calc_FP
        self.supercell = None
        self.xcut = 6.7
        self.scale = 0.2
        self.set_descriptor(descriptor)
        self.set_network_handler(network_handler)

    def set_network_handler(self,network_handler):
        if (isinstance(network_handler,NetworkHandler)):
            self.NetHandler = network_handler
        else:
            self.NetHandler = NetworkHandler(self.descriptor,train_dir=self.traindir)

    def set_descriptor(self,descriptor):

        if (isinstance(descriptor,fingerprints)):
            self.descriptor = descriptor

        else:
            print("descriptor passed was invalid, creating default fingerprints class")
            self.descriptor = fingerprints()
        self.bilength = self.descriptor.bilength
        self.powerlength = self.descriptor.powerlength

        return

    def setupNetwork(self):
        #check if Network is set up
        #if not check if the network can be loaded
        #if not try to train one
        if(self.NetHandler.trained or self.NetHandler.loaded):
            return
        else:
            self.NetHandler.load()
            if (not self.NetHandler.loaded and self.trainNet):
                self.NetHandler.get_data()
                self.NetHandler.train()
            elif (not self.NetHandler.loaded and not self.trainNet):
                print("Warning: Network can't be loaded and hasn't been trained results will be bad")
            return


    def get_frac_coords(self,at_posns,cell):
        #coordinates are absolute cartesians, want to put them in as fractional
        #at_posns[atom_idx,cartesian_element]
        frac_posns = np.zeros((at_posns.shape))
        inv_cell = np.linalg.inv(cell)
        for i in range(at_posns.shape[0]):
            frac_posns[i,:] = np.dot(inv_cell,at_posns[i,:])
        return frac_posns

    def taper(self,x):
        x_prime = (self.xcut - x)/self.scale
        #zeros = np.zeros((x_prime.shape))
        x_prime4 = (np.where(x_prime>=0,x_prime,0))**4
        taper = x_prime4/(1+x_prime4)
        return taper

    def get_cell_data(self,filename,cell_dir="./Cell_data/",fp_dir="./FP_data/",include_density=False):
        cell_keys=["cell","at_posns","grid","fin_density"]
        fp_keys=["fingerprints","density"]
        if (os.path.isdir(cell_dir) and os.path.isfile("{}{}".format(cell_dir,filename))):
            f = open("{}{}".format(cell_dir,filename),'rb')
            cell_dict = pickle.load(f)
            f.close()
            if (all([key in cell_dict.keys() for key in cell_keys])):
                print("cell data loaded")
            else:
                print("cell data incomplete")
            if (os.path.isdir(fp_dir) and os.path.isfile("{}{}".format(fp_dir,filename))):
                f = open("{}{}".format(fp_dir,filename),'rb')
                fp_dict = pickle.load(f)
                f.close()
                if (all([key in fp_dict.keys() for key in fp_keys])):
                    print("fp data loaded")
                else:
                    print("fp data incomplete")

                self.set_supercell(cell_dict["cell"],cell_dict["at_posns"],cell_dict["grid"],fp_dict["fingerprints"])
                if(include_density):
                    self.supercell.train_density=fp_dict["density"]
                    self.supercell.fin_density=cell_dict["fin_density"]

            elif(self.calc_FP):
                print("not implemented yet")
            else:
                print("fp data not loaded")

        elif (os.path.isdir(self.datadir)):
            print("not implemented yet")

        else:
            print("no cell data, aborting load")

        return

    def set_supercell(self,cell,at_posns,grid,FP):
        #use Andrew's supercell class as the storage for calculations on each unit cell
        self.supercell = supercell()
        self.supercell.set_cell(cell)
        frac_coords = self.get_frac_coords(at_posns,cell)
        self.supercell.cart_coords=at_posns
        self.supercell.set_positions(frac_coords)
        at_species = ['H' for i in range(at_posns.shape[0])]
        self.supercell.set_species(at_species)
        if(len(grid) == FP.shape[0]):
            self.supercell.grid = grid
            self.supercell.FP = FP
        else:
            print("incompatible grid and fingerprints")
        return

    def ensemble_predict(self,X=None):

        if (X is None):
            ensemble_mean,ensemble_std = self.NetHandler.predict(self.supercell.FP,standard=False)

        else:
            ensemble_mean,ensemble_std = self.NetHandler.predict(X,standard=False)

        return ensemble_mean,ensemble_std

    def get_full_grid(self):
        #figures out the interval between grid points and creates a full grid for the cell
        #this only works for a cuboid cell

        #calculation of the interval assumes that grid points are far more likely to be in a clump than isolated
        diff = self.supercell.grid[1:,:]-self.supercell.grid[:-1,:]
        interval = []# stats.mode(diff,axis=0)[0]
        for i in range(3):
            nonzero = diff[np.nonzero(diff[:,i]),i]
            inter = stats.mode(nonzero,axis=1)[0]
            interval.append(inter[0][0])
        print("interval: ",interval)

        #check to see if the cell vectors are divisible by the integer
        num_ints = [self.supercell.cell[i,i]/interval[i] for i in range(3)]
        grid_lines = []
        for i in range(3):
            if (num_ints[i]-round(num_ints[i])<0.0001):
                num_ints[i] = round(num_ints[i])
                grid_lines.append(np.linspace(0.0,self.supercell.cell[i,i]-interval[i],num_ints[i]))
            else:
                print("invalid interval")

        mesh = np.meshgrid(grid_lines[0],grid_lines[1],grid_lines[2])
        newmesh = []
        for i in range(3):
            newmesh.append(np.asarray(mesh[i]).flatten())
        fullgrid = np.asarray(newmesh)
        return fullgrid

    def getH(self,std):
        H = np.mean(np.log(std))
        return H

    def getnonzero_index(self,big_grid,small_grid):
        shape = big_grid.shape
        grid_len = int(np.cbrt(shape[0]+1))
        frac_along = small_grid/big_grid[-1,1]
        grid_index = ((grid_len-1)*frac_along + 0.00001).astype(int)

        # print("assigning points...")
        # #figure out something that doesn't need a for loop, this is really slow
        # print(small_grid[921,:])
        # for i in range(small_grid.shape[0]):
        #     #print(i)
        #     idx.append(np.where(small_grid[i,:]==big_grid)[0][0])
        idx = (((grid_len**1)*grid_index[:,0]) + ((grid_len**2)*grid_index[:,1]) + ((grid_len**0)*grid_index[:,2])).astype(int)
        idx = (grid_len*grid_index[:,0] + grid_len**2*grid_index[:,1] + grid_index[:,2]).astype(int)
        #print(big_grid[idx[1089],:])
        #print(small_grid[-40:,:])

        return idx

    def setCellDensities(self,filename=None):
        #requires supercell to be set
        if (self.supercell is None):
            print("set the cell values before continuing")
            return

        ensemble_mean,ensemble_std = self.ensemble_predict()
        H = self.getH(ensemble_std)
        taper = self.taper(H)

        fullgrid = self.get_full_grid().T

        density = np.zeros((fullgrid.shape[0]))

        #dims = fullgrid.max(0) +1
        #nonzero_idx = np.where(np.in1d(np.ravel_multi_index(fullgrid.T,dims),np.ravel_multi_index(self.supercell.grid.T,dims)))[0]
        nonzero_idx = self.getnonzero_index(fullgrid,self.supercell.grid)

        if (len(nonzero_idx)==len(ensemble_mean)):
            density[nonzero_idx] = ensemble_mean[range(len(nonzero_idx))]
        else:
            print("messed up the assigning")
        density *= taper
        self.supercell.set_edensity({"xyz":fullgrid,"density":density})

        # print(self.supercell.grid[1089,:])
        # print(fullgrid[nonzero_idx[1089],:])
        # print(ensemble_mean[1089]*taper)
        # print(density[nonzero_idx[1089]])

        if (filename is not None):

            wrapper = wrap_inhouse(self.supercell)
            wrapper.write_cell(mp_grid_spacing=[2,2,2],fname="Hnet_test.cell")
            wrapper.write_unformatted_density(fname=filename)
        return
