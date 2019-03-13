from __future__ import absolute_import
import numpy as np
from MiniProject_DensityRegression.den_fmt_io import castep_data
#from sympy.physics.quantum.cg import CG as CG_py
from MiniProject_DensityRegression.fortran.f90_descriptor import f90wrap_get_cg_tensor as get_CG_tensor
from MiniProject_DensityRegression.fortran.f90_descriptor import f90wrap_getbispectrum as get_bispect
from MiniProject_DensityRegression.fortran.f90_descriptor import f90wrap_bispect_length as bispect_length
from MiniProject_DensityRegression.fortran.f90_descriptor import f90wrap_getpowerspectrum as get_powerspect
from MiniProject_DensityRegression.fortran.f90_descriptor import f90wrap_invsqrtoverlap as invSqrtOverlap
import time
import os
import pickle
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D




def get_f90_array(array):
    #fortran has opposite ordering for arrays to python
    #fortran : [fast, slow]
    #python : [slow, fast]
    #so you need to transpose and order things properly
    return np.asarray(array.T,dtype=np.float64,order='F')

def get_python_array(array):

    return np.asarray(array.T,order='C')

class fingerprints():
    def __init__(self,nmax=6,lmax=6,r_c=6.0,globalFP=True,descriptor="bispectrum",train_dir='./FP_data/',cell_dir='./Cell_data/'):
        self.nmax = nmax
        self.lmax = lmax
        self.Rc = r_c
        self.includeGlobal = globalFP
        self.descriptor = descriptor#not that I have any other kind to use
        self.bilength = bispect_length(self.nmax,self.lmax)
        self.powerlength = nmax*(lmax+1)
        # self.filename = '{}.pckl'.format(self.descriptor)
        self.train_dir = train_dir
        self.cell_dir = cell_dir
        self.max_r = 0.0
        #self.readlength = 100#just some default value for debugging, assigned later

    def Load_data(self,datapath):
        casData = castep_data(dataPath=datapath,datalist=None)
        casData.load_castepData()
        self.rawInput = casData.data
        num_files = len(self.rawInput)
        self.readlength = [len(self.rawInput[i]['xyz']) for i in range(len(self.rawInput))]
        print('readlength: ',self.readlength)

        self.Cells=[]

        self.density = np.zeros(sum(self.readlength))
        for i in range(len(self.rawInput)):
            density = self.rawInput[i]['density']
            fin_density = self.rawInput[i]['fin_density']
            grid = self.rawInput[i]['xyz']
            cell = self.rawInput[i]['cell']
            posns = self.rawInput[i]['positions']
            elements = self.rawInput[i]['elements']

            at_dist = self.get_at_dist(posns)
            name = "Hcell{}".format(str(at_dist))
            cell_class = Cell_data(self.bilength,self.powerlength,grid,cell,posns,elements,name)
            cell_class.set_density(density)
            cell_class.set_final_density(fin_density)
            self.Cells.append(cell_class)
            max_r = self.get_max_r(cell,grid,posns)
            if (max_r>self.max_r):
                self.max_r = max_r
            self.density[sum(self.readlength[:i]):sum(self.readlength[:i+1])] = density
        return

    def getInitialVariables(self):
        self.W = invSqrtOverlap(self.nmax,self.nmax,self.nmax)
        if (self.descriptor=='bispectrum'):
            nl = self.lmax +1
            nm = 2*self.lmax +1
            self.cg_tensor = get_CG_tensor(self.lmax,nl,nl,nl,nm,nm,nm)
        return

    def get_at_dist(self,posns):
        if(posns.shape[1]==3):
            if (posns.shape[0]==2):
                dist = np.linalg.norm(posns[0,:]-posns[1,:])
                dist = round(dist,2)
            else:
                dist=0.0
        else:
            print("Bad atom position shape")
            dist = 0.0
        return dist

    def bispectrum(self):
        if (self.rawInput is None and self.Cells is None):
            print('No data loaded')
            return

        self.getInitialVariables()


        num_files = len(self.Cells)
        print("max_r: ",self.max_r)

        self.fingerprint = np.zeros((sum(self.readlength[:num_files]),self.bilength+self.powerlength))

        # for i in range(num_files):
        #     cell = (self.rawInput[i]['cell'])
        #     elements = self.rawInput[i]['elements']
        #     at_posns = (self.rawInput[i]['positions'])
        #     grid = (self.rawInput[i]['xyz'])
        #     density = (self.rawInput[i]['density'])
        start = time.time()
        for i, cell_class in enumerate(self.Cells):
            cell = cell_class.cell
            elements = cell_class.at_species
            at_posns = cell_class.at_posns
            grid = cell_class.grid
            density = cell_class.density
            natoms = len(elements)

            if (natoms not in at_posns.shape):
                print("invalid atom positions")
            # elif (natoms == at_posns.shape[0]):
            #     at_posns = at_posns.T

            glob_powerspectrum = get_powerspect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[0,:]),self.W,self.powerlength,False,self.powerlength)
            glob_powerspectrum = np.asarray(glob_powerspectrum.T,order="C")
            for j in range(self.readlength[i]):
                #print('Iteration: {}'.format(j))
                bispectrum = get_bispect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[j,:]),self.W,self.cg_tensor,self.bilength,True,self.bilength)
                bispectrum = np.asarray(bispectrum.T,order='C')
                self.fingerprint[sum(self.readlength[:i])+j,:self.bilength] = bispectrum
                self.fingerprint[sum(self.readlength[:i])+j,self.bilength:] = glob_powerspectrum
            cell_class.set_fingerprints(self.fingerprint[sum(self.readlength[:i]):sum(self.readlength[:i+1]),:])
            cell_class.Save_train_data(self.train_dir)
            print('Calculation of file {} complete, time taken:{}'.format(i,time.time()-start))
        end = time.time()
        print('{}{}'.format('Calculation complete, time taken: ',end-start))

    def Save_FP(self):
        #saves the cell class in two files: one contains the fingerprints and different densities, the other contains the rest of the cell information for testing
        for cell_class in self.Cells:
            cell_class.Save_train_data(self.train_dir)
            cell_class.Save_cell_data(self.cell_dir)


        # if (self.fingerprint is not None and self.density is not None):
        #     print("Saving Fingerprints")
        #     save_dict = {"fingeprints":self.fingerprint,"densities":self.density}
        #     f = open(self.filename,'wb')
        #     pickle.dump(save_dict,f)
        #     f.close()
        # else:
        #     print("no fingerprint to save")

        return

    def Load_FP(self):
        #check files in self.train_dir
        if(os.path.isdir(self.train_dir)):
            files = os.listdir(self.train_dir)
            self.fingerprint = np.zeros((1,self.bilength+self.powerlength))
            counter = 0
            self.readlength = []
            for i, file in enumerate(files):
                f = open('{}{}'.format(self.train_dir,file),'rb')
                savedict = pickle.load(f)
                f.close()
                tmp_fp = savedict['fingerprints']
                self.readlength.append(tmp_fp.shape[0])
                if (tmp_fp.shape[1]==(self.bilength+self.powerlength) and len(tmp_fp.shape)==2):
                    self.fingerprint=np.vstack((self.fingerprint,tmp_fp))
                    counter +=1
                else:
                    print("fingerprint array has the wrong shape")
            self.fingerprint = self.fingerprint[1:,:]
            self.density = np.zeros(self.fingerprint.shape[0])
            for i, file in enumerate(files):
                f = open('{}{}'.format(self.train_dir,file),'rb')
                savedict = pickle.load(f)
                f.close()
                tmp_den = savedict['density']
                self.density[sum(self.readlength[:i]):sum(self.readlength[:i+1])]=tmp_den
                counter +=1

            if(counter==2*len(files)):
                return True
            else:
                return False


        else:
            print("File directory not present")
            return False

        # if (not os.path.isfile(self.filename)):
        #     print("File not present")
        #     return False
        # else:
        #     print("File present, loading")
        #     f = open(self.filename,'rb')
        #     save_dict = pickle.load(f)
        #     f.close()
        #     tmp = save_dict["fingeprints"]
        #     if (tmp.shape[1] == self.fplength*2 and len(tmp.shape) == 2):
        #         self.fingerprint = tmp
        #         self.density = save_dict["densities"]
        #         return True
        #     else:
        #         print("file array has wrong shape")
        #         return False

    def get_FP(self,load=True,save=True):
        if (self.descriptor == 'bispectrum'):
            loaded = False
            if (load):
                loaded = self.Load_FP()
            if (loaded):
                print('Loaded fingerprints')
                return
            else:
                print("Didn't load fingerprints, calculating:")
                self.bispectrum()
                if (save):
                    self.Save_FP()
                return
        else:
            print('Only Bispectrum implemented at the moment, sorry. Set fingerprints.descriptor = "bispectrum"')
            return

    def standardise_FP(self):
        #centers each of the bispectrum elements about zero and sets them so their standard deviation is 1
        #stores mean and standard deviation for use in any other data.
        if (self.fingerprint is not None):
            #format is [environment,bispectrum_element]
            self.mean = self.fingerprint.mean(axis=0)
            self.fingerprint -= self.mean
            self.standev = np.std(self.fingerprint,axis=0)
            self.standev.reshape((len(self.standev),1))
            if (all(self.standev[self.bilength:]<1e-15)):
                print("there should only be one cell training")
                self.fingerprint[:,:self.bilength] = self.fingerprint[:,:self.bilength]/(self.standev[None,:self.bilength]+1e-8)
            else:
                self.fingerprint = self.fingerprint/(self.standev[None,:]+1e-8)

    def get_max_r(self,cell,grid,posns):
        #gets the maximum distance from any atom for which the density is non zero
        #since the atoms are all the same, it's probably safe to just find the closed atom and assume that since net electronegativity is zero electrons belong to the closest atom
        iter = [[-1],[0],[1]]
        periodic_posns = [[posns+i*cell[0,:]+j*cell[1,:]+k*cell[2,:]]for i in iter for k in iter for j in iter]
        dists = np.ones((grid.shape[0],len(iter)**3,posns.shape[0]))
        min_dists = np.zeros(grid.shape[0])

        for i in range(len(periodic_posns)):
            at_posns = periodic_posns[i][0]
            for at_idx in range(posns.shape[0]):
                dists[:,i,at_idx] = np.sqrt(np.sum(np.square(grid-at_posns[None,at_idx,:]),axis=1))
        min_dists = np.min(np.min(dists,axis=2),axis=1)
        max_r = np.max(min_dists)
        return max_r

    def analyse_densities(self):
        if (self.density is not None):
            mean = np.mean(self.density)
            std = np.std(self.density)
            return mean, std
        else:
            print("density not loaded")
            return


    def plot2D(self):
        if (self.fingerprint is not None):
            #format is [environment,bispectrum element]
            #mean = self.fingerprint.mean(axis=0)
            #self.fingerprint -= mean
            #print(self.fingerprint[3,:])
            grid = (self.rawInput[0]['xyz'])
            flat_idx = list(np.where(grid[:,2]==1.5)[0])
            plotgrid = grid[flat_idx,:2]

            plotZ = self.fingerprint[flat_idx,:self.bilength]
            plotZ = normalize(plotZ,axis=0,norm='l2',copy=False)
            Z = np.linalg.norm(plotZ,axis=1)
            print(Z.shape,plotgrid.shape)

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_trisurf(plotgrid[1:,0],plotgrid[1:,1],Z[1:],linewidth = 0.1, antialiased=True)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            plt.show()

        else:
            print("no fingerprints")
        return

class Cell_data():
    def __init__(self,bilength,powerlength,grid,cell_vects,at_posns,at_species,name):
        self.cell = cell_vects
        self.at_posns = at_posns#want to be natoms by 3 array
        self.name = name
        self.grid = grid
        self.at_species = at_species
        self.bilength = bilength
        self.powerlength = powerlength
        self.fingerprints = None
        self.fin_density = None
        self.density = None#difference between initial and final densities
        self.check_init()

    def check_init(self):
        if(not isinstance(self.cell,np.ndarray) or not self.cell.shape == (3,3)):
            print("Invalid cell shape and/or variable type is not a numpy nd array")
        if(not isinstance(self.at_posns,np.ndarray) or len(self.at_posns.shape)!=2):
            print("Invalid variable type for atom positions")
        elif(self.at_posns.shape[1]!=3 and self.at_posns.shape[0]!=3):
            print("atom positions aren't in 3D")
        elif(self.at_posns.shape[1]!=3 and self.at_posns.shape[0]==3):
            print("Transposing atom positions")
            self.at_posns = self.at_posns.T
        return

    def check_FPnD(self):
        if (len(self.density)!=self.fingerprints.shape[0]):
            print('Incompatible density and fingerprints')
        return

    def set_fingerprints(self,fingerprints):
        if(isinstance(fingerprints,np.ndarray) and (self.bilength+self.powerlength) in fingerprints.shape):
            if(fingerprints.shape[0]==(self.bilength+self.powerlength)):
                self.fingerprints = fingerprints.T
            else:
                self.fingerprints = fingerprints
            if (self.density is not None):
                self.check_FPnD()
        else:
            print("Invalid type or shape for fingerprints")
            print(type(fingerprints))
            print(fingerprints.shape)
        return

    def set_density(self,density):
        self.density = density
        if(self.fingerprints is not None):
            self.check_FPnD()
        return

    def set_final_density(self,fin_density):
        self.fin_density = fin_density
        if(self.fingerprints is not None):
            self.check_FPnD()
        return

    def Save_train_data(self,savedir):
        if (self.fingerprints is not None and self.density is not None):
            if (not os.path.isdir(savedir)):
                os.mkdir(savedir)
            filename = "{}{}.pckl".format(savedir,self.name)
            if (not os.path.isfile(filename)):
                print("saving fingerprints and densities from cell class {} in {}".format(self.name,savedir))
                dict={"fingerprints":self.fingerprints,"density":self.density}
                f = open(filename,'wb')
                pickle.dump(dict,f)
                f.close()
        else:
            print("Don't have fingerprints or densities")
        return

    def Save_cell_data(self,savedir):
        if (not os.path.isdir(savedir)):
            os.mkdir(savedir)
        filename = "{}{}.pckl".format(savedir,self.name)
        if (not os.path.isfile(filename)):
            dict={"cell":self.cell,"at_posns":self.at_posns,"grid":self.grid,"fin_density":self.fin_density}
            f = open(filename,'wb')
            pickle.dump(dict,f)
            f.close()
        return



# fp = fingerprints(lmax = 4,nmax = 5,r_c = 6.0)
# fp.Load_data('../CastepCalculations/DenNCells/')
# fp.get_FP(load=True,save=True)
# fp.Load_FP()
# fp.plot2D()
