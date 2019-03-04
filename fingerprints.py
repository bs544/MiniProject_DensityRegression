from __future__ import absolute_import
import numpy as np
from den_fmt_io import castep_data
#from sympy.physics.quantum.cg import CG as CG_py
from fortran.f90_descriptor import f90wrap_get_cg_tensor as get_CG_tensor
from fortran.f90_descriptor import f90wrap_getbispectrum as get_bispect
from fortran.f90_descriptor import f90wrap_bispect_length as bispect_length
from fortran.f90_descriptor import f90wrap_invoverlap as invOverlap
from fortran.f90_descriptor import f90wrap_invsqrtoverlap as invSqrtOverlap
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
    def __init__(self,nmax=6,lmax=6,r_c=6.0,globalFP=True,descriptor="bispectrum"):
        self.nmax = nmax
        self.lmax = lmax
        self.Rc = r_c
        self.includeGlobal = globalFP
        self.descriptor = descriptor#not that I have any other kind to use
        self.fplength = bispect_length(self.nmax,self.lmax)
        self.filename = '{}.npy'.format(self.descriptor)
        #self.readlength = 100#just some default value for debugging, assigned later

    def Load_data(self,datapath):
        casData = castep_data(dataPath=datapath)
        casData.load_castepData()
        self.rawInput = casData.data
        num_files = len(self.rawInput)
        grid_size = len(self.rawInput[0]['xyz'])
        self.readlength = grid_size

        self.density = np.zeros(num_files*self.readlength)
        for i in range(len(self.rawInput)):
            density = self.rawInput[i]['density']
            self.density[i*self.readlength:(i+1)*self.readlength] = density
        return

    def getInitialVariables(self):
        self.W = invSqrtOverlap(self.nmax,self.nmax,self.nmax)
        self.inv_S = invOverlap(self.nmax,self.nmax,self.nmax)
        if (self.descriptor=='bispectrum'):
            nl = self.lmax +1
            nm = 2*self.lmax +1
            self.cg_tensor = get_CG_tensor(self.lmax,nl,nl,nl,nm,nm,nm)
        return

    def bispectrum(self):
        if (self.rawInput is None):
            print('No data loaded')
            return

        self.getInitialVariables()

        grid_size = len(self.rawInput[0]['xyz'])
        num_files = len(self.rawInput)
        self.readlength = grid_size

        self.fingerprint = np.zeros((num_files*self.readlength,self.fplength*2))

        for i in range(num_files):
            cell = (self.rawInput[i]['cell'])
            elements = self.rawInput[i]['elements']
            at_posns = (self.rawInput[i]['positions'])
            grid = (self.rawInput[i]['xyz'])
            density = (self.rawInput[i]['density'])
            natoms = len(elements)

            glob_bispectrum = get_bispect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[0,:]),self.W,self.inv_S,self.cg_tensor,self.fplength,False,self.fplength)
            start = time.time()
            for j in range(self.readlength):
                #print('Iteration: {}'.format(j))
                bispectrum = get_bispect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[j,:]),self.W,self.inv_S,self.cg_tensor,self.fplength,True,self.fplength)
                bispectrum = np.asarray(bispectrum.T,order='C')
                self.fingerprint[(i+1)*j,:self.fplength] = bispectrum
                self.fingerprint[(i+1)*j,self.fplength:] = glob_bispectrum
            end = time.time()
            print(end-start)

    def Save_FP(self):
        #saves the fingerprint array as a numpy file
        if (self.fingerprint is not None):
            np.save(self.filename,self.fingerprint)
        else:
            print("no fingerprint to save")

        return

    def Load_FP(self):
        if (not os.path.isfile(self.filename)):
            print("File not present")
            return False
        else:
            tmp = np.load(self.filename)
            if (tmp.shape[1] == self.fplength*2 and len(tmp.shape) == 2):
                self.fingerprint = tmp
                return True
            else:
                print("file array has wrong shape")
                return False



    def get_FP(self,load=False,save=False):
        if (self.descriptor == 'bispectrum'):
            loaded = False
            if (load):
                loaded = self.Load_FP()
            if (loaded):
                #print(self.fingerprint[3,:])
                return
            else:
                self.bispectrum()
                if (save):
                    self.Save_FP()
                return
        else:
            print('Only Bispectrum implemented at the moment, sorry. Set fingerprints.descriptor = "bispectrum"')
            return

    def plot2D(self):
        if (self.fingerprint is not None):
            #format is [environment,bispectrum element]
            #mean = self.fingerprint.mean(axis=0)
            #self.fingerprint -= mean
            #print(self.fingerprint[3,:])
            grid = (self.rawInput[0]['xyz'])
            plotgrid = grid[:10000,:2]

            plotZ = self.fingerprint[:10000,:135]
            plotZ = normalize(plotZ,axis=0,norm='l2',copy=False)
            Z = np.linalg.norm(plotZ,axis=1)
            print(Z.shape,plotgrid.shape)

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_trisurf(plotgrid[1:,0],plotgrid[1:,1],Z[1:],linewidth = 0.1, antialiased=True)
            plt.show()

        else:
            print("no fingerprints")
        return



# fp = fingerprints(lmax = 4,nmax = 5,r_c = 6.0)
# fp.Load_data('../CastepCalculations/DenNCells/')
# fp.get_FP(load=True,save=True)
# fp.Load_FP()
# fp.plot2D()
