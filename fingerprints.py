from __future__ import absolute_import
import numpy as np
from den_fmt_io import castep_data
#from sympy.physics.quantum.cg import CG as CG_py
from fortran.f90_descriptor import f90wrap_get_cg_tensor as get_CG_tensor
from fortran.f90_descriptor import f90wrap_getbispectrum as get_bispect
from fortran.f90_descriptor import f90wrap_bispect_length as bispect_length
import time
import os
import pickle



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

    def Load_data(self,datapath):
        casData = castep_data(dataPath=datapath)
        casData.load_castepData()
        self.rawInput = casData.data

    def bispectrum(self):
        if (self.rawInput is None):
            print('No data loaded')
            return

        nl = self.lmax +1
        nm = 2*self.lmax +1
        cg_tensor = get_CG_tensor(self.lmax,nl,nl,nl,nm,nm,nm)

        grid_size = len(self.rawInput[0]['xyz'])
        num_files = len(self.rawInput)

        self.fingerprint = np.zeros((num_files*grid_size,self.fplength*2))

        for i in range(num_files):
            cell = (self.rawInput[i]['cell'])
            elements = self.rawInput[i]['elements']
            at_posns = (self.rawInput[i]['positions'])
            grid = (self.rawInput[i]['xyz'])
            density = (self.rawInput[i]['density'])
            natoms = len(elements)

            glob_bispectrum = get_bispect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[0,:]),cg_tensor,self.fplength,False,self.fplength)
            start = time.time()
            for j in range(1000):
                #print('Iteration: {}'.format(j))
                bispectrum = get_bispect(self.nmax,self.lmax,self.Rc,get_f90_array(at_posns),get_f90_array(cell),natoms,get_f90_array(grid[4,:]),cg_tensor,self.fplength,True,self.fplength)
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
            print(tmp.shape)
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
                return
            else:
                self.bispectrum()
                if (save):
                    self.Save_FP()
                return
        else:
            print('Only Bispectrum implemented at the moment, sorry. Set fingerprints.descriptor = "bispectrum"')
            return


fp = fingerprints(lmax = 6,nmax = 5,r_c = 3.0)
fp.Load_data('../CastepCalculations/DenNCells/')
fp.get_FP(load=False,save=True)
