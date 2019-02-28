from __future__ import absolute_import
import numpy as np
from den_fmt_io import castep_data
#from sympy.physics.quantum.cg import CG as CG_py
from fortran.f90_descriptor import f90wrap_get_cg_tensor as get_CG_tensor

class fingerprints():
    def __init__(self,nmax=6,lmax=6,r_c=6.0,globalFP=True,descriptor="bispectrum"):
        self.nmax = nmax
        self.lmax = lmax
        self.Rc = r_c
        self.includeGlobal = globalFP
        self.descriptor = descriptor#not that I have any other kind to use
        self.fplength = nmax**2 * lmax**2

    def Load_data(self,datapath):
        casData = castep_data(dataPath=datapath)
        casData.load_castepData()
        self.rawInput = casData.data

    def bispectrum(self):
        if (self.rawInput is None):
            print('No data loaded')
            return
        lmax=1
        n0 = lmax +1
        n1 = n0
        n2 = n0
        n3 = 2*lmax +1
        n4 = n3
        n5 = n3
        cg_tensor = get_CG_tensor(lmax,n0,n1,n2,n3,n4,n5)
        ang_momenta = [0,0,0,0,0,0]
        adjustment = [0,0,0,lmax,lmax,lmax]
        print(cg_tensor[0,1,1,lmax,lmax,lmax])

    def get_FP(self):
        if (self.descriptor == 'bispectrum'):
            return self.bispectrum()
        else:
            print('Only Bispectrum implemented at the moment, sorry. Set fingerprints.descriptor to be bispectrum.')
            return


fp = fingerprints()
fp.Load_data('../CastepCalculations/DenNCells/')
fp.get_FP()
