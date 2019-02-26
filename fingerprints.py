from __future__ import absolute_import
import numpy as np
from den_fmt_io import castep_data


class fingerprints():
    def __init__(self,nmax=6,lmax=6,r_c=6.0,globalFP=True,descriptor="bispectrum"):
        self.nmax = nmax
        self.lmax = lmax
        self.Rc = r_c
        self.includeGlobal = globalFP
        self.descriptor = "bispectrum"#not that I have any other kind to use
        self.fplength = nmax**2 * lmax**2

    def Load_data(self,datapath):
        casData = castep_data(dataPath=datapath)
        casData.load_castepData()
        self.rawInput = casData.data

    def bispectrum(self):
        if (self.rawInput is None):
            print('No data loaded')
            return
        print(self.rawInput[0].keys())


fp = fingerprints()
fp.Load_data('../CastepCalculations/DenNCells/')
fp.bispectrum()
