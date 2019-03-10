import tensorflow as tf
from network import NetworkHandler
from fingerprints import fingerprints
from castep_density import Castep_density
import matplotlib.pyplot as plt
import numpy as np

fp = fingerprints(lmax = 4,nmax = 5,r_c =4.0)
N = NetworkHandler(fp,nEpochs=50,train_dir='../CastepCalculations/practice/')
N.get_data()
N.load()
C = Castep_density(fp,N)
C.setupNetwork()
C.get_cell_data('Hcell1.41.pckl')
C.setCellDensities()


#N.load()
#for i in range(N.nNetworks):
#	rmse = N.get_rmse(i,N.X_test,N.y_test)
#	print(rmse)
