import tensorflow as tf
from network import NetworkHandler
from fingerprints import fingerprints
import matplotlib.pyplot as plt
import numpy as np

fp = fingerprints(lmax = 4,nmax = 5,r_c =4.0)
N = NetworkHandler(fp,train_dir='../CastepCalculations/practice/')
N.get_data()
#N.train(save=True)
#N.plot_loss()
#N.plot_rmse()

N.load()
for i in range(N.nNetworks):
	rmse = N.get_rmse(i,N.X_test,N.y_test)
	print(rmse)
