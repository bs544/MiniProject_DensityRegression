import tensorflow as tf
from network import NetworkHandler
from fingerprints import fingerprints
import matplotlib.pyplot as plt
import numpy as np

fp = fingerprints(lmax = 4,nmax = 5,r_c =4.0)
N = NetworkHandler(fp,train_dir='../CastepCalculations/practice/')
N.get_data()
N.train(save=True)
N.plot_loss()
N.plot_rmse()
# def taper(x,x_cut,scale):
#     x_prime = (x_cut - x)/scale
#     #zeros = np.zeros((x_prime.shape))
#     x_prime = np.where(x_prime>=0,x_prime,0)
#     x_prime4 = x_prime**4
#     taper = x_prime4/(1+x_prime4)
#     return taper
#
# x = np.linspace(0.0,0.15,1500)
# x_cut = 0.1
# scale = np.linspace(0.001,0.01,10)
# for x_scale in scale:
#     tap = taper(x,x_cut,x_scale)
#     plt.plot(x,tap)
#     print(x_scale)
#     plt.show()
