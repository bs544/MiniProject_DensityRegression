import tensorflow as tf
from network import NetworkHandler
from fingerprints import fingerprints

fp = fingerprints(lmax = 4,nmax = 5,r_c = 6.0)
N = NetworkHandler(fp)
N.get_data()
N.train(save=True)
