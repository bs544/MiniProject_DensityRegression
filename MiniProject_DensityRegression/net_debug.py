from network import NetworkHandler
from fingerprints import fingerprints
import numpy as np
import matplotlib.pyplot as plt
import pickle

train_params={"learning_rate":0.001,"optimizer":"rmsprop"}
fp = fingerprints(lmax = 4,nmax = 5,r_c =4.0)
N = NetworkHandler(fp,nodes=[10,10],nNetworks=1,train_params=train_params,batch_size=5,nEpochs=30,activation="sigmoid",name="debug")

# y = y[:,0]
#
f = open("./FP_data/Hcell2.0.pckl",'rb')
dict = pickle.load(f)
f.close()
fps = dict["fingerprints"]
den = dict["density"]
print(len(den))
f = open("./Cell_data/Hcell2.0.pckl",'rb')
dict = pickle.load(f)
f.close()
grid = dict["grid"]
idx = np.where(grid[1,1]==grid[:,1])[0]
idx = np.where(grid[1,2]==grid[idx,2])[0]
line = grid[idx,0]
fpline = fps[idx,:]
denline = den[idx]
# plt.plot(line,np.linalg.norm(fpline,axis=1),'r.')
# plt.show()
# plt.close()
# plt.plot(line,denline,'b.')
# plt.show()
# plt.close()


fp.fingerprint = fps[:,:134]#line[:,:134]
fp.density = den
fp.standardise_FP()
stdfp = fp.fingerprint
N.datamean = fp.mean
N.datastd = fp.standev
N.fplength = 67#135
fp.fplength = 67#135
N.X_train = fp.fingerprint
N.y_train = fp.density
N.X_train, N.y_train = N.shuffle_data_pairs(N.X_train,N.y_train)
N.X_test = fp.fingerprint[0:-1:10,:]
N.y_test = fp.density[0:-1:10]
N.y_mean = np.mean(fp.density)
N.y_std = 0.10#np.std(fp.density)
N.train()

y_pred,y_std = N.predict(stdfp[idx,:],standard=True)


plt.plot(line,denline)
plt.errorbar(line,y_pred,np.sqrt(y_std))
plt.plot(line,y_pred,'r')
plt.show()
plt.close()


# X = np.linspace(-10.0,10,10000).reshape(-1,1)
# y = np.zeros(X.shape)
# y[5000:,:] += np.power(X[5000:,:],2) + np.random.normal(scale=0.5,size=X[5000:,:].shape)
# y[:5000,:] += np.random.normal(scale=0.5,size=X[:5000,:].shape)
# y = y/1000
# fp.fplength=0.5
# fp.fingerprint = X
# fp.density = y[:,0]
# fp.standardise_FP()
#
# N.datamean = fp.mean
# N.datastd = fp.standev
# N.fplength=0.5
# N.train_fraction=0.95
#
# N.X_train = fp.fingerprint
# N.y_train = fp.density
# # N.X_train, N.y_train = N.shuffle_data_pairs(N.X_train,N.y_train)
# N.X_test = X
# N.y_test = y
# N.y_mean=0.0#np.mean(N.y_train)
# N.y_std=0.001#np.std(N.y_test)
#
# N.train()
#
# y_pred,y_std = N.predict(X,standard=False)
#
# plt.plot(X,y)
# plt.errorbar(X,y_pred,np.sqrt(y_std))
# plt.plot(X,y_pred,'r')
# plt.show()
# plt.close()
