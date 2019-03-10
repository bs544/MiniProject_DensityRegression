import tensorflow as tf
from fingerprints import fingerprints
import numpy as np
import pickle
import os
import time
import matplotlib.pyplot as plt

#hefty chunks copied outright from Andrew Fowler's code:
#    https://github.com/andrew31416/
#as well as Anirudh Vemula's code:
#   https://github.com/vvanirudh/deep-ensembles-uncertainty.git


class NetworkHandler():
    def __init__(self,descriptor,nodes=[50,50],nNetworks=3,means_per_network=1,train_params=None,train_fraction=0.9,batch_size=50,nEpochs=150,activation="sigmoid",name="DensityEnsemble",train_dir=None):
        self.nodes = nodes
        self.set_train_params(train_params)
        self.set_activation(activation)
        self.nNetworks = nNetworks
        self.nMeans_per_network = means_per_network
        self.batch_size = 500
        self.nEpochs = nEpochs
        self.name=name
        self.trained = False
        self.loaded = False
        if (train_dir is not None):
            self.train_dir=train_dir
        else:
            self.train_dir = './train_data/'
        self.set_descriptor(descriptor)#should already have n_max and l_max and all the rest initialised
        self.train_fraction = train_fraction

    def set_train_params(self,train_params):
        #sets parameters for the NN training run
        default_params = {"learning_rate":0.00005,"optimizer":"rmsprop"}#can't think of any more at the moment
        self.train_params = default_params
        if (train_params is not None and isinstance(train_params,dict)):
            for key in train_params.keys():
                self.train_params[key] = train_params[key]
        elif (train_params is not None and not isinstance(train_params,dict)):
            print("please input train_params as a dictionary")
        return

    def set_activation(self,activation):
        #sets the activation function for the NN nodes
        supported_activations = ["relu","tanh","sigmoid"]
        self.activation = "relu"
        if (activation not in supported_activations):
            print("chosen activation not supported")
        else:
            self.activation = activation
        return

    def set_descriptor(self,descriptor):

        if (isinstance(descriptor,fingerprints) or "lmax" in descriptor.__dict__.keys()): #figure out how to get the isinstance thing working when the descriptor is imported through a module
            self.descriptor = descriptor
            self.fplength = self.descriptor.fplength
        else:
            print("descriptor passed was invalid, creating default fingerprints class")
            self.descriptor = fingerprints()
            self.fplength = 20



        return

    def get_data(self,clear=True):
        #get the fingerprints and densities, randomly shuffle them and assign them to training or testing according to the train fraction

        #get the densities from the castep output
        self.descriptor.Load_data(self.train_dir)

        #get fingerprints
        self.descriptor.get_FP(save=True)

        #center data on zero and standardise
        self.descriptor.standardise_FP()
        self.y_mean,self.y_std = self.descriptor.analyse_densities()

        #randomly shuffle densities and fingerprints in the same way
        indices = np.random.choice(range(len(self.descriptor.density)),len(self.descriptor.density),replace=False)

        train_choice = np.random.choice([1,0],size=(len(indices)),p=[self.train_fraction,1-self.train_fraction])
        #train_indices = np.extract(train_choice==1,np.linspace(0,len(train_choice)-1,len(train_choice)))

        train_data = np.sum(train_choice)
        test_data = len(train_choice) - train_data
        self.X_train = np.zeros((train_data,self.descriptor.fplength*2))
        self.y_train = np.zeros(train_data)

        self.X_test = np.zeros((test_data,self.descriptor.fplength*2))
        self.y_test = np.zeros(test_data)

        train_counter = 0
        test_counter = 0
        for idx in indices:
            if (train_choice[idx]==1):
                self.X_train[train_counter,:] = self.descriptor.fingerprint[idx,:]
                self.y_train[train_counter] = self.descriptor.density[idx]
                train_counter +=1
            else:
                self.X_test[test_counter,:] = self.descriptor.fingerprint[idx,:]
                self.y_test[test_counter] = self.descriptor.density[idx]
                test_counter +=1

        if (clear):
            self.descriptor.fingerprint = None
            self.descriptor.densities = None

        return

    def shuffle_data_pairs(self,X,y):
        #if data needs to be shuffled while keeping link between X and y
        if (isinstance(X,np.ndarray) and isinstance(y,np.ndarray) and X.shape[0] == len(y)):
            indices = np.random.choice(range(len(y)),len(y),replace=False)
            X = X[indices,:]
            y = y[indices]
        else:
            print("X and y are invalid inputs")

        return X, y

    def get_batches(self,X,y,nBatches=None):
        if (X is not None and y is not None):
            X,y = self.shuffle_data_pairs(X,y)
            num_train_pairs = len(y)
            if (nBatches is None):
                nBatches = int(np.floor(num_train_pairs/self.batch_size))
                batch_size = self.batch_size
            else:
                batch_size = int(np.floor(num_train_pairs/nBatches))
            X_batches = []
            y_batches = []
            for i in range(nBatches):
                X_batches.append(X[i*batch_size:(i+1)*batch_size,:])
                y_batches.append(y[i*batch_size:(i+1)*batch_size])


            return X_batches, y_batches, nBatches
        else:
            print("need training data to get batches")
            return


    def setup_ensemble(self):

        if (self.nMeans_per_network == 1):
            out_nodes=[2]
        elif (self.nMeans_per_network > 1):
            out_nodes = [3*self.nMeans_per_network]
        else:
            print("need at least one mean per network, setting value to 1")
            self.nMeans_per_network = 1
            out_nodes=[2]

        self.session["ensemble"] = [Network(self.nodes,out_nodes,self.fplength,self.train_params,index=i,activation=self.activation) for i in range(self.nNetworks)]

        return

    def train(self,save=False):
        start = time.time()
        self.session = {"tf_session":None,"ensemble":None,"saver":None}
        self.setup_ensemble()

        self.session["saver"] = tf.train.Saver([_v for _v in tf.global_variables() if "RMSProp" not in _v.name])

        self.session["tf_session"] = tf.Session()
        self.session["tf_session"].run(tf.global_variables_initializer())

        self.loss = [[] for i in range(self.nNetworks)]
        self.train_rmse = [[] for i in range(self.nNetworks)]
        self.test_rmse = [[] for i in range(self.nNetworks)]

        X_batches, y_batches, nbatches = self.get_batches(self.X_train,self.y_train,nBatches=self.nNetworks)


        for net_index, network in enumerate(self.session["ensemble"]):

            if (True):
                self.session["tf_session"].run(tf.assign(network.output_mean,self.y_mean))
                self.session["tf_session"].run(tf.assign(network.output_std,self.y_std))

            print('Time: {}, '.format(time.time()-start),'network {}:'.format(net_index+1))

            net_X_batch = X_batches[net_index]
            net_y_batch = y_batches[net_index]


            counter = 0
            for epoch in range(self.nEpochs):
                print('epoch {}'.format(epoch+1))
                net_X_minibatches, net_y_minibatches, nbatches = self.get_batches(net_X_batch,net_y_batch)

                for minibatch_idx in range(nbatches):

                    input = {network.in_vect:net_X_minibatches[minibatch_idx],network.target:net_y_minibatches[minibatch_idx].reshape(-1,1)}

                    trainrun,loss = self.session["tf_session"].run([network.train_op,network.loss_val],input)

                    if (np.mod(counter,10) == 0):
                        self.loss[net_index].append(loss)

                    counter += 1
                train_rmse = self.get_rmse(net_index,net_X_batch,net_y_batch)
                test_rmse = self.get_rmse(net_index,self.X_test,self.y_test)
                self.train_rmse[net_index].append(train_rmse)
                self.test_rmse[net_index].append(test_rmse)
                print('test rmse: {}'.format(train_rmse))
                print('loss: {}'.format(loss))
            self.loss[net_index] = np.asarray(self.loss[net_index])
            self.train_rmse[net_index] = np.asarray(self.train_rmse[net_index])
            self.test_rmse[net_index] = np.asarray(self.test_rmse[net_index])
        print('Training finished, total time: {}'.format(time.time()-start))
        self.trained = True

        if(save):
            self.save()

        return

    def save(self):

        if not os.path.isdir(self.name):
            os.mkdir('./{}'.format(self.name))
        network_dict = {}
        unwanted_keys = ["session","X_batches","X_train","X_test","y_batches","y_train","y_test"]
        for key in self.__dict__:
            if (key not in unwanted_keys):
                network_dict.update({key:getattr(self,key)})
        f = open('./{}/{}.pckl'.format(self.name,'network_dict'),'wb')
        pickle.dump(network_dict,f)

        self.session["saver"].save(self.session["tf_session"],'./{}/{}'.format(self.name,self.name))

        return

    def load(self):
        self.loaded = False

        if (not os.path.isdir(self.name)):
            print("Ensemble Directory Not present")
            return

        self.session = {"tf_session":None,"ensemble":None,"saver":None}

        self.setup_ensemble()
        self.session["tf_session"] = tf.Session()

        self.session["saver"] = tf.train.Saver([_v for _v in tf.global_variables() if "RMSProp" not in _v.name])

        self.session["tf_session"] = tf.Session()
        self.session["tf_session"].run(tf.global_variables_initializer())

        self.session["saver"].restore(self.session["tf_session"],'./{}/{}'.format(self.name,self.name))

        print("Network Ensemble Loaded")

        self.loaded = True
        return

    def predict(self,X,nNetworks=None,standard=True):
        #standard is set to false if the incoming data wasn't centered on zero and its standard deviation set to 1 with the initial training data
        if (nNetworks is None):
            nNetworks = self.nNetworks

        if (not standard):
            X -= self.descriptor.mean
            X = X/(self.descriptor.standev[None,:]+1e-8)

        net_idx = range(nNetworks)

        mean = np.zeros((X.shape[0],nNetworks))
        std = np.zeros((X.shape[0],nNetworks))


        for idx,network in enumerate(self.session["ensemble"]):
            if (idx in net_idx):
                input = {network.in_vect:X}
                mean_,std_ = self.session["tf_session"].run([network.mean,network.std],input)
                mean[:,idx] = mean_[:,0]
                std[:,idx] = std_[:,0]

        ensemble_mean = mean.mean(axis=1)
        ensemble_std = np.sqrt((np.sum(np.square(mean)-np.square(ensemble_mean[:,None]),axis=1)+np.sum(np.square(std),axis=1))/self.nNetworks)


        return ensemble_mean,ensemble_std

    def get_rmse(self,network_index,X,y):
        mean = np.zeros((X.shape[0],1))
        std = np.zeros((X.shape[0],1))

        network = self.session["ensemble"][network_index]
        input = {network.in_vect:X}
        mean, std = self.session["tf_session"].run([network.mean,network.std],input)

        error = mean-y.reshape(-1,1)
        mse = np.mean(np.square(error))
        rmse = np.sqrt(mse)
        return rmse

    def plot_loss(self):
        x = np.linspace(1/len(self.loss[0]),self.nEpochs,len(self.loss[0]))
        for idx in range(self.nNetworks):
            plt.plot(x,self.loss[idx])
            plt.xlabel("Number of Epochs")
            plt.ylabel("Loss Value")
        plt.savefig('lossplot.pdf')
        plt.close()
        return

    def plot_rmse(self):
        x = np.linspace(1,self.nEpochs,len(self.train_rmse[0]))
        for idx in range(self.nNetworks):
            plt.plot(x,self.train_rmse[idx],'r')
            plt.plot(x,self.test_rmse[idx],'b')
            plt.xlabel("Number of Epochs")
            plt.ylabel("RMSE")
            plt.legend(["train RMSE","test RMSE"])
        plt.savefig('rmseplot.pdf')
        plt.close()
        return



class Network():
    def __init__(self,hidden_nodes=[100,100],out_nodes=[2],fplength=200,train_params=None,index=0,activation="relu",name="densityNetwork"):
        self.nodes = [fplength*2]+list(hidden_nodes)+out_nodes #input,hidden,output. The output is mean and uncertainty
        print("network {} created".format(index))
        self.netNumber = index
        if (out_nodes[0] == 2):
            self.numNetworks=1#number of networks treated in this class
        elif(out_nodes[0]%3 == 0):
            self.numNetworks = out_nodes[0]/3
        else:
            print("invalid number of output nodes, either have 2 if a single network is being used by the Network class or 3*n if n networks are being trained in a single instance")
        if (index is not None):
            self.name = '{}{}'.format(name,index)
        else:
            self.name = name
        self.dtype = tf.float64
        self.train_params = train_params
        if (True):
            with tf.variable_scope('{}{}'.format(self.name,"targets")):
                self.output_mean = tf.Variable(0.0,trainable=False,dtype=self.dtype)
                self.output_std = tf.Variable(0.1,trainable=False,dtype=self.dtype)
        self.set_activation(activation)
        self.setup()

    def set_activation(self,activation):
        #sets the activation function for the NN nodes
        supported_activations = ["relu","softmax","sigmoid"]
        if (activation not in supported_activations):
            print("chosen activation not supported, setting to relu")
            self.activation = tf.nn.relu
        else:
            if (activation == 'relu'):
                self.activation = tf.nn.relu
            elif(activation == 'softplus'):
                self.activation = tf.nn.softplus
            elif(activation == "sigmoid"):
                self.activation = tf.nn.sigmoid
        return

    def setup(self):

        def lossFn(target,mean,std,mixweight,numNetworks):
            #loss function for heteroskedastic regression.
            #function is of the form \sum_{n}ln(\sum_{k} \pi_{nk} * \frac{e^{-(mu_{nk}-t_{n})^2/2*\sigma_{nk}}}{\sqrt{\sigma_{nk}}})batch_idx

            #get term in exponential
            mu_ = mean - tf.reshape(tf.tile(target,[1,numNetworks]),tf.shape(mean))

            exp = tf.exp(-0.5*tf.divide(tf.square(mu_),std))

            #get term in sum
            lnloss = tf.reduce_sum(tf.multiply(mixweight,tf.divide(exp,tf.sqrt(std))),axis=1)
            #add small value to reduce noise
            lnloss += tf.fill(dims=tf.shape(lnloss),value=tf.cast(1e-8,dtype=lnloss.dtype))
            #get total
            loss = -tf.reduce_sum(tf.log(lnloss))

            return loss

        if (self.numNetworks > 1):
            print("Haven't got this part working yet")
            return
        else:
            #placeholders for NN and loss functions respectively
            self.in_vect = tf.placeholder(self.dtype,[None,self.nodes[0]])
            self.target = tf.placeholder(self.dtype,[None,1])

            #set up weights and biases as a list of matrices and vectors
            self.weights=[]
            self.biases=[]

            with tf.variable_scope('{}{}'.format(self.name,"NNparams")):
                for i in range(1,len(self.nodes)):
                    self.weights.append(tf.Variable(tf.random_normal([self.nodes[i-1],self.nodes[i]],stddev=0.1,dtype=self.dtype),name='{}{}'.format('weightMatrix_',str(i-1)),dtype=self.dtype))
                    self.biases.append(tf.Variable(tf.random_normal([self.nodes[i]],stddev=0.1,dtype=self.dtype),name='{}{}'.format('biasVector_',str(i-1)),dtype=self.dtype))

            #set up NN function
            x = self.in_vect
            for i in range(0,len(self.nodes)-2):
                x = self.activation(tf.add(tf.matmul(x,self.weights[i]),self.biases[i]))
            #output layer
            out  = tf.add(tf.matmul(x,self.weights[-1]),self.biases[-1])

            if (self.numNetworks>1):
                out_sizes = tf.convert_to_tensor([self.numNetworks,self.numNetworks,self.numNetworks])
                self.mean, self.std, self.mixweight = tf.split(out,out_sizes,axis=1)

                #process mixweights so that they are positive and normalised
                self.mixweight = tf.add(tf.nn.softplus(self.mixweight), tf.fill(dims=tf.shape(self.mixweight),value=tf.cast(1e-8,dtype=self.dtype)))
                #self.mixweight = tf.divide(self.mixweight,tf.cast(tf.add(tf.sqrt(tf.reduce_sum(self.mixweight))),tf.cast(1e-8,dtype=self.dtype)))
                self.mixweight = tf.cast(tf.math.l2_normalize(self.mixweight),self.dtype)
            elif(self.numNetworks==1):
                out_sizes = tf.convert_to_tensor([1,1])
                self.mean,self.std = tf.split(out,out_sizes,axis=1)
                self.mixweight = tf.cast(1.0,dtype=self.dtype)

            self.std = tf.add(tf.nn.softplus(self.std),tf.fill(dims=tf.shape(self.std),value=tf.cast(1e-8,dtype=self.dtype)))
            if (True):
                self.mean = tf.add(tf.multiply(self.mean,tf.fill(dims=tf.shape(self.mean),value=self.output_std)),tf.fill(dims=tf.shape(self.mean),value=self.output_mean))
                self.std = tf.multiply(self.std,tf.fill(dims=tf.shape(self.mean),value=self.output_std**2))

            self.loss_val = lossFn(self.target,self.mean,self.std,self.mixweight,self.numNetworks)

            tfvariables = tf.trainable_variables()

            self.gradients = tf.gradients(self.loss_val,tfvariables)

            if (self.train_params["optimizer"]=="rmsprop"):
                self.optimizer = tf.train.RMSPropOptimizer(self.train_params["learning_rate"])
            elif (self.train_params["optimizer"] == "adam"):
                self.optimizer = tf.train.AdamOptimizer(learning_rate=self.train_params["learning_rate"])
            else:
                print("optimizer not supported yet, using default rmsprop")
                self.optimizer = tf.train.RMSPropOptmizer(self.train_params["learning_rate"])

            self.train_op = self.optimizer.apply_gradients(zip(self.gradients,tfvariables))

            return
