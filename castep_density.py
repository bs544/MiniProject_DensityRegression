from parsers.structure_class import supercell
from fingerprints import fingeprints
from network import NetworkHandler
import numpy as np
import pickle


class Castep_density(NetworkHandler):
    def __init__(self,descriptor=None,network_handler=None,trainNet=False,data_dir='./data/',train_dir='../CastepCalculations/DenNCells/'):
        self.traindir = train_dir
        self.datadir = data_dir
        self.trainNet = trainNet
        self.supercell = None
        self.set_descriptor(descriptor)
        self.set_network_handler(network_handler)

    def set_network_handler(self,network_handler):
        if (isinstance(network_handler,NetworkHandler)):
            self.NetHandler = network_handler
        else:
            self.NetHandler = NetworkHandler(self.descriptor,train_dir=self.traindir)

    def set_descriptor(self,descriptor):

        if (isinstance(descriptor,fingerprints)):
            self.descriptor = descriptor
            self.fplength = self.descriptor.fplength
        else:
            print("descriptor passed was invalid, creating default fingerprints class")
            self.descriptor = fingerprints()
            self.fplength = 20
        return

    def setupNetwork(self):
        #check if Network is set up
        #if not check if the network can be loaded
        #if not try to train one
        if(self.NetHandler.trained or self.NetHandler.loaded):
            return
        else:
            self.NetHandler.load()
            if (not self.NetHandler.loaded and self.trainNet):
                self.NetHandler.get_data()
                self.NetHandler.train()
            elif (not self.NetHandler.loaded and not self.trainNet):
                print("Warning: Network can't be loaded and hasn't been trained results will be bad")
            return


    def get_frac_coords(self,at_posns,cell):
        #coordinates are absolute cartesians, want to put them in as fractional
        #at_posns[atom_idx,cartesian_element]
        frac_posns = np.zeros((at_posns.shape))
        inv_cell = np.linalg.inv(cell)
        for i in range(at_posns.shape[0]):
            frac_posns[i,:] = np.dot(inv_cell,at_posns[i,:])
        return frac_posns

    def taper(self,x,x_cut,scale):
        x_prime = (x_cut - x)/scale
        #zeros = np.zeros((x_prime.shape))
        x_prime4 = (np.where(x_prime>=0,x_prime,0))**4
        taper = x_prime4/(1+x_prime4)
        return taper

    def get_cell_data(self,filename,cell_dir="./Cell_data/",fp_dir="./FP_data/"):
        cell_keys=["cell","at_posns","grid","fin_density"]
        fp_keys=["fingerprints","density"]
        if (os.path.isdir(cell_dir) and os.path.isfile("{}{}".format(cell_dir,filename))):
            f = open("{}{}".format(cell_dir,filename))
            cell_dict = pickle.load(f)
            f.close()
            if ((key in cell_dict.keys() for key in cell_keys).all())
                print("cell data loaded")
            else:
                print("cell data incomplete")


            else:
                print("no cell data")
                return



        elif (os.path.isdir(self.datadir)):


    def set_supercell(self,cell,at_posns,at_species,grid,FP):
        #use Andrew's supercell class as the storage for calculations on each unit cell
        self.supercell = supercell()
        self.supercell.set_cell(cell)
        frac_coords = self.get_frac_coords(at_posns,cell)
        self.supercell.set_positions(frac_coords)
        self.supercell.set_species(at_species)
        if(len(grid) == FP.shape[0]):
            self.supercell.grid = grid
            self.supercell.FP = FP
        else:
            print("incompatible grid and fingerprints")
        return

    def ensemble_predict(self):
        #requires supercell to be set
        if (self.supercell is None):
            print("set the cell values before continuing")
            return
