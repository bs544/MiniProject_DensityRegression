from parsers.structure_class import supercell
from fingerprints import fingeprints
from network import NetworkHandler
import numpy as np


class Castep_density(NetworkHandler):
    def __init__(self,descriptor=None,network_handler=None,data_dir='./data/',train_dir='../CastepCalculations/DenNCells/'):
        self.traindir = train_dir
        self.datadir = data_dir
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

    def check_train(self):
        if (self.NetHandler.trained):
            return
        else:
            self.NetHandler.get_data()
            self.NetHandler.train()
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
        x_prime = np.where(x_prime>=0,x_prime,0)
        x_prime4 = x_prime**4
        taper = x_prime4/(1+x_prime4)
        return taper

    def set_supercell(self,cell,at_posns,at_species,grid):
        #use Andrew's supercell class as the storage for calculations on each unit cell
        self.supercell = supercell()
        self.supercell.set_cell(cell)
        frac_coords = self.get_frac_coords(at_posns,cell)
        print(at_posns)
        print(frac_coords)
        self.supercell.set_positions(frac_coords)
        self.supercell.set_species(at_species)
