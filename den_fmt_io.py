import numpy as np
from parsers.dft.parser_castep import parse
#parser code courtesy of Andrew Fowler
from parsers.structure_class import supercell
import os
import copy

path_ = '../CastepCalculations/DenNCells/'



class castep_data():
    def __init__(self,dataPath):
        self.path = dataPath
        self.data = []
        self.relevantFileTypes = ['castep','den_fmt']

    def get_densities(self,name):
        #returns a dictionary containing densities and positions and cell parameters taken from the .den_fmt file given
        dendict = {}
        filename = '{}{}'.format(self.path,name)
        type = 'den_fmt'
        l = len(type)
        if (type not in name[-l:]):
            print('need to pass a .den_fmt file to get_densities()')
            return
        den_parser = parse(filename,type)
        den_parser.run()
        den_cell = den_parser.get_supercells()[0]
        den = den_cell.get_edensity()
        dendict['xyz'] = den['xyz']
        dendict['density'] = den['density']

        return dendict

    def get_atoms(self,name):
        #returns an dictionary of positions, elements and cell parameters, the index of the element corresponds to its index in the position array
        celldict = {}
        type = 'castep'
        l = len(type)
        filename = '{}{}'.format(self.path,name)

        if (type not in name[-l:]):
            print('need to pass a .castep file to get_positions()')
            return
        pos_parser = parse(filename,type)
        pos_parser.run()
        pos_cell = pos_parser.get_supercells()[0]
        posns = pos_cell.get_positions()
        species = pos_cell.get_species()
        cell = pos_cell.get_cell()

        #positions are relative coordinates when they come out
        #multiply by cell values to get coordinates in Angstroms
        for i in range(posns.shape[0]):
            posns[i,:] = np.dot(cell,posns[i,:])

        celldict['cell'] = cell
        celldict['elements'] = species
        celldict['positions'] = posns

        return celldict

    def load_castepData(self):
        #returns list of dictionaries. Each dictionary contains the data of a single file
        files = os.listdir(self.path)
        if (files is None):
            print('Invalid Path')
            return
        casnames = []
        dennames = []
        pairednames = []

        #get names of .castep files and .den_fmt files
        for f in files:
            fname,type = f.split('.')
            if (type == 'castep'):
                casnames.append(fname)
            elif (type == 'den_fmt'):
                dennames.append(fname)
            else:
                print('Unsupported file of type {} in file'.format(type))
        #check that each .castep file has a corresponding .den_fmt file
        for name in casnames:
            if (name not in dennames):
                print('Unpaired file: {}.castep'.format(name))
            else:
                pairednames.append(name)

        for name in pairednames:
            celldict = self.get_atoms('{}{}'.format(name,'.castep'))
            dendict = self.get_densities('{}{}'.format(name,'.den_fmt'))
            celldict.update(dendict)
            self.data.append(celldict)



#
# getter = castep_data(dataPath=path_)
# getter.load_castepData()
# getter.get_FPsandTargets()
