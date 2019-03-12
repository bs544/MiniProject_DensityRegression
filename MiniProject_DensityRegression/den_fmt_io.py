import numpy as np
from parsers.dft.parser_castep import parse
#parser code courtesy of Andrew Fowler
from parsers.structure_class import supercell
import pickle
import os
import copy

path_ = '../CastepCalculations/DenNCells/'



class castep_data():
    def __init__(self,dataPath,filename="castep_density_data",savedir="./data/",datalist=None):
        self.path = dataPath
        self.data = []
        self.relevantFileTypes = ['castep','den_fmt']
        self.threshold = 1000 # non zero threshold: if a density multiplied by this is less than the maximum density, it's pretty much zero and you don't need to include it
        self.filename = filename
        self.datadir = savedir
        self.datalist = datalist
        #self.set_datalist(datalist)

    # def set_datalist(self,datalist):
    #     if (datalist is None):
    #         files = os.listdir(self.datadir)
    #         self.datalist = [i for i in range(len(files))]
    #     else:
    #         self.datalist = datalist

    def get_densities(self,name,final=True,removeZero=True):
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

        if (not removeZero or (not final and removeZero)):
            if(removeZero):
                print("keeping zero densities, they aren't saved anyway")
            dendict['xyz'] = den['xyz']
            dendict['density'] = den['density']
            dendict['nonzero_indices'] = [i for i in range(len(den['density']))]

        elif(final and removeZero):
            maxden = np.amax(den['density'])
            nonzero_indices = [i for i in range(len(den['density'])) if maxden < self.threshold*den['density'][i]]
            dendict['xyz'] = den['xyz'][nonzero_indices]
            dendict['density'] = den['density'][nonzero_indices]
            dendict['nonzero_indices'] = nonzero_indices


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

    def get_diff(self,init_dict,fin_dict):
        diff_dict = {}
        nonzero_indices = fin_dict['nonzero_indices']
        #check the positions match
        if(np.array_equal(fin_dict['xyz'],init_dict['xyz'][nonzero_indices])):
            diff_dict['xyz'] = fin_dict['xyz']
            diff_dict["fin_density"] = fin_dict["density"]
            diff_dict['density'] = fin_dict['density']-init_dict['density'][nonzero_indices]
        else:
            print('Incompatible inital and final density dictionaries')
        return diff_dict


    def load_castepData(self,load=True,save=True):
        #returns list of dictionaries. Each dictionary contains the data of a single file

        if (load):
            print("loading data from pickle file")
            if (os.path.isdir(self.datadir)):
                datafiles = os.listdir(self.datadir)
                if (datafiles is None):
                    print("Files not present, loading data from castep files")
                else:
                    saved_data=[]
                    for i, file in enumerate(datafiles):
                        if (self.datalist is None or i in self.datalist):
                            f = open('{}{}'.format(self.datadir,file),'rb')
                            saved_dict = pickle.load(f)
                            f.close()
                            keys = ["density","xyz","fin_density","cell","elements","positions"]
                            if all(key in saved_dict.keys() for key in keys):
                                saved_data.append(saved_dict)
                    if (len(saved_data)==len(datafiles)):
                        self.data = saved_data
                        print('good load')
                        return
                    else:
                        print("Invalid saved data, loading from castep files")

        files = os.listdir(self.path)
        if (files is None):
            print('Invalid Path')
            return
        casnames = []
        fin_dennames = []
        init_dennames = []
        completenames = []

        #get names of .castep files and .den_fmt files
        #there should be initial and final density files included
        #names look like Hmol_dist1.00.den_fmt
        #or look like Hmol_dist1.00initial.den_fmt
        for f in files:
            fname,decimal,type = f.split('.')
            if (type == 'castep'):
                fname = '{}.{}'.format(fname,decimal)
                casnames.append(fname)
            elif (type == 'den_fmt'):
                if (len(decimal)>2):
                    fname = '{}.{}'.format(fname,decimal[:-7])#-7 to remove initial
                    init_dennames.append(fname)
                else:
                    fname = '{}.{}'.format(fname,decimal)
                    fin_dennames.append(fname)
            else:
                print('Unsupported file of type {} in file'.format(type))
        #check that each .castep file has a corresponding .den_fmt file
        for name in casnames:
            if (name not in fin_dennames or name not in init_dennames):
                print('Unpaired file: {}.castep'.format(name))
            else:
                completenames.append(name)

        for i,name in enumerate(completenames):
            if (self.datalist is None or i in self.datalist):
                celldict = self.get_atoms('{}{}'.format(name,'.castep'))
                fin_dendict = self.get_densities('{}{}'.format(name,'.den_fmt'),final=True,removeZero=True)
                init_dendict = self.get_densities('{}{}{}'.format(name,'initial','.den_fmt'),final=False,removeZero=True)
                dendict = self.get_diff(init_dendict,fin_dendict)
                celldict.update(dendict)
                if (save):
                    if (not os.path.isdir(self.datadir)):
                        os.mkdir(self.datadir)
                    print("Saving castep data into {}{}{}.pckl".format(self.datadir,self.filename,i))
                    f = open('{}{}{}.pckl'.format(self.datadir,self.filename,i),'wb')
                    pickle.dump(celldict,f)
                    f.close()
                self.data.append(celldict)
        # if(save):
        #     print("Saving castep data in {}{}".format(self.filename,'.pckl'))
        #     f = open('{}{}'.format(self.filename,'.pckl'),'wb')
        #     pickle.dump(self.data,f)
        #     f.close()
        # return






#
# getter = castep_data(dataPath=path_)
# getter.load_castepData()
# getter.get_FPsandTargets()
