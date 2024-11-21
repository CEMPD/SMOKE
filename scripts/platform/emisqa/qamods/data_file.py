from __future__ import division
from builtins import range
from builtins import object
import numpy as np

class DataFile(object):
    """
    Input file wrapper
    """

    def __init__(self, infile, verbosity = False, informat = 'NCF', ptsr = False, zip_dict = {}):
        if informat == 'UAM':
            from qamods.camx.read_uam import CAMxFile
            self.infile = CAMxFile(infile, verbosity, ptsr, zip_dict)
        elif informat == 'CSV':
            from qamods.csv.read_csv import CSVFile
            self.infile = CSVFile(infile, verbosity, zip_dict)
        elif informat == 'NCF':
            from qamods.cmaq.read_ncf import NCFFile
            self.infile = NCFFile(infile, verbosity, zip_dict)
        else:
            raise ValueError('Wrong input format. Specify NCF, UAM, or CSV with -c.')

        self.infile_name = self.infile.infile_name
        self.species_list = self.infile.species_list
        self.sdate = self.infile.sdate
        self.ptsr = ptsr
        self.informat = informat

    def __str__(self):
        return self.infile_name

    def __call__(self):
        return self.infile

    def get_species(self, species_name, grid = '', grid_desc = '', ignore_spec = False, inln = False, ptsr = False, stacks= '', informat = 'NCF'):
        if informat == 'NCF':
            species = self.infile.get_species(species_name, grid, grid_desc, ignore_spec, inln, stacks)
        else:
            species = self.infile.get_species(species_name, ptsr)

        return species

    def dump_val(self, species_name, all_hours = False, grid = '', grid_desc = '', ignore_spec = False, inln = False, interpolate = False, layer = '', stacks = ''):
        '''
        Returns an array of the hourly data of the selected species in the open NCF.
        Optionally takes all_hours as T/F.  If true, all hours in the NCF are dumped.  If false, just the first 24/1 day (ie. 0-23).
        Flattens all layers unless a single layer to use is specified.
        '''
        species = self.get_species(species_name, grid, grid_desc, ignore_spec, inln, self.ptsr, stacks, self.informat)

        if len(species.shape) == 2:
            '''
            Assume that all 2D netCDFs are just ROWxCOL
            '''
            species = species[np.newaxis,np.newaxis,:,:]
        elif len(species.shape) == 3:
            '''
            Assume that all 3D netCDFs are HOURxROWxCOL
            '''
            species = species[:,np.newaxis,:,:]

        if len(species.shape) != 4:
            raise ValueError('Input variable arrays must be of size 2D, 3D, or 4D')

        if all_hours:
            hours = species.shape[0]
        else:
            hours = 24
        data = np.zeros([hours, species.shape[-2], species.shape[-1]], '>f4')

        if layer:
            layer = int(layer) - 1
            if layer not in list(range(species.shape[1])): 
                raise IndexError('The specified layer is out of range.')
            for hour in range(hours):
                data[hour] += species[hour][layer][:][:]
        else:
            for hour in range(hours):
                for layer in range(species.shape[1]):
                    if interpolate:
                        data[hour] += (species[hour][layer][:][:] + species[hour+1][layer][:][:])/2
                    else:
                        data[hour] += species[hour][layer][:][:]
        return data

    def sum_val(self, species_name, all_hours = False, grid = '', grid_desc = '', ignore_spec= False, inln = False, interpolate = False, layer = '', stacks = ''):
        '''
        Returns an array of the summed hourly data of the selected species in the open NCF.
        Optionally takes all_hours as T/F.  If true, all hours in the NCF are summed.  If false, just the first 24/1 day (ie. 0-23).
        Flattens all layers unless a single layer to use is specified.
        '''
        species = self.get_species(species_name, grid, grid_desc, ignore_spec, inln, self.ptsr, stacks, self.informat)

        if all_hours:
            hours = species.shape[0]
        else:
            hours = 24

        data = np.zeros([1, species.shape[2], species.shape[3]], '>f')
        if layer:
            layer = int(layer) - 1
            if layer not in list(range(species.shape[1])):
                raise IndexError('The specified layer is out of range.')
            for hour in range(hours):
                data[0] += species[hour][layer][:][:]
        else:
            for hour in range(hours):
                for layer in range(species.shape[1]):
                    if interpolate == True:
                        data[0] += ((species[hour][layer][:][:] + species[hour+1][layer][:][:])/2)
                    else:
                        data[0] += species[hour][layer][:][:]
        return data

    def close_file(self):
        self.infile.close_file()

