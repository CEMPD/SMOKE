from __future__ import print_function
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray
import sys

def get_dict(species_list, fileName, all_hours = False, grid = '', grid_desc = '', ignore_spec = False, inln = False, interpolate = False, layer = '', region = '', 
    stacks = False, ptsr = False, inFormat = 'NCF', verbosity = False, zipDict = {}):
    '''
    Sums up every grid cell for a day for each species.
    '''
    if not fileName: 
        raise ValueError('You must specify an input filename using the -f argument.')
    if region: 
        raise ValueError('This run type does not support grid to fips conversion.  Please remove -e argument from command line.')

    file1 = DataFile(fileName, verbosity, inFormat, ptsr, zipDict)

    print('Running max/min')
    outDict = dict()

    for species_name in species_list:
        if verbosity: 
            print('Running max/min for species: %s' %species_name)

        if species_name not in file1.species_list and ignore_spec:
            print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(species_name, file1))
            continue    

        array1 = SpeciesArray(file1.sum_val(species_name, all_hours, grid, grid_desc, ignore_spec, inln, interpolate, layer, stacks), species_name)
        mmVals = array1.maxMin()
        print('For species %s\n' %species_name)
        print('Min Value: %s   Max Value: %s\n' %(mmVals[0], mmVals[1]))

    sys.exit(0)
    outdict[species_name] = DV
    return outdict


