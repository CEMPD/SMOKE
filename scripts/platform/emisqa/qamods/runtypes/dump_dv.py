from __future__ import print_function
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray

def get_dict(species_list, fileName, all_hours = False, grid = '', grid_desc = '', ignore_spec = False, inln = False, interpolate = False, layer = '', region = '', 
	stacks = False, ptsr = False, inFormat = 'NCF', verbosity = False, zipDict = {}):
	'''
	Dumps all summed column and rows for the day.
	'''
	if not fileName: 
		raise ValueError('You must specify an input filename using the -f argument.')

	file1 = DataFile(fileName, verbosity, inFormat, ptsr, zipDict)

	print('Writing Domain...')
	out_dict = dict()

	for species_name in species_list:
		if verbosity: 
			print('Creating summed dump for species: %s' %species_name)

		if species_name not in file1.species_list and ignore_spec:
			print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(species_name, file1))
			continue	

		DV = SpeciesArray(file1.sum_val(species_name, all_hours, grid, grid_desc, ignore_spec, inln, interpolate, layer, stacks), species_name)

		out_dict[species_name] = DV

	return out_dict

