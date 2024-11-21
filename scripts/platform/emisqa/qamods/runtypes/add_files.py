from __future__ import print_function
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray
def get_dict(species_list, file_name, all_hours = False, grid = '', grid_desc = '', ignore_spec = False, inln = False, interpolate = False, layer = '', region = '', 
	stacks = False, ptsr = False, in_format = 'NCF', verbosity = False, zip_dict = {}):
	'''
	Adds together two NCF files
	'''
	if len(file_name.split(',')) != 2: 
		raise ValueError('You must specify two input filenames using the -f argument.')

	file1 = DataFile(file_name.split(',')[0])
	file2 = DataFile(file_name.split(',')[1])

	print('Adding files...')
	out_dict = {}  # Create the output dictionary that will be of shape { row: { col: { speciesname: val ... } ... } ... }

	for species_name in species_list:

		if verbosity:
			print('Creating domain total for species: %s' %species_name)

		if species_name not in file1.species_list and ignore_spec:
			print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(species_name, file1))

		if verbosity: 
			print('Adding for species: %s' %species_name)
		array1 = SpeciesArray(file1.dump_val(species_name, all_hours, grid, grid_desc, ignore_spec, inln, interpolate, layer, stacks), species_name)
		array2 = SpeciesArray(file2.dump_val(species_name, all_hours, grid, grid_desc, ignore_spec, inln, interpolate, layer, stacks), species_name)

		array1.add_array(array2())

		out_dict[species_name] = array1

	return out_dict

