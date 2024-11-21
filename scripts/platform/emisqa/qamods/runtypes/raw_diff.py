from __future__ import print_function
from sys import exit
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray

def getDict(speciesList, fileName, grid = '', gridDesc = '', ignoreSpec = False, inln = False, interpolate = False, layer = '', region = '', \
        stacks = False, ptsr = False, inFormat = 'NCF', verbosity = False, zipDict = {}):
	'''
	Gets raw difference of two NCF files
	'''
	if len(fileName.split(',')) != 2: 
		exit('ERROR: You must specify two input filenames using the -f argument.')

	file1 = DataFile(fileName.split(',')[0])
	file2 = DataFile(fileName.split(',')[1])

	print('Adding files...')
	outDict = {}  # Create the output dictionary that will be of shape { row: { col: { speciesname: val ... } ... } ... }

	for speciesName in speciesList:

		if verbosity:
			print('Creating domain total for species: %s' %speciesName)

		if speciesName not in file1.speciesList and ignoreSpec:
			print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(speciesName, file1))

		if verbosity: 
			print('Adding for species: %s' %speciesName)
		array1 = SpeciesArray(file1.sumVal(speciesName, allHours, grid, gridDesc, ignoreSpec, inln, interpolate, layer, stacks), speciesName)
		array2 = SpeciesArray(file2.sumVal(speciesName, allHours, grid, gridDesc, ignoreSpec, inln, interpolate, layer, stacks), speciesName)

		DIFF = array1.diffArray(array2)

		outDict[speciesName] = DIFF 

	return outDict

