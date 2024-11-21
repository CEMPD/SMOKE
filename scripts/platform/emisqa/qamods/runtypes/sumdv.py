from __future__ import print_function
from builtins import range
from sys import exit
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray
from qamods.dateloop.inday import InDay
from qamods.default_paths import *
import os.path

def get_dict(species_list, all_hours = False, grid = '', gsdate = '', case = '', sector = '', inPath = '', spec = '', repDays = '', runDays = '', gridDesc = '', 
	ignoreSpec = False, inln = False, interpolate = False, layer = '', region = '', stacks = False, ptsr = False, inFormat = 'NCF', verbosity = False, zipDict = {}):
        
	'''
	Sums the daily values of NCF files from a start date through the number of run dates.
	'''
	if not grid or not gsdate or not case or not sector or not inPath or not spec or not runDays: 
		exit('ERROR: You must specify grid, gsdate, case, sector, speciation, input path, and rundays.')

	print('Summing files...')
	out_dict = {}  # Create the output dictionary that will be of shape { row: { col: { speciesname: val ... } ... } ... }

	# Create a percent completion counter
#	if verbosity == False:
#		ctr = 0
#		cMax = float(len(species_list)) * runDays
#		sys.stdout.write('\r%.2f%%' %(float(ctr)/cMax))
	
	for species_name in species_list:
		if verbosity: 
			print('Creating daily sum for species: %s' %species_name)

		current_day = InDay(gsdate, repDays, runDays, smkDatesPath)

		for day in range(runDays):
			day_mult = current_day.current_mult()

			# Skip days that have already been represented.
			if day_mult == 0:
				if day != (runDays - 1): 
					current_day.iterday()
				continue

			# Set the input file name prefix for inline versus 2D
			if inln:
				inPrefix = 'inln'
			else:
				inPrefix = 'emis'

			if sector.lower() == 'mrggrid':
				inFileName = os.path.join(inPath, 'emis_mole_all_%s_%s_%s_%s.ncf' %(current_day, grid, spec, case))
			elif sector.lower() == 'mrggrid_withbeis':
				inFileName = os.path.join(inPath, 'emis_mole_all_%s_%s_withbeis_%s.ncf' %(current_day, grid, case))
			elif sector.lower() == 'mrggrid_nobeis':
				inFileName = os.path.join(inPath, 'emis_mole_all_%s_%s_nobeis_%s.ncf' %(current_day, grid, case))
			else:
				inFileName = os.path.join(inPath, sector, '%s_mole_%s_%s_%s_%s_%s.ncf' %(inPrefix, sector, current_day, grid, spec, case))  # v5 directory structure 

			inFile = DataFile(inFileName, verbosity, inFormat, ptsr, zipDict)

			if species_name not in inFile.species_list and ignoreSpec:
#				sys.stdout.flush()
#				sys.stdout.write('\n')
				print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(species_name, inFile))
				break

			inArray = inFile.sum_val(species_name, all_hours, grid, gridDesc, ignoreSpec, inln, interpolate, layer, stacks) * day_mult

			if day == 0: 
				SUM = SpeciesArray(inArray, species_name)
			else:
				SUM.add_array(inArray)

			inFile.close_file()

			if day != (runDays - 1): 
				current_day.iterday()

			# Update completion counter		
#			if not verbosity:
#				ctr += 1	
#				sys.stdout.flush()
#				sys.stdout.write('\r%.2f%%' %((float(ctr)/cMax)*100))

		if species_name in inFile.species_list:
			out_dict[species_name] = SUM

#	sys.stdout.write('\n')	
	return out_dict

