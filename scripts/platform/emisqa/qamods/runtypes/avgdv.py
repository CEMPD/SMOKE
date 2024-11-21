from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from sys import exit
from qamods.data_file import DataFile
import numpy as np
from qamods.species_array import SpeciesArray
from qamods.dateloop.inday import InDay
import os.path

def get_dict(species_list, all_hours = False, grid = '', gsdate = '', case = '', sector = '', in_path = '', spec = '', run_days = '', grid_desc = '', \
	ignore_spec = False, inln = False, interpolate = False, layer = '', region = '', stacks = False, ptsr = False, in_format = 'NCF', verbosity = False, zip_dict = {}):
        
	'''
	Calculates the daily values of NCF files from a start date through the number of run dates.
	'''
	if not grid or not gsdate or not case or not sector or not in_path or not spec or not run_days: 
		exit('ERROR: You must specify grid, gsdate, case, sector, speciation, input path, and rundays.')

	print('Summing files...')
	out_dict = {}  # Create the output dictionary that will be of shape { row: { col: { speciesname: val ... } ... } ... }

	# Create a percent completion counter
#	if verbosity == False:
#		ctr = 0
#		cMax = float(len(species_list)) * run_days
#		sys.stdout.write('\r%.2f%%' %(float(ctr)/cMax))
	
	for species_name in species_list:
		if verbosity: 
			print('Creating daily sum for species: %s' %species_name)

		current_day = InDay(gsdate)

		for day in range(run_days):
			dayMult = current_day.currentMult()

			# Skip days that have already been represented.
			if dayMult == 0:
				if day != (run_days - 1): 
					current_day.iterday()
				continue

			if sector.lower() == 'mrggrid':
				infile_name = os.path.join(in_path, 'emis_mole_all_%s_%s_%s_%s.ncf' %(current_day, grid, spec, case))
			elif sector.lower() == 'mrggrid_withbeis':
				infile_name = os.path.join(in_path, 'emis_mole_all_%s_%s_withbeis_%s.ncf' %(current_day, grid, case))
			elif sector.lower() == 'mrggrid_nobeis':
				infile_name = os.path.join(in_path, 'emis_mole_all_%s_%s_nobeis_%s.ncf' %(current_day, grid, case))
			else:
				infile_name = os.path.join(in_path, sector, '%s_mole_%s_%s_%s_%s.ncf' %(in_prefix, sector, current_day, grid, spec, case))  # v5 directory structure 

			in_file = DataFile(infile_name)

			if species_name not in list(in_file.NCF.variables.keys()) and ignore_spec:
#				sys.stdout.flush()
#				sys.stdout.write('\n')
				print('WARNING: The species %s does not exist in the file %s.  Skipping.' %(species_name, in_file))
				break

			in_array = in_file.sum_val(species_name, layer, all_hours) * day_mult

			if day == 0: 
				SUM = SpeciesArray(in_array, species_name)
			else:
				SUM.add_array(in_array)

			in_file.close_file()

			if day != (run_days - 1): 
				current_day.iterday()

			# Update completion counter		
#			if verbosity == False:
#				ctr += 1	
#				sys.stdout.flush()
#				sys.stdout.write('\r%.2f%%' %((float(ctr)/cMax)*100))

		if species_name in list(in_file.NCF.variables.keys()):
			out_dict[species_name] = old_div(SUM(), run_days)

#	sys.stdout.write('\n')	
	return out_dict

