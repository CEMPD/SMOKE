import numpy as np
import gzip, random
from qamods.helpers import data_blocks, conv2jul
import os.path
from qamods.default_paths import tmp_dir

class CSVFile(object):
	"""
	Instance of a CSVFile.
	"""
	def __init__(self, infile_name, verbosity = False, zip_dict = {}):
		self.infile_name = infile_name
		self.verbosity = verbosity
		self.cols = 0
		self.rows = 0
		self.gsdate = ''
		self.hours = 0
		self.units = ''
		self.species_list = []

		self.csv_file = self.open_CSV(zip_dict)

		self.sdate = conv2jul(self.gsdate)

	def __str__(self):
		return self.infile_name

	def get_species(self, species_name, grid = '', grid_desc = '', ignore_spec = False, inln = False, stacks = ''):
		if species_name not in self.species_list:
			if ignore_spec:
				print 'WARNING: The species %s does not exist in the file %s.' %(species_name, self.infile_name)
				spec_shape = self.spec_dict[self.species_list[0]].shape
				species = np.zeros(spec_shape, '>f')
			else:
				raise ValueError, 'The species %s does not exist in the file %s.' %(species_name, self.infile_name)
		else:
			species = self.spec_dict[species_name]

		data_in = species[:]

		return data_in

	def _check_int(self, x, name, long_name):
		try:
			y = int(x)
		except ValueError:
			raise ValueError, '%s defined as "%s" in file header (#%s). Should be integer.' %(long_name, x, name)

	def _read_data(self, inFile):
		ln = 0
		for line in inFile:
			ln += 1

			# Read the meta data headers: #COLS, #ROWS, #HOURS, #GSDATE, #UNITS
			if line.startswith('#'):
				line = [cell.strip() for cell in line.strip('#').strip().split(' ')]

				metaKey = line[0].upper()
				if metaKey == 'COLS':
					self._check_int(line[1], 'COLS', 'Columns')
					self.cols = int(line[1])
				if metaKey == 'ROWS':
					self._check_int(line[1], 'ROWS', 'Rows')
					self.rows = int(line[1])
				elif metaKey == 'GSDATE':
					self._check_int(line[1], 'GSDATE', 'Starting date')
					self.gsdate = line[1]
				elif metaKey == 'UNITS':
					self.units = line[1]
				elif metaKey == 'HOURS':
					self._check_int(line[1], 'HOURS', 'Hours')
					self.hours = int(line[1])

				continue
			
			line = [cell.strip().upper() for cell in line.strip().split(',')]

			# Read column header: col, row, hour, speciesname1, speciesname2, ...
			if line[0] == 'COL':
				if self.cols == 0 or self.rows == 0 or self.hours == 0 or not self.gsdate:
					raise ValueError, 'Must set file header with lines of #COLS, #ROWS, #HOURS, #GSDATE.'

				head_dict = dict((col, x) for x, col in enumerate(line))
				self.species_list = line[3:]
				self.spec_dict = dict((spec, np.zeros([self.hours, 1, self.rows, self.cols], '>f')) for spec in self.species_list)

				continue

			if self.species_list == []:
				raise ValueError, 'No species found in column header. Set column names starting with col, row, and hour: eg. col,row,hour,SPEC1,SPEC2,SPEC3...'

			# Check current col, row, hour against the boundaries set in the header
			col = int(line[0]) - 1
			if (col + 1) > self.cols or col < 0:
				print 'WARNING: Column %s on line %s outside of column boundary as defined in header.' %(col + 1, ln)
			row = int(line[1]) - 1
			if (row + 1) > self.rows or row < 0:
				print 'WARNING: Row %s on line %s outside of row boundary as defined in header.' %(row + 1, ln)
			hour = int(line[2]) - 1
			if (hour + 1) > self.hours or hour < 0:
				print 'WARNING: Hour %s on line %s outside of hours maximum as defined in header.' %(hour + 1, ln)

			# Put the values in the np data dictionary
			for spec in self.species_list:
				try:
					val = float(line[head_dict[spec]])
				except IndexError:
					raise IndexError, 'Column for species %s not found on line %s' %(spec, ln)
				except ValueError:
					raise ValueError, 'Value for species %s on line %s cannot be converted to a float.' %(spec, ln)

				if self.spec_dict[spec][hour,0,row,col] == 0:
					self.spec_dict[spec][hour,0,row,col] = val
				else:
					print 'Duplicate column, row, and hour combination at line %s.' %ln
					print line
					raise ValueError, 'Check for additional entry for column, row, and hour.'

	def open_CSV(self, zip_dict = {}):
		'''
		Opens the CSV input file and returns an open file object.
		'''
		if self.verbosity: 
			print 'Opening file for reading: %s' %self.infile_name
		try: 
			file_in = open(self.infile_name, 'r')
		except IOError:
			print 'WARNING: %s not available for access.  Attempting to open zipped version.' %self.infile_name
			if self.verbosity: 
				print 'Opening file for reading: %s.gz' %self.infile_name
			# Try to read a zipped version of the input file
			try: 
				zip_in = gzip.open(self.infile_name+'.gz','rb')
			except IOError:
				raise IOError, '%s.gz not available for access.' %self.infile_name
			else:
				with zip_in:
					if self.infile_name in zip_dict:
						# Check if the file has already been unzipped
						print 'Previously unzipped version found.'
						file_in = open(zip_dict[self.infile_name], 'r')

					tmp_filename = os.path.join(tmp_dir, 'pyqa-tmp-%s.csv' %str(random.randint(100, 9999999)))
					tmp_file = open(tmp_filename, 'wb')
					for data in data_blocks(zip_in):
						tmp_file.write(data)
					tmp_file.close()
					zip_dict[self.infile_name] = tmp_filename  # Store the unzipped temporary file name in the zip dictionary
					file_in = open(tmp_filename, 'r')

		self._read_data(file_in)
		return file_in

