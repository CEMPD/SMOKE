from __future__ import print_function
from builtins import str
from builtins import object
import numpy as np
import gzip, random
from qamods.helpers import data_blocks, conv2jul
import os.path
from qamods.default_paths import tmp_dir

class CSVFile(object):
    """
    Instance of a CSVFile.
    The CSV File should contain gridded time-stepped species data. The "hours" in the file may
      represent any time step, but are typically hours.

    The format requires a header containing the number of columns, rows, hours, and the 
      8 character Gregorian date.
    Those are defined on header lines starting with "#" followed by COLS, ROWS, HOURS, and
      GSDATE respectively.
    The data should be formatted as column number (starting at 1), row number (starting at 1),
      hour, species name, and value. No header of column names are required for the data. If the
      column names are included, then the first column name should start with "col" or "COL".

    ex.
    #COLS 265
    #ROWS 300
    #HOURS 2
    #GSDATE 20110520
    col,row,hour,speciesname,value
    3,5,1,NOX,265.4
    231,35,2,CO,122.3
    ...
   
    """
    def __init__(self, infile_name, verbosity = False, zip_dict = {}):
        self.infile_name = infile_name
        self.verbosity = verbosity
        self.cols = 0
        self.rows = 0
        self.gsdate = ''
        self.hours = 0
        self.units = ''
        self.spec_dict = {}

        self.open_CSV(zip_dict)
        self.sdate = conv2jul(self.gsdate)

    def __str__(self):
        return self.infile_name

    def get_species(self, species_name, grid = '', grid_desc = '', ignore_spec = False, inln = False, stacks = ''):
        '''
        Function to read a species from the object and return the numpy array
        '''
        if species_name not in self.species_list:
            if ignore_spec:
                print('WARNING: The species %s does not exist in the file %s.' %(species_name, self.infile_name))
                spec_shape = self.spec_dict[self.species_list[0]].shape
                species = np.zeros(spec_shape, '>f')
            else:
                raise ValueError('The species %s does not exist in the file %s.' %(species_name, self.infile_name))
        else:
            species = self.spec_dict[species_name]
        data_in = species[:]
        return data_in

    def _check_int(self, x, name):
        try:
            y = int(x)
        except ValueError:
            raise ValueError('"%s" in file header (#%s). Should be integer.' %(x, name))
        else:
            return y

    def _read_data(self, infile):
        '''
        Read in the formatted CSV file
        '''
        with open(infile,'r') as csvfile:
            for ln, line in enumerate(csvfile):
                # Read the meta data headers: #COLS, #ROWS, #HOURS, #GSDATE, #UNITS
                if line.startswith('#'):
                    line = [cell.strip().upper() for cell in line.strip('#').strip().split(' ')]
                    if line[0] in ('COLS','ROWS','HOURS','GSDATE','UNITS'):
                        setattr(self, line[0].lower(), self._check_int(line[1], line[0]))
                else:
                    line = [cell.strip().upper() for cell in line.strip().split(',')]
                    if self.cols == 0 or self.rows == 0 or self.hours == 0 or not self.gsdate:
                        raise ValueError('Must set file header with lines of #COLS, #ROWS, #HOURS, #GSDATE.')
                    elif line[0].startswith('COL'):
                        # Ignore the header
                        pass
                    else:
                        # Check current col, row, hour against the boundaries set in the header
                        col = int(line[0]) - 1
                        if (col + 1) > self.cols or col < 0:
                            print('WARNING: Column %s on line %s outside of column boundary as defined in header.' %(col+1, ln+1))
                        row = int(line[1]) - 1
                        if (row + 1) > self.rows or row < 0:
                            print('WARNING: Row %s on line %s outside of row boundary as defined in header.' %(row+1, ln+1))
                        hour = int(line[2]) - 1
                        if (hour + 1) > self.hours or hour < 0:
                            print('WARNING: Hour %s on line %s outside of hours maximum as defined in header.' %(hour+1, ln+1))
                        species = line[3]
                        self.spec_dict.setdefault(species, np.zeros([self.hours, 1, self.rows, self.cols], '>f'))

                        # Put the values in species dictionary of numpy arrays
                        try:
                            val = float(line[4])
                        except ValueError:
                            raise ValueError('Value for species %s on line %s cannot be converted to a float.' %(spec, ln+1))
                        else:
                            if self.spec_dict[species][hour,0,row,col] == 0:
                                self.spec_dict[species][hour,0,row,col] = val
                            else:
                                print('Duplicate column, row, and hour combination at line %s.' %ln+1)
                                print(line)
                                raise ValueError('Check for additional entry for column, row, hour, and species.')
        self.species_list = list(self.spec_dict.keys())

    def open_CSV(self, zip_dict = {}):
        '''
        Finds the correct CSV file, unzips if necessary
        '''
        if self.verbosity: 
            print('Opening file for reading: %s' %self.infile_name)
        try: 
            file_in = open(self.infile_name, 'r')
        except IOError:
            print('WARNING: %s not available for access.  Attempting to open zipped version.' %self.infile_name)
            if self.verbosity: 
                print('Opening file for reading: %s.gz' %self.infile_name)
            # Try to read a zipped version of the input file
            try: 
                zip_in = gzip.open(self.infile_name+'.gz','rb')
            except IOError:
                raise IOError('%s.gz not available for access.' %self.infile_name)
            else:
                with zip_in:
                    if self.infile_name in zip_dict:
                        # Check if the file has already been unzipped
                        print('Previously unzipped version found.')
                        file_in = zip_dict[self.infile_name]
                    else:
                        tmp_filename = os.path.join(tmp_dir, 'pyqa-tmp-%s.csv' %str(random.randint(100, 9999999)))
                        tmp_file = open(tmp_filename, 'wb')
                        for data in data_blocks(zip_in):
                                tmp_file.write(data)
                        tmp_file.close()
                        zip_dict[self.infile_name] = tmp_filename  # Store the unzipped temporary file name in the zip dictionary
                        file_in = tmp_filename
        else:
            file_in = self.infile_name
        self._read_data(file_in)

