from __future__ import division
from __future__ import print_function
from builtins import chr
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import numpy as np
import gzip, random
import os.path
from qamods.default_paths import tmp_dir
from qamods.helpers import data_blocks

### FUTURE UPDATES - USE STRUCT GETSIZE TO GET OFFSET SIZES

class CAMxFile(object):

    def __init__(self, infile, verbosity = False, ptsr = False, zip_dict = {}):
        self.verbosity = verbosity
        self.infile_name = infile
        self.camx = self.open_uam(zip_dict)
        # Defining binary file formats
        #  format = [('dataname','[number of words][endian-ness][datatype]'), ...]
        #  ex: [('sep','>i'),('note','60>i')] - the first word in the binary file will be called "sep", it is a big endian 32 bit integers; 
        #	the next sixty words will be called "note", they are big endian 32 bit integers which will be converted into ASCII text
        # File header format definition
        self.h_dtype = [('sep','>i'),('name','10>i'),('note','60>i'),('one','>i'),('spec','>i'),('sdate','>i'),('stime','>f'),('edate','>i'),('etime','>f')] # 77 * 4 = 308
        # Grid information format definition
        self.g_dtype = [('sep','2>i'),('x-utm','>f'), ('y-utm','>f'), ('zone-utm','>i'), ('xorig','>f'), ('yorig','>f'), ('xcell','>f'), ('ycell','>f'),\
            ('cols','>i'), ('rows','>i'), ('z-cells','>i'), ('cell-break','>i'), ('cell-top','>i'), ('height-surf','>f'), ('h-break','>f'), ('h-cell','>f')]  # 17 * 4 = 68
        # Total number of words in file header (77) + grid information (17) + segment information (8)
        #  This will help determine the total header offset 
        self.camx_headlen = 77 + 17 + 8 
        self._read_head()
        if ptsr:
            self._read_stackdata() 

    def open_uam(self, zip_dict = {}):
        '''
        Opens UAM file for reading.
        '''
        if self.verbosity: 
            print('Opening file for reading: %s' %self.infile_name)
        try: 
            file_in = open(self.infile_name, 'rb')
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
                        file_in = open(zip_dict[self.infile_name], 'rb')
                    else:
                        tmp_filename = os.path.join(tmp_dir('pyqa-tmp-%s.uam' %str(random.randint(100, 9999999))))
                        tmp_file = open(tmp_filename, 'wb')
                        for data in data_blocks(zip_in):
                            tmp_file.write(data)
                        tmp_file.close()
                        zip_dict[self.infile_name] = tmp_filename  # Store the unzipped temporary file name in the zip dictionary
                        file_in = open(tmp_filename, 'rb')
        return file_in

    def _read_head(self):
        # Get Name, species number, and date/time information
        head = np.fromfile(self.camx, np.dtype(self.h_dtype), count = 1)
        self.specnum = head[0]['spec']
#		self.tstep = ((head[0]['edate'] - head[0]['sdate']) * 24) + (int(head[0]['etime'] - head[0]['stime']) + 1)
        self.tstep = ((head[0]['edate'] - head[0]['sdate']) * 24) + (int(head[0]['etime'] - head[0]['stime']))
        self.sdate = str(head[0]['sdate'])
        if self.sdate.startswith('9'):
            self.sdate = '19'+ self.sdate
        else:
            self.sdate = '20' + self.sdate
        grid = np.fromfile(self.camx, np.dtype(self.g_dtype), count = 1)
        self.ncols = grid[0]['cols']
        self.nrows = grid[0]['rows']
        self.nlays = 1
        self.cell = grid[0]['xcell']
        self.xorig = grid[0]['xorig']
        self.yorig = grid[0]['yorig']
        # Get species list	
        s_dtype = [('sep','8>i')] # 8 * 4 = 32  # Segment information, can skip
        s_dtype += [('spec%s' %x,'10>i') for x in range(self.specnum)] # specnum * 10 * 4
        species = np.fromfile(self.camx, np.dtype(s_dtype), count = 1)
        self.species_list = [self.int_to_str(species[0][x]).strip() for x in range(1,len(species[0]))]
        self.camx_headlen += len(self.species_list) * 10

    def int_to_str(self, int_list):
        """
        Get characters from stored integers
        Output as a string
        """
        str_out=''
        for i in int_list:
            try:
                str_out += chr(old_div(((old_div(((old_div((i-32),256))-32),256))-32),256))
            except ValueError:
                print('Warning: Could not convert binary data to string: %s' %i)
        return str_out

    def _read_stackdata(self):
        # Get the stack data for converting from stack to grid
        self.stack_list = []
        self.camx.seek(self.camx_headlen*4)
        stack_head = np.fromfile(self.camx, np.dtype([('sep1','2>i'),('seg','>i'),('npmax','>i'),('sep2','2>i')]), count = 1)
        self.camx_headlen += 6
        self.stacks = stack_head[0]['npmax']
        stack_data = np.fromfile(self.camx, np.dtype([('x','>f'),('y','>f'),('col','>i'),('row','>i'),('h','>f'),('dia','>f')]), count = self.stacks)
        self.camx_headlen += 6 * self.stacks
        for stack in range(self.stacks):
            self.stack_list.append({'col': int(abs(old_div((stack_data[stack]['x'] - self.xorig), self.cell))), 'row': int(abs(old_div((stack_data[stack]['y'] - self.yorig), self.cell)))})

    def get_species(self, species_name, ptsr = False):
        if ptsr:
            species = self.get_ptsrspecies(species_name)
        else:
            species = self.get_emisspecies(species_name)
        return species

    def get_emisspecies(self, species_name):
        """
        Get the species array in TSTEP,LAYER,ROW,COL format
        """
        try:
            spec_idx = self.species_list.index(species_name)
        except ValueError:
            raise ValueError('Species %s not available in %s' %(species_name, self.infile_name))
        species = np.zeros([self.tstep, self.nlays, self.nrows, self.ncols], 'f')
        for hour in range(self.tstep):
            # Location of species data = len of file header + ((time steps + 1) * length of times tep header) + ((species number + 1) * length of species header) + (species number * cols * rows * layers) \
            #      + (time step * (total number of species * (length of species header + (cols * rows * layers)))) 
            spec_loc = self.camx_headlen + ((hour + 1) * 6) + ((spec_idx + 1) * 13) + (spec_idx * self.ncols * self.nrows * self.nlays) + ( hour * (self.specnum * (13 + (self.ncols * self.nrows * self.nlays))) )
#			self.camx.seek((spec_loc - 10) * 4)
#			specHead = np.fromfile(self.camx, np.dtype([('spec','10>i')]), count = 1)
#			specName = self.int_to_str(specHead[0]['spec']).strip()
            self.camx.seek(spec_loc * 4)
            data = np.fromfile(self.camx, np.dtype(('>f',(self.nrows,self.ncols))), count = self.nlays) # (ncols*nrows*4) or datasizebytes
            species[hour][:] = data
        return species

    def get_ptsrspecies(self, species_name):
        """
        Get the species array in TSTEP,LAYER,ROW,COL format from a ptsr file
        """
        try:
            spec_idx = self.species_list.index(species_name)
        except ValueError:
            raise ValueError('Species %s not available in %s' %(species_name, self.infile_name))
        self.tstep = 24
        species_in = np.zeros([self.tstep, 1, self.stacks], 'f')
        species = np.zeros([self.tstep, self.nlays, self.nrows, self.ncols], 'f')
        for hour in range(self.tstep):
            # Location of species data = len of file header + ((time steps + 1) * length of times tep header) + ((species number + 1) * length of species header) + (species number * cols * rows * layers) \
            #      + (time step * (total number of species * (length of species header + (cols * rows * layers)))) 
            spec_loc = self.camx_headlen + ((hour + 1) * 12) + ((hour + 1) * 5 * self.stacks) + ((spec_idx + 1) * 13) + (spec_idx * self.stacks) + ( hour * (self.specnum * (13 + (self.stacks))) )

#			self.camx.seek((spec_loc - 10) * 4)
#			specHead = np.fromfile(self.camx, np.dtype([('spec','10>i')]), count = 1)
#			print int_to_str(specHead[0]['spec'])
            self.camx.seek(spec_loc * 4)
            data = np.fromfile(self.camx, np.dtype(('>f',(self.stacks,))), count = 1) # (ncols*nrows*4) or datasizebytes
            species_in[hour][:] = data
        for stack in range(self.stacks):
            col = self.stack_list[stack]['col']
            row = self.stack_list[stack]['row']
            if row not in list(range(species.shape[2])) or col not in list(range(species.shape[3])):
#				print('stack: %s at col: %s row: %s outside of bounds' %(stack + 1, col + 1, row + 1))
                continue
            try:
                species[:,0,row,col] += species_in[:,0,stack]
            except IndexError:
                raise IndexError('Inline to grid problem at: ROW %s COL %s STACK %s' %(row+1,col+1,stack+1))
        return species

