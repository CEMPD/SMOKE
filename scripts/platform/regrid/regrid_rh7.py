#!/usr/bin/env python
# Intelligently regrids from one domain to another
# Lambert to Lambert only at this time
# 4/12/13 James Beidler <beidler.james@epa.gov>
# Set dtype to >f4 for compatability on new garnet 7 Jul 2013

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from builtins import object
import numpy as np
#import netcdf_file as ncf 
#from scipy.io.netcdf import *
import netCDF4 as ncf
import sys, os
from datetime import datetime, timedelta

def checkEV(evName):
	"""
	Checks if an environment variable is set.  If not, exits.  If it is, returns the variable.
	Takes the name of the environment variable.
	"""
	try: 
		var = os.environ[evName]
	except:
		print("ERROR: Environment variable '%s' is not defined." %evName)
		sys.exit(1)
	else: 
		return var

# Set the default grid description file.  
def_griddesc = checkEV('GRIDDESC')
#def_griddesc = '/work/MOD3APP/ktalgo/ge_dat/gridding/griddesc_lambertonly_17apr2018_v2.txt'

########


class Grid(object):
    """
    Reads the grid description file and loads the grid information for the specified grid named.
    """
    def __init__(self, grid_name, grid_desc):
        self.name = grid_name
        self.grid_desc = grid_desc
        self._load_gd()

    def _parse_float(self, x):
        """
        Returns a floating point with the correct number of trailing zeros based on the .Dx
        """
        if 'D' in x:    
            num = x.strip().split('.')[0]
            zerofloat = x.strip().split('.')[1][1:]
            numlen = len(num) + int(zerofloat)
            while len(num) < numlen:
                num = num + '0'
            return np.float64(num)
        else:
            return np.float64(x)
    
    def _load_gd(self):
        """
        Read in the grid description file and store the grid data as object attributes
        """
        with open(self.grid_desc) as gd:
            state = 1
            for line in gd.readlines():
                if state == 1:
                    if self.name.strip().upper() in line:
                        state = 2
                    continue
                elif state == 2:
                    split_line = line.split(',')
                    self.xorig = self._parse_float(split_line[1])
                    self.yorig = self._parse_float(split_line[2])
                    self.cell = self._parse_float(split_line[3])
                    self.cols = int(split_line[5])
                    self.rows = int(split_line[6])
                    break
            if state == 1: 
                raise ValueError('Grid %s not found in grid description file.' %self.name)

class GridBounds(object):
    '''
    Set the boundaries of the new grid and the old grid based on the new grid
    '''
    def __init__(self, in_ncf, grid):
        self._set_bounds(in_ncf, grid)

    def _set_bounds(self, in_ncf, grid):
        # Calculate the raw difference for the distance between the origination points on the new and old grids
        x_dist = (grid.xorig - in_ncf.XORIG)
        y_dist = (grid.yorig - in_ncf.YORIG)

        # Calculate the starting points for copying on both grids
        # Keep in mind that the grid starts at the SW corner
        if x_dist < 0: # If the old x origination is inside the new grid
            self.in_col_orig = 0 # Start at the first old grid cell 
            self.out_col_orig = int(abs(old_div(x_dist, grid.cell))) # Offset some new grid cells
        else: # If the old x origination is outside of equal to the orig on new grid
            self.in_col_orig = int(abs(old_div(x_dist, in_ncf.XCELL)))
            self.out_col_orig = 0
        if y_dist < 0:
            self.in_row_orig = 0
            self.out_row_orig = int(abs(old_div(y_dist, grid.cell)))
        else:
            self.in_row_orig = int(abs(old_div(y_dist, in_ncf.XCELL)))
            self.out_row_orig = 0

        # Calculate the end points to form boundary fields
        in_x_end = in_ncf.XORIG + (in_ncf.NCOLS * in_ncf.XCELL)
        in_y_end = in_ncf.YORIG + (in_ncf.NROWS * in_ncf.XCELL)
        out_x_end = grid.xorig + (grid.cols * grid.cell)
        out_y_end = grid.yorig + (grid.rows * grid.cell)
        x_dist = out_x_end - in_x_end
        y_dist = out_y_end - in_y_end

        if x_dist > 0: # If the endpoint of the old is inside the new
            self.in_col_end = in_ncf.NCOLS
            self.out_col_end = int(grid.cols - abs(old_div(x_dist, grid.cell))) # Calculate the last grid cell to aggregate
            if self.out_col_end > int(self.out_col_end): # Check to see if there is anything after the decimal ie. a partial column calculated.
                self.out_col_end = int(self.out_col_end) + 1  # If there is a partial column then add a column to be processed to hold the extra input data
        else:
            self.in_col_end = int(in_ncf.NCOLS - abs(old_div(x_dist, in_ncf.XCELL)))
            self.out_col_end = grid.cols
        if y_dist > 0:
            self.in_row_end = in_ncf.NROWS
            self.out_row_end = int(grid.rows - abs(old_div(y_dist, grid.cell)))

            if self.out_row_end > int(self.out_row_end):
                self.out_row_end = int(self.out_row_end) + 1
        else:
            self.in_row_end = int(in_ncf.NROWS - abs(old_div(y_dist, in_ncf.XCELL)))
            self.out_row_end = grid.rows

    #        print "In start col: %s  End: %s    Out start: %s  End: %s" %(self.in_col_orig, self.in_col_end, self.out_col_orig, self.out_col_end) 
    #        print "In start row: %s  End: %s    Out start: %s  End: %s" %(self.in_row_orig, self.in_row_end, self.out_row_orig, self.out_row_end)
        self.chunk_dim = old_div(float(grid.cell), float(in_ncf.XCELL))  # How many times larger is the output cell versus the input

class NCF(ncf.Dataset):
    """
    NCF subclass
    """
    def __init__(self, file_name, mode='r'):
        print('Opening %s' %file_name)
        ncf.Dataset.__init__(self, file_name, mode, format='NETCDF3_64BIT')

    def _set_dims(self, in_ncf, grid):
        '''
        Set up the outfile dimensions and attributes based on the new grid
        and the input file
        '''
        # Set outfile dimensions 
        out_dims = { 'VAR': len(in_ncf.dimensions['VAR']), 'TSTEP': len(in_ncf.dimensions['TSTEP']), 
            'DATE-TIME': 2, 'ROW': grid.rows, 'COL': grid.cols, 'LAY': len(in_ncf.dimensions['LAY']) }
        for dim, value in list(out_dims.items()): 
            self.createDimension(dim, int(value))

        # Ignore automatically created attributes, plus the grid specific ones and history because that doesn't play well for some reason.
        for att_name in in_ncf.ncattrs():
            if att_name == att_name.upper(): 
                # Checks to see if the attribute is in caps.  All the ones that we create should be.
                att_val = getattr(in_ncf, att_name)
                setattr(self, att_name, att_val) 
        self.NROWS = grid.rows
        self.NCOLS = grid.cols
        self.XORIG = grid.xorig
        self.YORIG = grid.yorig
        self.GDNAM = grid.name
        self.XCELL = grid.cell
        self.YCELL = grid.cell
        self.HISTORY = ' '

    def regrid(self, in_ncf, grid_info):
        """
        the main regridding section
        sets loop over variables and decides how to regrid
        """ 
        grid = grid_info
        self._set_dims(in_ncf, grid)
        var_list = getattr(in_ncf, 'VAR-LIST')
        species_list = [var_list[x*16:(x+1)*16].strip() for x in range(int(in_ncf.NVARS))] 
        bounds = GridBounds(in_ncf, grid)
        print('Regridding %s to %s %s' %(in_ncf.filepath(), self.filepath(), grid.name))
        cpl = 0
        c_max = float(len(species_list))
        sys.stdout.write("\r%.2f%%" %(old_div(float(cpl),c_max)))

        self._write_tflag(in_ncf)
        # Loop through and subset each species variable
        for species_name in species_list:
            species_in = in_ncf.variables[species_name]
            data_in = species_in[:]
            species_out = self.createVariable(species_name, np.float32, ('TSTEP', 'LAY', 'ROW', 'COL'))
            species_out.long_name = species_in.long_name
            species_out.units = species_in.units
            species_out.var_desc = species_in.var_desc
            if in_ncf.XCELL < grid.cell:
                # Aggregate
                data_out = self._agg(data_in, grid, bounds)
            elif in_ncf.XCELL > grid.cell:
                # Subdivide
                data_out = self._sub_div(data_in, grid, bounds)
            else:
                # Window or expand
                data_out = self._subset(data_in, grid, bounds)
            species_out[:] = data_out

            cpl += 1
            sys.stdout.flush()
            sys.stdout.write("\r%.2f%%" %((old_div(float(cpl),c_max))*100.))
        sys.stdout.write("\n")

    def _subset(self, data_in, grid, bounds):
        """
        Window to a grid or expand to a grid of equal cell size
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.rows, grid.cols], np.float32)
        for col_cnt,col in enumerate(range(bounds.out_col_orig, bounds.out_col_end)):
            for row_cnt,row in enumerate(range(bounds.out_row_orig, bounds.out_row_end)):
                i_row = bounds.in_row_orig + row_cnt
                i_col = bounds.in_col_orig + col_cnt
                if i_row not in list(range(bounds.in_row_orig, bounds.in_row_end)) or i_col not in list(range(bounds.in_col_orig, bounds.in_col_end)):
#                    print 'WARNING: _output grid cell outside of range of input'
#                    sys.exit('IN: %s,%s   OUT: %s,%s' %(i_col,i_row,col,row))
                    continue
                data_out[:,:,row,col] = data_in[:,:,i_row,i_col]
        return data_out
        
    def _agg(self, data_in, grid, bounds):
        """
        Aggregate from a finer to a coarser grid
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.rows, grid.cols], np.float32)
        for col_cnt,col in enumerate(range(bounds.out_col_orig, bounds.out_col_end)):
            i_col = bounds.in_col_orig + (col_cnt * int(bounds.chunk_dim)) # Move over three old columns from the old grid origination point for each new column 
            for row_cnt,row in enumerate(range(bounds.out_row_orig, bounds.out_row_end)):
                i_row = bounds.in_row_orig + (row_cnt * int(bounds.chunk_dim))
                cell_out = np.zeros([data_in.shape[0],data_in.shape[1]], np.float32)
                # Loop over three cells for each dimension on the old grid            
                for x in range(int(bounds.chunk_dim)):
                    x_blk = i_col + x
                    if x_blk not in list(range(bounds.in_col_orig, bounds.in_col_end)): 
                        break
                    for y in range(int(bounds.chunk_dim)):
                        y_blk = i_row + y 
                        if y_blk not in list(range(bounds.in_row_orig, bounds.in_row_end)): 
                            break
#                        print 'In (col,row): %s,%s   Out: %s,%s' %(x_blk, y_blk, col, row)
                        cell_out[:,:] += data_in[:,:,y_blk,x_blk]
                data_out[:,:,row,col] = cell_out
        return data_out

    def _sub_div(self, data_in, grid, bounds):
        """
        Divide and window to a grid of smaller size
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.rows, grid.cols], np.float32)
        for col_cnt,col in enumerate(range(bounds.out_col_orig, bounds.out_col_end)):
            for row_cnt,row in enumerate(range(bounds.out_row_orig, bounds.out_row_end)):
                i_row = int(bounds.in_row_orig + (row_cnt * bounds.chunk_dim))
                i_col = int(bounds.in_col_orig + (col_cnt * bounds.chunk_dim))
                if i_row not in list(range(bounds.in_row_orig, bounds.in_row_end)) or i_col not in list(range(bounds.in_col_orig, bounds.in_col_end)):
#                    print 'WARNING: _output grid cell outside of range of input'
#                    sys.exit('IN: %s,%s   OUT: %s,%s' %(i_col,i_row,col,row))
                    continue
                data_out[:,:,row,col] = data_in[:,:,i_row,i_col] * (bounds.chunk_dim * bounds.chunk_dim)
        return data_out

    def _write_tflag(self, in_ncf):
        '''
        Copy time flag information from input file
        '''
        species_in = in_ncf.variables['TFLAG']
        data_in = species_in[:]
        species_out = self.createVariable('TFLAG', np.int32, ('TSTEP','VAR','DATE-TIME'))
        species_out[:] = data_in
        species_out.long_name = 'TFLAG'
        species_out.units = '<YYYYDDD,HHMMSS>'
        species_out.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'
  
def main():
    if len(sys.argv) != 4:
            print("You must provide an input file name, out file name, output grid name")
            print("regrid.py <infile> <outfile> <grid name>")
            sys.exit()

    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]
    out_grid = sys.argv[3]

    out_grid = Grid(out_grid, def_griddesc)
    with NCF(infile_name) as in_ncf, NCF(outfile_name, 'w') as out_ncf:
        out_ncf.regrid(in_ncf, out_grid)

if __name__ == '__main__':
	main()
