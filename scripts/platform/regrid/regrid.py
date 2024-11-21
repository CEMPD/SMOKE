#!/work/EMIS/python/miniforge3/envs/smoke_python/bin/python
# Intelligently regrids from one domain to another for b3grd
# Lambert to Lambert only at this time
# 4/12/13 James Beidler <beidler.james@epa.gov>

import sys
from datetime import datetime, timedelta
from math import ceil
import numpy as np
import fauxioapi as io

# Set the default grid description file.  
def_griddesc = '/work/EMIS/em_v8.1/ge_dat/gridding/griddesc_lambertonly_05jun2023_nf_v29.txt'

########

def main():
    if len(sys.argv) != 4:
            print("You must provide an input file name, out file name, output grid name")
            print("regrid.py <infile> <outfile> <grid name>")
            sys.exit()
    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]
    out_grid = sys.argv[3]
    out_grid = io.Grid(out_grid, def_griddesc)
    with NCF(infile_name) as in_ncf, NCF(outfile_name, 'w') as out_ncf:
        out_ncf.regrid(in_ncf, out_grid)
        out_ncf.set_attributes(in_ncf.SDATE, out_grid)
        out_ncf.write_TFLAG()

class GridBounds:
    '''
    Set the boundaries of the new grid and the old grid based on the new grid
    '''
    def __init__(self, in_ncf, grid):
        self._set_bounds(in_ncf, grid)

    def _set_bounds(self, in_ncf, grid):
        # Calculate the raw difference for the distance between the origination points on the new and old grids
        x_dist = (grid.XORIG - in_ncf.XORIG)
        y_dist = (grid.YORIG - in_ncf.YORIG)
        # Calculate the starting points for copying on both grids
        # Keep in mind that the grid starts at the SW corner
        if x_dist < 0: # If the old x origination is inside the new grid
            self.in_col_orig = 0 # Start at the first old grid cell 
            self.out_col_orig = int(abs(x_dist/grid.XCELL)) # Offset some new grid cells
        else: # If the old x origination is outside of equal to the orig on new grid
            self.in_col_orig = int(abs(x_dist/in_ncf.XCELL))
            self.out_col_orig = 0
        if y_dist < 0:
            self.in_row_orig = 0
            self.out_row_orig = int(abs(y_dist/grid.YCELL))
        else:
            self.in_row_orig = int(abs(y_dist/in_ncf.YCELL))
            self.out_row_orig = 0
        # Calculate the end points to form boundary fields
        in_x_end = in_ncf.XORIG + (in_ncf.NCOLS * in_ncf.XCELL)
        in_y_end = in_ncf.YORIG + (in_ncf.NROWS * in_ncf.YCELL)
        out_x_end = grid.XORIG + (grid.NCOLS * grid.XCELL)
        out_y_end = grid.YORIG + (grid.NROWS * grid.YCELL)
        x_dist = out_x_end - in_x_end
        y_dist = out_y_end - in_y_end
        if x_dist > 0: # If the endpoint of the old is inside the new
            self.in_col_end = in_ncf.NCOLS
            self.out_col_end = int(grid.NCOLS - abs(x_dist/grid.XCELL)) # Calculate the last grid cell to aggregate
            if self.out_col_end > int(self.out_col_end): # Check to see if there is anything after the decimal ie. a partial column calculated.
                self.out_col_end = int(self.out_col_end) + 1  # If there is a partial column then add a column to be processed to hold the extra input data
        else:
            self.in_col_end = int(in_ncf.NCOLS - abs(x_dist/in_ncf.XCELL))
            self.out_col_end = grid.NCOLS
        if y_dist > 0:
            self.in_row_end = in_ncf.NROWS
            self.out_row_end = int(grid.NROWS - abs(y_dist/grid.YCELL))
            if self.out_row_end > int(self.out_row_end):
                self.out_row_end = int(self.out_row_end) + 1
        else:
            self.in_row_end = int(in_ncf.NROWS - abs(y_dist/in_ncf.YCELL))
            self.out_row_end = grid.NROWS
    #        print "In start col: %s  End: %s    Out start: %s  End: %s" %(self.in_col_orig, self.in_col_end, self.out_col_orig, self.out_col_end) 
    #        print "In start row: %s  End: %s    Out start: %s  End: %s" %(self.in_row_orig, self.in_row_end, self.out_row_orig, self.out_row_end)
        self.xchunk = float(grid.XCELL)/float(in_ncf.XCELL)  # How many times larger is the output cell versus the input
        self.ychunk = float(grid.YCELL)/float(in_ncf.YCELL)  # How many times larger is the output cell versus the input

    def set_slices(self):
        '''
        Set the numpy ranges
        '''
        self.row_cnt = len(range(self.out_row_orig, self.out_row_end))
        self.in_rows = (self.in_row_orig, self.in_row_orig + ceil(self.row_cnt * self.xchunk))
        self.out_rows = (self.out_row_orig, self.out_row_end)
        self.col_cnt = len(range(self.out_col_orig, self.out_col_end))
        self.in_cols = (self.in_col_orig, self.in_col_orig + ceil(self.col_cnt * self.ychunk))
        self.out_cols = (self.out_col_orig, self.out_col_end)

class NCF(io.IODataset):
    """
    NCF subclass
    """
    def __init__(self, file_name, mode='r'):
        print('Opening %s' %file_name)
        io.IODataset.__init__(self, file_name, mode)

    def regrid(self, in_ncf, grid_info):
        """
        the main regridding section
        sets loop over variables and decides how to regrid
        """ 
        grid = grid_info
        var_list = getattr(in_ncf, 'VAR-LIST')
        var_list = [var_list[x*16:(x+1)*16].strip() for x in range(int(in_ncf.NVARS))] 
        self.set_dimensions(LAY=in_ncf.NLAYS, ROW=grid.NROWS, COL=grid.NCOLS, VAR=len(var_list))
        bounds = GridBounds(in_ncf, grid)
        bounds.set_slices()
        if in_ncf.XCELL < grid.XCELL:
            # Aggregate
            regrid_method = self.agg
            print('  Aggregating %s to %s' %(in_ncf.GDNAM, grid.GDNAM))
        elif in_ncf.XCELL > grid.XCELL:
            # Subdivide
            regrid_method = self.subdiv
            print('  Sub-dividing %s to %s' %(in_ncf.GDNAM, grid.GDNAM))
        else:
            # Window or expand
            regrid_method = self.subset
            print('  Subsetting %s to %s' %(in_ncf.GDNAM, grid.GDNAM))
        print('Regridding %s to %s %s' %(in_ncf.filepath(), self.filepath(), grid.GDNAM))
        cpl = 0
        c_max = float(len(var_list))
        sys.stdout.write("\r%.2f%%" %(float(cpl)/c_max))
        # Loop through and subset each species variable
        for species_name in var_list:
            species_in = in_ncf.variables[species_name]
            data_in = species_in[:]
            species_out = self.create_variable(species_name, 'REAL', ('TSTEP', 'LAY', 'ROW', 'COL'))
            species_out.long_name = species_in.long_name
            species_out.units = species_in.units
            species_out.var_desc = species_in.var_desc
            species_out[:] = regrid_method(data_in, grid, bounds) 
            self.sync()
            cpl += 1
            sys.stdout.flush()
            sys.stdout.write("\r%.2f%%" %((float(cpl)/c_max)*100.))
        sys.stdout.write("\n")

    def subset(self, data_in, grid, bounds):
        """
        Window to a grid or expand to a grid of equal cell size
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.NROWS, grid.NCOLS], np.float32)
        data_out[:,:,slice(*bounds.out_rows),slice(*bounds.out_cols)] = data_in[:,:,slice(*bounds.in_rows),slice(*bounds.in_cols)]
        return data_out

    def agg(self, data_in, grid, bounds):
        """
        Aggregate from a finer to a coarser grid
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.NROWS, grid.NCOLS], np.float32)
        # reshape and aggregate the input
        out_shape = (data_in.shape[0], data_in.shape[1], 
          (bounds.out_rows[1]-bounds.out_rows[0]), int(bounds.xchunk),
          (bounds.out_cols[1]-bounds.out_cols[0]), int(bounds.ychunk))
        data_in = data_in[:,:,slice(*bounds.in_rows),slice(*bounds.in_cols)].reshape(out_shape).sum(axis=(3,5))
        data_out[:,:,slice(*bounds.out_rows),slice(*bounds.out_cols)] = data_in[:]
        return data_out

    def subdiv(self, data_in, grid, bounds):
        """
        Divide and window to a grid of a smaller size
        """
        data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.NROWS, grid.NCOLS], np.float32)
        data_in = data_in[:,:,slice(*bounds.in_rows),slice(*bounds.in_cols)] * bounds.xchunk * bounds.ychunk
        data_in = data_in.repeat(1 // bounds.xchunk, axis=2).repeat(1 // bounds.ychunk, axis=3)
        data_out[:,:,slice(*bounds.out_rows),slice(*bounds.out_cols)] = data_in[:,:,slice(*bounds.out_rows),slice(*bounds.out_cols)]
        return data_out

if __name__ == '__main__':
	main()
