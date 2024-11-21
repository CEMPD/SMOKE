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
import sys
import time
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from ioapi import *

# Set the default grid description file.  
def_griddesc = '/work/EMIS/users/bte/WO144.3_hemi/regrid_hemi/griddesc.txt'

########


def regrid(in_ncf, out_ncf, out_grid, matrix):
    """
    the main regridding section
    sets loop over variables and decides how to regrid
    Collapse layers
    """ 
    var_list = getattr(in_ncf, 'VAR-LIST')
    species_list = [var_list[x*16:(x+1)*16].strip() for x in range(int(in_ncf.NVARS))]
    out_ncf.set_dimensions(LAY=1, ROW=out_grid.NROWS, 
      COL=out_grid.NCOLS, VAR=len(species_list))
    print('Regridding %s to %s %s' %(in_ncf.filepath(), out_ncf.filepath(), out_grid.GDNAM))
    cpl = 0
    c_max = float(len(species_list))
    sys.stdout.write("\r%.2f%%" %(old_div(float(cpl),c_max)))
    # Loop through and subset each species variable
    for species_name in species_list:
        species_in = in_ncf.variables[species_name]
        data_in = species_in[:]
        species_out = out_ncf.create_variable(species_name, 'REAL', ('TSTEP', 'LAY', 'ROW', 'COL'))
        species_out.long_name = species_in.long_name
        species_out.units = species_in.units
        try:
            species_out.var_desc = species_in.var_desc
        except AttributeError:
            species_out.var_desc = species_in.long_name
        # Collapse layers
        data_out = np.sum(data_in, axis=1, keepdims=True)
        # Apply regridding matrix
        data_out = apply_matrix(data_out, out_grid, matrix)
        species_out[:] = data_out
        cpl += 1
        sys.stdout.flush()
        sys.stdout.write("\r%.2f%%" %((old_div(float(cpl),c_max))*100.))
    sys.stdout.write("\n")
    out_ncf.set_attributes(in_ncf.SDATE, out_grid)
    write_TFLAG(out_ncf)

def apply_matrix(data_in, grid, matrix):
    data_out = np.zeros([data_in.shape[0], data_in.shape[1], grid.NROWS, grid.NCOLS], np.float32)
    for grid_row in matrix:
        row = grid_row[0].astype(int) - 1
        col = grid_row[1].astype(int) - 1
        inrow = grid_row[2].astype(int) - 1
        incol = grid_row[3].astype(int) - 1
        data_out[:,:,row,col] += data_in[:,:,inrow,incol] * grid_row[4]
    return data_out

def get_matrix(matrix_file):
    matrix_df = pd.read_csv(matrix_file, sep=' ', comment='#', names=['row','col','inrow','incol','frac'], 
      skipinitialspace=True, dtype={'row':'i','col':'i','inrow':'i','incol':'i','frac':'f'})
    '''
    sum_df = matrix_df[['row','col','frac']].copy().groupby(['row','col'], as_index=False).sum()
    sum_df.rename(columns={'frac':'sum'},inplace=True)
    matrix_df = pd.merge(matrix_df, sum_df, on=['row','col'], how='left')
    # Normalize the fractions of input to output factors to get a weighting factor for the percentages
    matrix_df['frac'] = old_div(matrix_df['frac'],matrix_df['sum'])
    matrix_df.to_csv('matrix.csv', index=False)
    '''
    return matrix_df[['row','col','inrow','incol','frac']].as_matrix()

def main():
    if len(sys.argv) != 5:
            print("You must provide an input file name, out file name, output grid name")
            print("regrid.py <infile> <outfile> <grid name> <matrix file>")
            sys.exit()
    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]
    out_grid = sys.argv[3]
    matrix_file = sys.argv[4]
    out_grid = Grid(out_grid, def_griddesc)
    matrix = get_matrix(matrix_file)
    with IODataset(infile_name) as in_ncf, IODataset(outfile_name, 'w') as out_ncf:
        regrid(in_ncf, out_ncf, out_grid, matrix)

if __name__ == '__main__':
	main()



