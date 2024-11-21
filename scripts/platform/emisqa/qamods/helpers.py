from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
from datetime import timedelta, date, datetime
from qamods.chem_mechs import *
import os

# Helper procedures    
def checkEV(evName):
    """
    Checks if an environment variable is set.
    """
    try: 
        var = os.environ[evName]
    except KeyError: 
        return '' 
    else: 
        return var

def parse_float(x):
    """
    Returns a floating point with the correct number of trailing zeros based on the .Dx
    """
    from math import pow
    from numpy import float64
    if 'D' in x:    
        num = float(x.strip().split('D')[0])
        exp = int(x.strip().split('D')[1])
        return float64(num * pow(10, exp))
    else:
        return float64(x)

def conv2jul(gsdate):
    """
    Returns Julian date from Gregorian date.
    """
    gdate = datetime.strptime(str(gsdate), '%Y%m%d')
    return int(datetime.strftime(gdate, '%Y%j'))

def conv2greg(jul_date):
    """
    Returns Gregorian date from a Julian date.
    """
    jdate = datetime.strptime(str(jul_date), '%Y%j')
    return datetime.strftime(jdate, '%Y%m%d')

def moles2tons(val, species_name, in_format, mech = 'cmaq_cb05'):
    """
    Converts a value or array of values from moles/s to tons/hr
    """
    mechDct = molecDct[mech]

    if species_name in mechDct: 
        factor = mechDct[species_name]
    else: 
        factor = 1

    val = val * factor   # Convert moles per second to grams per second

    if in_format != 'UAM':
        val = val * 3600   # Convert tons per second to tons per hour

    val = val * (0.00000110231131)  # Convert grams per hour to tons per hour 

    return val  

def data_blocks(infile, size=1024):
    """
    Reads in binary files by blocks
    """
    while True:
        block = infile.read(size)
        if len(block) == 0:
            break
        yield block

def parse_ratio(region, grid, srg_file):
    """
    Parses in the county to grid conversion ratios.
    Returns a conversion dictionary of format { fips: { cell: ratio, ...}, ...}.
    """
    with open(srg_file) as infile:
        cell_size = grid.XCELL
        ratio_table = {}    
        for line in infile:
            line = [cell.strip() for cell in line.split('\t') if cell and cell != '!']
            if line[0].startswith('#'): 
                if line[0].strip() == '#GRID':
                    # Calculate grid cell offset from surrogate grid and cell range for our grid
                    xorig = float(line[2])
                    yorig = float(line[3])
                    col_offset = abs(int(old_div((xorig - grid.XORIG), cell_size)))
                    row_offset = abs(int(old_div((yorig - grid.YORIG), cell_size)))
                    col_range = range(1 + col_offset, grid.NCOLS + 1 + col_offset)
                    row_range = range(1 + row_offset, grid.NROWS + 1 + row_offset)
            else:
                row_dict = dict(list(zip(['code', 'fips', 'col', 'row', 'fac', 'cellarea', 'ctyarea', 'fac2'], line )))
                if region == 'state': 
                    fips = row_dict['fips'][:2]
                elif region == 'county': 
                    fips = row_dict['fips']

                # Check to see if the surrogate grid col and row is within the range of our grid
                if int(row_dict['col'])in col_range and int(row_dict['row']) in row_range:
                    # Offset the columns and rows.
                    col = int(row_dict['col']) - col_offset
                    row = int(row_dict['row']) - row_offset

                    cell = '%s,%s' %(col - 1, row - 1)
                    ratio = old_div(float(row_dict['cellarea']), (cell_size * cell_size))
                    ratio_table.setdefault(fips, {})
                    ratio_table[fips].setdefault(cell, 0)
                    ratio_table[fips][cell] += ratio
    return ratio_table

def clean_temp(zip_dict):
    """
    Cleans up the temporary zip output files at the end of a run.
    """
    if len(zip_dict) == 0: 
        return
    for name in zip_dict:
        os.remove(zip_dict[name])


