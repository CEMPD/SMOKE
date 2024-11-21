from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import numpy as np
import netCDF4 as ncf
from datetime import date, datetime, timedelta
from qamods.helpers import conv2jul, get_grid, parse_ratio

class data_out(object):
    """
    Outfile class.
    """

    def __init__(self, outfile_name, out_type, gsdate, verbosity = False):
        self.outfile_name = outfile_name
        self.out_type = out_type.upper()
        self.gsdate = gsdate
        self.verbosity = verbosity    
    
        if self.verbosity == True: 
            print('Opening %s file for writing: %s' %(self.out_type, self.outfile_name))
        if self.out_type == 'NCF':
            self.outfile = self._open_NCF()
        elif self.out_type == 'CSV':
            self.outfile = self._open_CSV()
        else:
            raise ValueError('Wrong outfile type specified.')

    def __str__(self):
        return self.outfile_name

    def __call__(self):
        return self.outfile

    def write_outfile(self, out_dict, grid = '', grid_desc = '', region = '', tons = False, units = '',
       srg_file=''):
        species_list = list(out_dict.keys())
        while '0' in species_list:  # Dunno why this gets inserted?
            species_list.remove('0')
        self.species_list = sorted(species_list)

        if self.out_type == 'NCF':
            self._write_NCF(out_dict, self.species_list, grid, grid_desc, tons, units)
        elif self.out_type == 'CSV':
            if region == 'state' or region == 'county':
                self._write_FIPS_CSV(out_dict, self.species_list, grid, grid_desc, tons, region, srg_file)
            else: 
                self._write_grid_CSV(out_dict, self.species_list, tons)

    def _open_CSV(self):
        try: 
            outfile = open(self.outfile_name, 'w')
        except IOError:
            raise IOError('%s not available for access.' %self.outfile_name)
        else: 
            return outfile 

    def _write_grid_CSV(self, out_dict, species_list, tons):
        '''
        Writes the dictionary to an output file in csv format by grid cell.  Takes output dictionary and the species_list.
        '''
        if out_dict[species_list[0]]().shape[0] == 1:
            hour_type = 'sum'
        else:
            hour_type = 'hourly' 
        self.outfile.write('hour,row,col,' + ','.join(species_name for species_name in species_list) + '\n')
        for hour in range(out_dict[species_list[0]]().shape[0]):
            if hour_type == 'sum':
                out_hour = 'sum'
            else:
                out_hour = hour
            for row in range(out_dict[species_list[0]]().shape[1]):
                srow = str(row)
                for col in range(out_dict[species_list[0]]().shape[2]):
                    scol = str(col) 
                    out_line = '%s,%s,%s' %(out_hour, row + 1, col + 1)
                    for species_name in species_list:
                        out_line = '%s,%s' %(out_line, out_dict[species_name]()[hour,row,col])
                    self.outfile.write(out_line + '\n')

    def _write_FIPS_CSV(self, out_dict, species_list, grid, grid_desc, tons, region, srg_file):
        '''
        Writes the dictionary to an output file in csv format by fips.  Takes output dictionary and the species_list.
        '''
        if not grid: 
            raise ValueError('No grid specified.  Grid needed to write state or county based csv.')
        grid_info = get_grid(grid, grid_desc)

        ratio_table = parse_ratio(region, grid, grid_desc, srg_file)
        self.outfile.write('hour,fips,' + ','.join(species_name for species_name in species_list) + '\n')
        fips_list = sorted(ratio_table.keys())
        if out_dict[species_list[0]]().shape[0] == 1:
            hour_type = 'sum'
        else:
            hour_type = 'hourly' 
        for hour in range(out_dict[species_list[0]]().shape[0]):
            if hour_type == 'sum':
                out_hour = 'sum'
            else:
                out_hour = hour
            for fips in fips_list:
                out_line = '%s,%s' %(out_hour, fips)
                line_dict = dict((species_name,float(0)) for species_name in species_list)
                for species_name in species_list:
                    for cell in ratio_table[fips]:
                        col = int(cell.split(',')[0])
                        row = int(cell.split(',')[1])
                        if row > grid_info['rows'] or col > grid_info['cols']:
                            continue  # Skip ratios that are outside the domain
                        line_dict[species_name] = line_dict[species_name] + (out_dict[species_name]()[hour,row - 1,col - 1] * ratio_table[fips][cell])
                    out_line = '%s,%s' %(out_line, line_dict[species_name])
                self.outfile.write(out_line + '\n') 

    def _open_NCF(self):
        '''
        Opens the netCDF input file and returns an open file object.
        '''
        try: 
            outfile = ncf.Dataset(self.outfile_name, 'w', format='NETCDF3_64BIT')
        except TypeError:
            raise IOError('%s not available for access.' %self.outfile_name)
        else: 
            return outfile 

    def _write_NCF(self, out_dict, species_list, grid, grid_desc, tons = False, units = ''):
        '''
        Writes the dictionary to an output file in NCF format.  Takes output dictionary and the species_list.
        '''
        hours = out_dict[species_list[0]]().shape[0]
        self._outfile_settings(hours, species_list, grid, grid_desc, self.gsdate)
        self._create_TFLAG(hours, self.gsdate)
        for species_name in species_list:
            self._write_species(out_dict, species_name, tons, units)
            
    def _write_species(self, out_dict, species_name, tons = False, units = '', long_name = '', var_desc = ''):
        """
        Takes the output dictionary, species name, and optionally long name, units, and variable description.
        Creates a species of name species_name with the standard smoke shape of TSTEP, LAY, ROW, and COL.
        Returns a multidimensional array of shape TSTEP, LAY, ROW, COL.
        Must set global file attributes before running.
        """
        if not units:
            if tons:
                units = 'tons/day'
            else:
                units = 'moles/s'

        if not long_name: 
            long_name = species_name
        if not var_desc: 
            var_desc = 'Model species ' + species_name
        d_shape = [out_dict[species_name]().shape[0],1,out_dict[species_name]().shape[1],out_dict[species_name]().shape[2]]
        species_out = self.outfile.createVariable(species_name, np.float32, ('TSTEP','LAY','ROW','COL'))
        data_out = np.zeros(d_shape, '>f4')
        species_out.long_name = long_name
        species_out.units = units
        species_out.var_desc = var_desc
        for lay in range(data_out.shape[1]):
                try:
                    data_out[:,lay,:,:] = out_dict[species_name]()
                except ValueError:
                    raise ValueError('Array size mismatch. Please check that input domain matches the size of the intended output domain (-G [GRID]).')
        species_out[:] = data_out

    def _outfile_settings(self, hours, species_list, grid, grid_desc, gsdate):
        '''
        Set the output file dimensions and IOAPI metadata
        '''
        import time
        grid_info = get_grid(grid, grid_desc)
        out_dims = {'VAR': len(species_list), 'TSTEP': hours, 'DATE-TIME': 2, 'ROW': grid_info['rows'], 
                    'COL': grid_info['cols'], 'LAY': 1}
        for dim, value in list(out_dims.items()): 
            self.outfile.createDimension(dim, value)
        # Set the time step based on how many hours are in the DS
        if hours > 1:
            hstep = 10000
        else:
            hstep = 240000
        esdate = conv2jul(gsdate)
        var_list = ''.join([species_name.ljust(16)[:16] for species_name in species_list])
        out_atts = {'IOAPI_VERSION': '$Id: @(#) ioapi library version 3.1 $'.ljust(80),
                    'EXEC_ID': '????????????????'.ljust(80), 'VAR-LIST': var_list, 'FTYPE': 1, 
                    'CDATE': int(time.strftime('%Y%j')), 'CTIME': int(time.strftime('%H%M%S')),
                    'WDATE': int(time.strftime('%Y%j')), 'WTIME': int(time.strftime('%H%M%S')),
                    'SDATE': esdate, 'STIME': 0, 'TSTEP': hstep, 'NTHIK': grid_info['nthik'], 
                    'NCOLS': grid_info['cols'], 'NROWS': grid_info['rows'],
                    'NLAYS': 1, 'NVARS': len(species_list), 'GDTYP': grid_info['gdtyp'],
                    'P_ALP': grid_info['palp'], 'P_BET': grid_info['pbet'], 'P_GAM': grid_info['pgam'], 
                    'XCENT': grid_info['xcent'], 'YCENT': grid_info['ycent'],
                    'XCELL': grid_info['xcell'], 'YCELL': grid_info['ycell'], 
                    'XORIG': grid_info['xorig'], 'YORIG': grid_info['yorig'],
                    'VGTYP': -1, 'VGTOP': np.zeros([1], '>f'),'VGLVLS': np.zeros([2], '>f'),
                    'GDNAM': grid_info['name'].ljust(16), 'UPNAM': 'METXTRACT'.ljust(16), 
                    'FILEDESC':  'EMISQA', 'HISTORY': ' '}
        for att_name, att in list(out_atts.items()): 
            setattr(self.outfile, att_name, att)

    def _create_TFLAG(self, hours, gsdate):
        """
        Create a new TFLAG and adjust to size for number of variables in the outfile based on an infile.
        """
        data_out = np.zeros([hours, len(self.species_list), 2], '>i')
        species_out = self.outfile.createVariable('TFLAG', np.int32, ('TSTEP', 'VAR', 'DATE-TIME'))
        cur_time = datetime.strptime(gsdate, '%Y%m%d')
        for tStep in range(data_out.shape[0]):
            for tVar in range(data_out.shape[1]):
                data_out[tStep,tVar,0] = datetime.strftime(cur_time, '%Y%j')      # Set date for TFLAG timestep
                data_out[tStep,tVar,1] = datetime.strftime(cur_time, '%H%M%S')    # Set time for TFLAG timestep
            cur_time += timedelta(hours=1)
        species_out[:] = data_out
        species_out.long_name = 'TFLAG'
        species_out.units = '<YYYYDDD,HHMMSS>'
        species_out.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

