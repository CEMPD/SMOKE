from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import numpy as np
import netCDF4 as ncf
import time
from datetime import date, datetime, timedelta

class IOVar(ncf.Variable):
    '''
    I/OAPI variable
    '''
    def __init__(self, ds, vname, dtype, dims, **kwargs):
        '''
        IOAPI dataset variable
        Variable Name
        Dtype : INT, REAL or DBLE
        Dimensions: ie ['TSTEP','LAY','ROW','COL']
        Keyword arguments to set attributes: long_name, units or var_desc 
        '''
        ioapi_dtypes = {'INT': np.int32, 'REAL': np.float32, 'DBLE': np.float64}
        np_dtype = ioapi_dtypes[dtype] 
        ncf.Variable.__init__(self, ds, vname, np_dtype, dims)
        self.long_name = vname.ljust(80)
        self.units = ''.ljust(80)
        self.var_desc = ''.ljust(80)
        for att, val in list(kwargs.items()):
            if att in ('long_name','units','var_desc'):
                setattr(self, att, val)

class IODataset(ncf.Dataset):
    '''
    I/OAPI dataset
    '''
    def __init__(self, fname, mode='r'):
        ncf.Dataset.__init__(self, fname, mode, format='NETCDF3_CLASSIC')

    def create_variable(self, vname, dtype, dims, **kwargs):
        '''
        Create an IOAPI variable for this DS
        '''
        new_var = IOVar(self, vname, dtype, dims, **kwargs)
        self.variables[vname] = new_var
        self.sync()
        return new_var

    def set_dimensions(self, ftype='GRID', **kwargs):
        '''
        Set the file dimensions
        Use None (nonetype) for Unlimited
        '''
        if ftype == 'GRID':
            dim_list = ('TSTEP','LAY','ROW','COL','VAR')
        else:
            raise ValueError('Unknown I/OAPI filetype. GRID only type currently supported.')
        for dim, val in list(kwargs.items()):
            if dim in dim_list:
                self.createDimension(dim, val)
        if 'TSTEP' not in list(kwargs.keys()):
            self.createDimension('TSTEP', None)
        self.createDimension('DATE-TIME', 2)

    def set_attributes(self, sdate, grid, **kwargs):
        # Set grid-based attributes from grid desc
        for att in grid.grid_atts:
            if att == 'GDNAM':
                val = grid.GDNAM.ljust(16)
            else:
                val = getattr(grid, att)
            setattr(self, att, val)
        var_list = [var for var in self.variables if var != 'TFLAG']
        self.NVARS = len(var_list)
        setattr(self, 'VAR-LIST', ''.join([var.ljust(16) for var in var_list]))
        # Currently only gridded / type 1 is supported
        self.FTYPE = 1
        self.STIME = 0
        self.SDATE = int(sdate)
        self.TSTEP = 10000  # One hour per step
        self.IOAPI_VERSION = 'FAKE IOAPI'.ljust(80)
        self.EXEC_ID = '?'.ljust(80)
        self.CDATE = self.WDATE = int(time.strftime('%Y%j'))
        self.CTIME = self.WTIME = int(time.strftime('%H%M%S'))
        # Vertical levels. Assume 1.
        self.NLAYS = len(self.dimensions['LAY'])
        self.VGTYP = -1
        self.VGTOP = np.zeros([1], np.float32)
        self.VGLVLS = np.zeros([2], np.float32)
        self.UPNAM = 'FI'.ljust(16)
        self.FILEDESC = 'FAKE IOAPI'.ljust(80)
        self.HISTORY = ' '
        # Override any attribute with a keyword arg
        for att, val in list(kwargs.items()):
            setattr(self, att, val) 

def write_TFLAG(ds):
    tflag = ds.create_variable('TFLAG', 'INT', ('TSTEP','VAR','DATE-TIME'),
      long_name='TFLAG',
      units='<YYYYDD,HHMMSS>',
      var_desc='Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS')
    arr = np.zeros(tflag.shape, np.int32)
    cur_time = datetime.strptime(str(ds.SDATE), '%Y%j')
    for tstep in range(arr.shape[0]):
        for tvar in range(arr.shape[1]):
            arr[tstep,tvar,0] = datetime.strftime(cur_time, '%Y%j')      # Set date for TFLAG timestep
            arr[tstep,tvar,1] = datetime.strftime(cur_time, '%H%M%S')    # Set time for TFLAG timestep
        cur_time += timedelta(hours=1)
    tflag[:] = arr 

class Grid(object):
    """
    Reads the grid description file and loads the grid information for the specified grid named.
    """
    def __init__(self, grid_name, grid_desc):
        self.GDNAM = grid_name
        self.load_gridinfo(grid_desc)
        self.grid_atts = ['GDTYP','GDNAM','GDTYP','P_ALP','P_BET','P_GAM','XCENT','YCENT','XORIG',
          'YORIG','XCELL','YCELL','NCOLS','NROWS','NTHIK']

    def _parse_float(self, x):
        """
        Returns a floating point with the correct number of trailing zeros based on the .Dx
        """
        x = x.replace('D','E') 
        return np.float64(x)

    def _split_line(self, line):
        return [cell.strip().strip("'") for cell in line.strip().split('!')[0].split(',')]

    def load_gridinfo(self, grid_desc):
        """
        Read in the grid description file and store the grid data as object attributes
        """
        with open(grid_desc) as gd:
            state = 'proj'
            proj_table = dict()
            for line in gd:
                s_line = self._split_line(line)
                if state == 'proj':
                    if s_line[0]:
                        if s_line[0] == ' ':
                            state = 'grid'
                        else:
                            proj_name = line.strip().strip("'")
                            line = next(gd)
                            s_line = self._split_line(line)
                            proj_table[proj_name] = {'GDTYP': int(s_line[0]),
                                'P_ALP': self._parse_float(s_line[1]),
                                'P_BET': self._parse_float(s_line[2]),
                                'P_GAM': self._parse_float(s_line[3]),
                                'XCENT': self._parse_float(s_line[4]),
                                'YCENT': self._parse_float(s_line[5])}
                else:
                    if s_line[0] == self.GDNAM:
                        line = next(gd)
                        s_line = self._split_line(line)
                        proj_name = s_line[0]
                        self.XORIG, self.YORIG, self.XCELL, self.YCELL = \
                          [self._parse_float(x) for x in s_line[1:5]]
                        self.NCOLS, self.NROWS, self.NTHIK = [int(x) for x in s_line[5:8]]
                        for k, v in list(proj_table[proj_name].items()):
                            setattr(self, k, v)
                        break
            if state == 1: 
                raise ValueError('Grid %s not found in grid description file.' %self.gdnam)


