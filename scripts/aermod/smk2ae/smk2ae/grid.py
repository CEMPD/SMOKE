from math import floor
import numpy as np
import pandas as pd
from mpl_toolkits.basemap.pyproj import Proj

class Grid(object):
    """
    Reads the grid description file and loads the grid information for the specified grid named.
    """
    def __init__(self, grid_name, grid_desc):
        self.gdnam = grid_name
        self.grid_desc = grid_desc
        self.proj_dict = {'POL_HEMI': {'gdtyp': 6, 'p_alp': 1., 'p_bet': 45., 'p_gam': -98., 
            'xcent': -98., 'ycent': 90.},
            'LAM_40N97W': {'gdtyp': 2, 'p_alp': 33., 'p_bet': 45., 'p_gam': -97., 'xcent': -97., 'ycent': 40.},
            '12CONUS1': {'gdtyp': 2, 'p_alp': 33., 'p_bet': 45., 'p_gam': -97., 'xcent': -97., 'ycent': 40.}}
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
                    if self.gdnam.strip().upper() in line.split('!')[0].upper():
                        state = 2
                    continue
                elif state == 2:
                    split_line = line.split(',')
                    self.proj_name = split_line[0].strip("'")
                    self.xorig = self._parse_float(split_line[1])
                    self.yorig = self._parse_float(split_line[2])
                    self.xcell = self._parse_float(split_line[3])
                    self.ycell = self._parse_float(split_line[4])
                    self.ncols = int(split_line[5])
                    self.nrows = int(split_line[6])
                    self._set_proj()
                    break
            if state == 1: 
                raise ValueError, 'Grid %s not found in grid description file.' %self.gdnam

    def _set_proj(self):
        '''
        Define the projection
        '''
        if self.proj_name in self.proj_dict:
            for k,v in self.proj_dict[self.proj_name].iteritems():
                setattr(self, k, v)
        if self.gdtyp == 6:
            self.proj = Proj(proj='stere', lat_ts=self.p_bet, lat_0=self.ycent, lon_0=self.p_gam)
        elif self.gdtyp == 2:
            self.proj = Proj(proj='lcc', lat_1=self.p_alp, lat_2=self.p_bet, lon_0=self.p_gam, lat_0=self.ycent, a=6370000, b=6370000, no_defs=True)

    def get_coords(self, lon, lat):
        '''
        Get the projected (x,y) coordinates
        '''
        return self.proj(lon, lat)

    def calc_dist(self, lon_0, lat_0, lon_1, lat_1):
        '''
        Get the projected distance between two lat and lon pairs in native projection units 
        '''
        x_0, y_0 = self.get_coords(lon_0, lat_0)
        x_1, y_1 = self.get_coords(lon_1, lat_1)
        return (abs(x_1-x_0)**2. + abs(y_1-y_0)**2.)**(1./2.)

    def get_cell(self, x, y):
        '''
        Get the grid col and row of the projected (x,y) coordinates
        '''
        col = int(floor((x-self.xorig)/self.xcell) + 1)
        row = int(floor((y-self.yorig)/self.ycell) + 1)
        return col, row

    def ll_to_colrow(self,lon,lat):
        '''
        Abstraction for going directly from lat/lon to grid col and row
        '''
        x,y = self.get_coords(lon,lat)
        return self.get_cell(x,y)

    def inverse_proj(self, x, y):
        '''
        Apply inverse projection to get lat lon
        '''
        return self.proj(x, y, inverse=True)

    def colrow_to_ll(self, col, row):
        '''
        Get the lon and lat SW corner from the column and row
        '''
        cells = pd.DataFrame()
        cells['x'] = ((col.astype('f') - 1) * self.xcell) + self.xorig
        cells['y'] = ((row.astype('f') - 1) * self.ycell) + self.yorig
        df = cells.apply(lambda row: self.inverse_proj(*row), axis=1)
        return pd.DataFrame(df.tolist(), columns=['lon','lat'], index=df.index)
        
    def colrow_to_coords(self, col, row):
        '''
        Get the lambert SW corner from the column and row
        '''
        x = ((col - 1) * self.xcell) + self.xorig
        y = ((row - 1) * self.ycell) + self.yorig
        return x,y

