from __future__ import division
from builtins import object
from math import floor
import numpy as np
import pandas as pd
from pyproj import Proj

class Grid(object):
    """
    Reads the grid description file and loads the grid information for the specified grid named.
    """
    def __init__(self, grid_name, grid_desc):
        self.GDNAM = grid_name
        self.load_gridinfo(grid_desc)
        self._set_proj()

    def _parse_float(self, x):
        """
        Returns a floating point with the correct number of trailing zeros based on the .Dx
        """
        x = x.replace('D','E') 
        return np.float64(x)

    def _split_line(self, line):
        '''
        The grid description file is terminated with a ! and separated with ","
        '''
        return [cell.strip().strip("'") for cell in line.strip().split('!')[0].split(',')]

    def load_gridinfo(self, grid_desc):
        """
        Read in the grid description file and store the grid data as object attributes
        The grid description file is multi-part containing inconsistently ordered information for
          grids and projections.
        See the IOAPI or SMOKE documentation for information on the grid description file
        """
        with open(grid_desc) as gd:
            # Start out by reading all of the projections in the file
            state = 'proj'
            proj_table = dict()
            for line in gd:
                s_line = self._split_line(line)
                if state == 'proj':
                    if s_line[0]:
                        if s_line[0] == ' ':
                            state = 'grid'
                        else:
                            # Read the projection parameters
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
                    # Read through the grids until the specified grid name is found
                    if s_line[0] == self.GDNAM:
                        line = next(gd)
                        s_line = self._split_line(line)
                        proj_name = s_line[0]
                        self.XORIG, self.YORIG, self.XCELL, self.YCELL = \
                          [self._parse_float(x) for x in s_line[1:5]]
                        self.NCOLS, self.NROWS, self.NTHIK = [int(x) for x in s_line[5:8]]
                        for k, v in list(proj_table[proj_name].items()):
                            setattr(self, k, v)
                        state = 'done'
                        break
            if state == 'grid':
                raise ValueError('Grid %s not found in grid description file.' %self.GDNAM)

    def _set_proj(self):
        '''
        Define the projection
        Only included LCC and polar sterographic. Really should only use LCC for these.
        LCC parameters follow typical SMOKE grid definitions including a radius of 6370 km
        '''
        if self.GDTYP == 6:
            self.proj = Proj(proj='stere', lat_ts=self.P_BET, lat_0=self.YCENT, lon_0=self.P_GAM)
        elif self.GDTYP == 2:
            self.proj = Proj(proj='lcc', lat_1=self.P_ALP, lat_2=self.P_BET, lon_0=self.P_GAM, lat_0=self.YCENT, a=6370000, b=6370000, no_defs=True)

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
        return (abs(x_1-x_0)**2. + abs(y_1-y_0)**2.)**(0.5)

    def get_cell(self, x, y):
        '''
        Get the grid col and row of the projected (x,y) coordinates
        '''
        col = int(floor((x-self.XORIG)/self.XCELL) + 1)
        row = int(floor((y-self.YORIG)/self.YCELL) + 1)
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
        cells['x'] = ((col.astype('f') - 1) * self.XCELL) + self.XORIG
        cells['y'] = ((row.astype('f') - 1) * self.YCELL) + self.YORIG
        df = cells.apply(lambda row: self.inverse_proj(*row), axis=1)
        return pd.DataFrame(df.values.tolist(), columns=['lon','lat'], index=df.index)
        
    def colrow_to_coords(self, col, row):
        '''
        Get the lambert SW corner from the column and row
        '''
        x = ((col - 1) * self.XCELL) + self.XORIG
        y = ((row - 1) * self.YCELL) + self.YORIG
        return x,y


