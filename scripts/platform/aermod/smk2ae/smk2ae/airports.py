from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import os.path
from math import ceil, log10
import pandas as pd
from smk2ae.utm import UTM
from smk2ae.temporal import match_temporal

class Airports(object):
    def __init__(self, airports):
        self.locs = self._read_airport_locs(airports)
        self.refresh_fac_list(self.locs['facility_id'].drop_duplicates())

    def strip(self, x):
        '''
        Basic return to strip leading/trailing whitespace from a string
        '''
        try:
            return x.strip()
        except AttributeError:
            return x

    def _read_airport_locs(self, airports):
        '''
        Read the airport locations file
        '''
        usecols = ['width','facility_id','start_x','start_y','end_x','end_y']
        df = pd.read_csv(airports, usecols=usecols, dtype={'facility_id': '|S16'})
        for fac in list(df['facility_id'].drop_duplicates()):
            n_srcs = len(df[df['facility_id'] == fac])
            srcs = ['AP%s' %str(int(x)).zfill(int(ceil(log10(n_srcs)))) for x in range(1,n_srcs+1)]
            df.loc[df['facility_id'] == fac, 'src_id'] = srcs
        return df
 
    def refresh_fac_list(self, fac_list):
        '''
        Update the list of airport facilities in the object
        '''
        self.fac_list = fac_list

    def refresh_locs(self):
        '''
        Refresh the runway locations
        '''
        self.locs = self.locs[self.locs['facility_id'].isin(self.fac_list)].copy()
        self.runways = list(self.locs['facility_id'].drop_duplicates())

    def match_airport_inv(self, inv):
        '''
        Match the airport facilities to the inventory
        Only keep the airport facilities that are also in the inventory
        '''
        self.inv = inv.stk[inv.stk['fac_source_type'] == '100'].copy()
        inv_fac_list = list(self.inv['facility_id'].drop_duplicates())
        self.refresh_fac_list(inv_fac_list)
        print('Total airport facilities: %s' %len(inv_fac_list))
        self.emis = inv.emis.loc[inv.emis['facility_id'].isin(self.fac_list), 
          ['facility_id','smoke_name','ann_value']].groupby(['facility_id','smoke_name'], 
          as_index=False).sum()
        self.inv = self.inv[self.inv['facility_id'].isin(self.fac_list)].copy() 
        self.refresh_locs()

    def write_runways(self, grid, temp):
        '''
        Write the files for airports with runways 
        '''
        inv_cols = ['facility_id','facility_name','fac_col','fac_row','state']
        df = pd.merge(self.locs, self.inv[inv_cols].drop_duplicates(), on='facility_id', how='left')
        if not df.empty:
            print('Airports facilities with runways: %s' %len(df['facility_id'].drop_duplicates()))
            df['length'] = df.apply(lambda row: grid.calc_dist(\
                row['start_x'],row['start_y'],row['end_x'],row['end_y']), axis=1)
            df['area'] = df['length'] * df['width']
            df['lnemis'] = old_div(1000., df['area'])
            df['relhgt'] = 3.
            df['szinit'] = 3.
            df['type'] = 'LINE'
            df['centroid_lon'] = grid.colrow_to_ll(df['fac_col'].astype('f') + 0.5, 
              df['fac_row'].astype('f') + 0.5)['lon'].astype('f')
            utm = UTM()
            df['utm_zone'] = df['centroid_lon'].apply(utm.get_zone)
            df[['xs1','ys1']] = df[['start_x','start_y','utm_zone']].apply(lambda x: pd.Series(\
              utm.get_coords(x['start_x'],x['start_y'],x['utm_zone'])), axis=1)
            df[['xs2','ys2']] = df[['end_x','end_y','utm_zone']].apply(lambda x: pd.Series(\
              utm.get_coords(x['end_x'],x['end_y'],x['utm_zone'])), axis=1)
            self._write_runway_locs(df)
            self._write_runway_params(df)
            self._write_temp_prof(df, temp, 'line')

    def _write_runway_locs(self, df):
        '''
        Write the locations file for airports with runways
        '''
        cols = ['state','facility_id','facility_name','src_id','xs1','ys1','xs2','ys2','utm_zone',
          'fac_col','fac_row']
        df = df[cols].copy().drop_duplicates(['facility_id','src_id']).sort()
        fname = os.path.join(os.environ['WORK_PATH'],'locations','airport_line_locations.csv') 
        df.to_csv(fname, index=False, quotechar=' ')

    def _write_runway_params(self, df):
        '''
        Write the release parameters file for airports with runways
        This is the "params" file
        '''
        cols = ['facility_id','facility_name','src_id','type','area','relhgt','width','szinit']
        df = df[cols].copy().sort()
        runway_cnt = df[['facility_id','src_id']].drop_duplicates().groupby('facility_id').size()
        runway_cnt = pd.DataFrame(runway_cnt, columns=['fract']).reset_index()
        df = pd.merge(df, runway_cnt, on='facility_id', how='left')
        df['fract'] = old_div(1., df['fract'])
        df.drop_duplicates(['facility_id','src_id'], inplace=True)
        fname = os.path.join(os.environ['WORK_PATH'],'parameters','airport_line_params.csv')
        df.to_csv(fname,  index=False, quotechar=' ',
          columns=['facility_id','facility_name','src_id','type','area','fract','relhgt','width','szinit'])

    def _write_temp_prof(self, df, temp, airport_type):
        '''        
        Write the temporal factors by source ID for each facility
        Adjust the output columns based on the max number of scalars used within the facility
        '''
        cols = ['facility_id','facility_name','src_id']
        df = df[cols].copy()
        inv_df = self.inv.loc[self.inv['facility_id'].isin(list(df['facility_id'].drop_duplicates())),
            ['region_cd','scc','facility_id','ann_value']].copy()
        inv_df.sort('ann_value', ascending=False, inplace=True)
        inv_df.drop_duplicates('facility_id', inplace=True)
        scalar_cols = [col for col in temp.columns if col.startswith('Scalar')]
        hierarchy = [['region_cd','scc','facility_id'],['scc','facility_id'],['region_cd','facility_id'],
        ['facility_id',], ['region_cd','scc'], ['scc',], ['region_cd',]]
        temp = match_temporal(inv_df, temp, ['qflag',]+scalar_cols, hierarchy)
        qflag_list = list(temp['qflag'].drop_duplicates())
        if 'MHRDOW7' in qflag_list:
            scalar_cols = ['Scalar%s' %x for x in range(1,2017)]
        elif 'MHRDOW' in qflag_list:
            scalar_cols = ['Scalar%s' %x for x in range(1,865)]
        elif 'HROFDY' in qflag_list:
            scalar_cols = ['Scalar%s' %x for x in range(1,25)]
        else:
            scalar_cols = ['Scalar%s' %x for x in range(1,13)]
        df = pd.merge(df, temp[['facility_id','qflag']+scalar_cols], on='facility_id', how='left') 
        df.drop_duplicates(['facility_id','src_id'], inplace=True)
        fname = os.path.join(os.environ['WORK_PATH'],'temporal','airport_%s_temporal.csv' %airport_type)
        df.to_csv(fname, index=False, quotechar=' ')

    def write_no_runways(self, grid, temp):
        '''
        Write the files for airports without runways 
        '''
        inv_cols = ['facility_id','facility_name','fac_col','fac_row','state','latitude',
          'longitude']
        df = self.inv.loc[~ self.inv['facility_id'].isin(self.runways), 
          inv_cols].copy()
        df.drop_duplicates(inplace=True)
        if not df.empty:
            print('Airports facilities without runways: %s' %len(df['facility_id'].drop_duplicates()))
            df['len_x'] = 10.
            df['len_y'] = 10.
            df['src_id'] = 'AP1'
            df['relhgt'] = 3.
            df['szinit'] = 3.
            df['angle'] = 0
            df['centroid_lon'] = grid.colrow_to_ll(df['fac_col'].astype('f') + 0.5, 
              df['fac_row'].astype('f') + 0.5)['lon'].astype('f')
            utm = UTM()
            df['utm_zone'] = df['centroid_lon'].apply(utm.get_zone)
            df[['utm_x','utm_y']] = df[['longitude','latitude','utm_zone']].apply(lambda x: pd.Series(\
              utm.get_coords(x['longitude'],x['latitude'],x['utm_zone'])), axis=1)
            df[['grid_x','grid_y']] = df[['longitude','latitude']].apply(lambda x: pd.Series(\
              grid.get_coords(x['longitude'],x['latitude'])), axis=1)
            self._write_nonrunway_locs(df)
            self._write_nonrunway_params(df)
            self._write_temp_prof(df, temp, 'nonrunway')

    def _write_nonrunway_locs(self, df):
        '''
        Write the unique locations of each source ID by facility
        '''
        cols = ['state','facility_id','facility_name','src_id','grid_x','grid_y','longitude',
            'latitude','utm_x','utm_y','utm_zone','fac_col','fac_row']
        df = df[cols].copy()
        fname = os.path.join(os.environ['WORK_PATH'], 'locations', 'airport_nonrunway_locations.csv')
        df.to_csv(fname, columns=cols, index=False, quotechar=' ')

    def _write_nonrunway_params(self, df):
        '''
        Write the unique params of each source ID by facility
        '''
        cols = ['facility_id','facility_name','src_id','relhgt','len_x','len_y','angle','szinit']
        df = df[cols].copy()
        fname = os.path.join(os.environ['WORK_PATH'], 'parameters', 'airport_nonrunway_params.csv')
        df.to_csv(fname, columns=cols, index=False, quotechar=' ')

    def calc_cell_xref(self, met_grid):
        '''
        Calculate the airport to met grid cell x-ref
        '''
        air_to_cell = calc_predominate_cell(self.inv, met_grid)
        air_to_cell = air_to_cell.drop_duplicates().sort()
        self.inv = pd.merge(self.inv, air_to_cell, on='facility_id', how='left')
        self.write_emis_xwalk()

    def write_emis_xwalk(self):
        '''
        Write the emissions x-walk
        '''
        df = pd.merge(self.emis, self.inv[['facility_id','facility_name','state']], 
          on='facility_id', how='left')
        df = df.groupby(['state','facility_id','facility_name','smoke_name'], as_index=False).sum()
        fname = os.path.join(os.environ['WORK_PATH'],'xwalk','airport_srcid_emis.csv')
        df.to_csv(fname, columns=['state','facility_id','facility_name','smoke_name','ann_value'],
            index=False, quotechar=' ')

def calc_predominate_cell(inv, grid_info):
    '''
    Calculate the predominate grid cell for each facility by weighting to the grid cell with
        the most emissions from a single facility
    '''
    df = inv[['facility_id','latitude','longitude','ann_value']].copy()
    df['fac_col'] = df.apply(lambda x: grid_info.ll_to_colrow(x['longitude'],x['latitude'])[0],
        axis=1)
    df['fac_row'] = df.apply(lambda x: grid_info.ll_to_colrow(x['longitude'],x['latitude'])[1],
        axis=1)
    df = df[['facility_id','fac_col','fac_row','ann_value']].groupby(['facility_id','fac_col',
        'fac_row'], as_index=False).sum() 
    df.sort(['facility_id','ann_value'], ascending=False, inplace=True)
    df.drop_duplicates('facility_id', inplace=True)
    return df[['facility_id','fac_col','fac_row']].copy() 

