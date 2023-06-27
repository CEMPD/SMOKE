from __future__ import print_function
from builtins import object
import pandas as pd

class AnnualFF10(object):
    '''
    Import the annual FF10 inventory files into emissions and stack dataframes 
    '''
    def __init__(self, sector, inv_list, invtable, st_fips=''):
        self.st_fips = st_fips
        self.state = 'all'
        self._keys = ['facility_id','facility_name','rel_point_id','longitude','latitude','scc',
          'erptype','stkhgt','stkdiam','stktemp','stkflow','stkvel','fug_height','fug_width_ydim',
          'fug_length_xdim','fug_angle','region_cd','fac_source_type','unit_id','process_id']
        self._usecols = self._keys+['ann_value','poll']
        self._dtype = {'facility_id': str, 'facility_name': str, 'rel_point_id': str, 
          'region_cd': str, 'scc': str, 'poll': str, 'fac_source_type': str, 
          'process_id': str, 'unit_id': str}
        self.fips_list = []
        self.emis = pd.DataFrame()
        self.stk = pd.DataFrame()
        self.sector = sector
        for inv_file in inv_list:
            self.load_inv(inv_file, invtable)
        self.refresh_fac_list()
        self.scc_list = list(self.emis['scc'].drop_duplicates())

    def refresh_fac_list(self):
        '''
        Update the facilities list
        '''
        self.fac_list = list(self.stk['facility_id'].drop_duplicates())

    def refresh_fips_list(self, fips_list):
        '''
        Update the fips list
        '''
        for fips in fips_list:
            if fips not in self.fips_list:
                self.fips_list.append(fips)

    def load_inv(self, inv_file, invtable):
        '''
        Import the individual FF10 inventory
        '''
        print('Loading %s' %inv_file)
        df = pd.read_csv(inv_file, comment='#', usecols=self._usecols, dtype=self._dtype)
        df = pd.merge(df, invtable, on='poll', how='left')
        df = df[df['smoke_name'].notnull()].copy()
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        df = df[df['ann_value'] > 0].copy()
        df.drop(['poll','spec_factor'], inplace=True, axis=1)
        df['region_cd'] = df['region_cd'].str.zfill(5)
        if self.st_fips:
            df = df[df['region_cd'].str.startswith(self.st_fips)].copy()
        self.refresh_fips_list(df['region_cd'].drop_duplicates())
        # Pandas cannot aggregate with a NaN key. Assign all NaNs to -99.
        df.fillna(-99., inplace=True)
        df = df.groupby(self._keys+['smoke_name',], as_index=False).sum()
        df.loc[(df['scc'].isin(('2275050011','2275060011'))) & (df['smoke_name'] == 'LEAD'), 'ann_value'] = \
          df.loc[(df['scc'].isin(('2275050011','2275060011'))) & (df['smoke_name'] == 'LEAD'), 'ann_value'] * 0.5
        '''
        Populate the emissions dataframe by facility_id, rel_point, scc, and pollutant
        '''
        self.emis = pd.concat((self.emis, 
            df[['facility_id','unit_id','rel_point_id','scc','smoke_name','ann_value']]))
        '''
        Populate the stack dataframe with the stack parameters and the total HAPs emissions
        '''
        df.drop('smoke_name',axis=1,inplace=True)
        df = df.groupby(self._keys, as_index=False).sum()
        self.stk = pd.concat((self.stk,df))

    def drop_facilities(self, fac_list):
        '''
        Drop facilities
        '''
        self.emis = self.emis[~ self.emis['facility_id'].isin(fac_list)]
        self.stk = self.stk[~ self.stk['facility_id'].isin(fac_list)]
        self.refresh_fac_list()


