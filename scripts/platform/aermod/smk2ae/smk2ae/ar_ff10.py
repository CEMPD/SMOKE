from __future__ import print_function
from builtins import object
import numpy as np
import pandas as pd

class AnnualFF10(object):
    '''
    Import the annual FF10 inventory files into emissions and stack dataframes 

    self.emis contains the emissions by region_cd, scc, and pollutant
    '''
    def __init__(self, inv_list, invtable, st_fips='', use_shapes=False):
        self.emis = pd.DataFrame() 
        self._st_fips = st_fips
        # Specifications of monthly data fields in the FF10
        self.mon_vals = ['jan_value','feb_value','mar_value','apr_value','may_value','jun_value',
          'jul_value','aug_value','sep_value','oct_value','nov_value','dec_value']
        # Specification of possible months
        self.mons = ['january','february','march','april','may','june','july','august','september',
          'october','november','december']
        # Roll up for the data to region_cd and scc. Use the shape_id too if the use_shapes flag is true.
        self._keys = ['region_cd','scc']
        if use_shapes:
            self._keys.append('shape_id')
        # Read in the columns that we think that we need
        self._usecols = self._keys+['ann_value',] + self.mon_vals
        self._dtype = {'region_cd': str, 'scc': str, 'poll': str, 'ann_value': 'f'}
        for mon in self.mon_vals:
            self._dtype[mon] = 'f'
        for inv_file in inv_list:
            self.load_inv(inv_file, invtable)
        self.emis.drop_duplicates(inplace=True)
        # Simple check to see if there are monthly values in the inventory
        month_sum = sum(self.emis[self.mons].sum(axis=1))
        if month_sum > 0:
            self.monthly = True
        else:
            self.monthly = False

    def merge_rungroups(self, xref):
        '''
        Merge in the rungroups by SCC
        '''
        self.emis = pd.merge(self.emis, xref[['scc','run_group']], on='scc', how='left')
        if not self.emis[(self.emis['run_group'].isnull()) | (self.emis['run_group'] == '')].empty:
            print('WARNING: Unmatched SCCs to run groups')
            print(self.emis.loc[self.emis['run_group'].isnull(), 'scc'].drop_duplicates())

    def get_emis(self, run_group):
        '''
        Get the emissions sum by region_cd and scc based on the run group 
        '''
        if 'run_group' in self.emis.index:
            run_emis = self.emis[self.emis['run_group'] == run_group].copy() 
            run_emis = run_emis[['run_group','region_cd','scc','ann_value']].copy().groupby(\
              ['run_group','region_cd','scc'], as_index=False, sort=False).sum()
        else:
            print('WARNING: Missing run groups. Merge run groups to emis first.')
            run_emis = pd.DataFrame()
        return run_emis

    def load_inv(self, inv_file, invtable):
        '''
        Import the individual FF10 inventory
        '''
        print('Loading %s' %inv_file)
        df = pd.read_csv(inv_file, comment='#', usecols=self._usecols+['poll',], dtype=self._dtype)
        df.rename(columns=dict(zip(self.mon_vals, self.mons)), inplace=True)
        # Zero pad the FIPS to 5 characters
        df['region_cd'] = df['region_cd'].str.zfill(5)
        # Split out the mode types
        df.loc[df['poll'].str.contains('__'), 'poll'] = df.loc[df['poll'].str.contains('__'), 
          'poll'].str.split('__').str[1]
        # Keep only the pollutants that are marked as kept in the inventory table
        df = pd.merge(df, invtable, on='poll', how='left')
        df = df[df['smoke_name'].notnull()].copy()
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        # Only keep emissions greater than 0. This reduces the size of helper files
        df = df[df['ann_value'] > 0].copy()
        df.drop(['poll','spec_factor'], inplace=True, axis=1)
        if self._st_fips:
            fips_list = [st.strip() for st in self._st_fips.strip().split(',')]
            df = df[df['region_cd'].str[:2].isin(fips_list)].copy()
        if not df.empty:
            # Pandas cannot aggregate with a NaN key. Assign all NaNs to -99.
            df[self._keys+['smoke_name',]] = df[self._keys+['smoke_name',]].fillna(-99)
            df = df.groupby(self._keys+['smoke_name',], as_index=False, sort=False).sum()
            '''
            Populate the emissions dataframe by facility_id, rel_point, scc, and pollutant
            '''
            self.emis = pd.concat((self.emis, df[self._keys+['ann_value','smoke_name']+self.mons]))
              

