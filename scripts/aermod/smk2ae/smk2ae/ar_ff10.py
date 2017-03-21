import numpy as np
import pandas as pd

class AnnualFF10(object):
    '''
    Import the annual FF10 inventory files into emissions and stack dataframes 
    '''
    def __init__(self, inv_list, haps_list, st_fips='', use_shapes=False):
        self._st_fips = st_fips
        self.season_map = {'winter': ['jan_value','dec_value','feb_value'],
          'spring': ['mar_value','apr_value','may_value'], 
          'summer': ['jun_value','jul_value','aug_value'],
          'fall': ['sep_value','oct_value','nov_value']}
        self.mon_vals = [mon for season in self.season_map.itervalues() for mon in season]
        self._keys = ['region_cd','scc']
        if use_shapes:
            self._keys.append('shape_id')
        self._usecols = self._keys+['ann_value',] + self.mon_vals
        self._dtype = {'region_cd': '|S6', 'scc': '|S10', 'poll': '|S16', 'ann_value': 'f'}
        for mon in self.mon_vals:
            self._dtype[mon] = 'f'
        self.haps = pd.read_csv(haps_list, usecols=['poll','smoke_name','spec_factor'], dtype={'poll':'|S16'})
        self.inv_emis = pd.DataFrame()
        for inv_file in inv_list:
            self._load_inv(inv_file)
        self.inv_emis.drop_duplicates(inplace=True)
        season_sum = sum(self.inv_emis['summer']+self.inv_emis['winter']+self.inv_emis['spring']+\
          self.inv_emis['fall'])
        if season_sum > 0:
            self.seasons = True
        else:
            self.seasons = False
        self.emis = self.inv_emis[['region_cd','scc','ann_value']].copy().groupby(\
          ['region_cd','scc'], as_index=False, sort=False).sum()
        self.update_lists()

    def update_lists(self):
        self.scc_list = list(self.emis['scc'].drop_duplicates())
        self.fips_list = list(self.emis['region_cd'].drop_duplicates())
       
    def _load_inv(self, inv_file):
        '''
        Import the individual FF10 inventory
        '''
        print 'Loading %s' %inv_file
        df = pd.read_csv(inv_file, comment='#', usecols=self._usecols+['poll',], dtype=self._dtype)
        df['region_cd'] = df['region_cd'].str.zfill(6)
        df.ix[df['poll'].str.contains('__'), 'poll'] = df.ix[df['poll'].str.contains('__'), 
          'poll'].str.split('__').str[1]
        df = pd.merge(df, self.haps, on='poll', how='left')
        df = df[df['smoke_name'].notnull()].copy()
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        # Only keep emissions greater than 0
        df = df[df['ann_value'] > 0].copy()
        df.drop(['poll','spec_factor'], inplace=True, axis=1)
        if self._st_fips:
            fips_list = [st.strip() for st in self._st_fips.strip().split(',')]
            df = df[df['region_cd'].str[:3].isin(fips_list)].copy()
        if not df.empty:
            # Pandas cannot aggregate with a NaN key. Assign all NaNs to -99.
            df[self._keys+['smoke_name',]] = df[self._keys+['smoke_name',]].fillna(-99)
            df = df.groupby(self._keys+['smoke_name',], as_index=False, sort=False).sum()
            df = self._fill_seasons(df)
            '''
            Populate the emissions dataframe by facility_id, rel_point, scc, and pollutant
            '''
            self.inv_emis = pd.concat((self.inv_emis, 
              df[self._keys+['ann_value','smoke_name']+self.season_map.keys()]))

    def _fill_seasons(self, df):
        '''
        Fill the seasonal values
        '''
        for season, mon_list in self.season_map.iteritems():
            df[season] = df[mon_list].sum(axis=1, skipna=True)
            df.drop(mon_list, axis=1, inplace=True)
        return df
