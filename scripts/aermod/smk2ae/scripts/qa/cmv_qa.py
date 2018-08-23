#!/usr/bin/env python
'''
QA for AERMOD CMV helper files
requires python 3.5, pandas, and pyproj
'''
from __future__ import print_function
from past.utils import old_div 

from builtins import str
from builtins import object
import os
import string
from math import floor
import pandas as pd
from pyproj import Proj

class UTM(object):
    '''
    Define various functions for finding UTM coordinates form the lat and lon
    '''
    def __init__(self):
        self.projs={}

    def _define_proj(self, utm_zone):
        '''
        Define the UTM projection from specified UTM zone. Store it for later use.
        '''
        self.projs[utm_zone] = Proj(proj='utm',zone=utm_zone,ellps='WGS84',datum='WGS84',units='m')

    def get_coords(self, lon, lat, zone):
        '''
        Get the UTM coordinates from the lat/lon
        '''
#        zone = self.get_zone(lon,lat)
        if zone not in list(self.projs.keys()):
            self._define_proj(zone)
        return self.projs[zone](lon,lat)

    def get_zone(self, lon):
        '''
        Get the UTM zone from the lat
        '''
        zone = str(int(floor((old_div((lon+180),6))+1)))
        if zone not in list(self.projs.keys()):
            self._define_proj(zone)
        return zone

def get_inv_list():
    '''
    Get a list of all the inventory files from the environment variables
    '''
    inv_list = [] 
    for inv_var in ['EMISINV_%s' %alpha for alpha in string.ascii_uppercase]:
        if inv_var in os.environ:
            inv_list.append(os.environ[inv_var])
    return inv_list

def check_locs(run_group, work_path):
    '''
    Check the locations files for col/row ranges, bad UTM zones and wrong run groups
    return a df of facid and src_id
    '''
    fn = os.path.join(work_path, 'locations', '%s_locations.csv' %run_group)
    locs = pd.read_csv(fn, skipinitialspace=True)
    cols = list(locs['col'].drop_duplicates().sort_values())
    rows = list(locs['row'].drop_duplicates().sort_values())
    print('\nLow Col: %s  High Col: %s  Low Row: %s  High Row: %s' %(cols[0], cols[-1], rows[0], rows[-1]))
    # Check that all sources within a met fac are within the same UTM zone
    multi_zone = locs[['facid','src_id','utm_zone']].drop_duplicates(['facid','utm_zone'])
    print('\nFacilities with multiple UTM zones:')
    print(multi_zone[multi_zone.duplicated('facid')])
    # Check UTM locations
    utmd = locs[['facid','src_id','utmx','utmy','utm_zone','lat','lon']].drop_duplicates()
    utmp = UTM()
    utmd['zone_check'] = utmd['lon'].apply(lambda lon: utmp.get_zone(lon))
    print('\nNon matching UTM zones:')
    print(utmd[utmd['utm_zone'].astype(int) != utmd['zone_check'].astype(int)])
    utmd[['check_x','check_y']] = utmd[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utmp.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    print('\nNon matching UTM coords:')
    print(utmd[(utmd['utmx'].astype('f').round(0) != utmd['check_x'].astype('f').round(0)) \
      | (utmd['utmy'].astype('f').round(0) != utmd['check_y'].astype('f').round(0))])
    return locs[['facid','src_id']].copy()

def check_params(run_group, work_path):
    '''
    Check the source parameters for correct release params and for fac alignment/res
    '''
    params = pd.DataFrame()
    usecols = ['facid','src_id','rel_ht','sz']
    fn = os.path.join(work_path, 'parameters', '%s_area_params.csv' %run_group)
    params = pd.read_csv(fn, skipinitialspace=True, usecols=usecols)
    params.loc[params['src_id'].str.endswith('c1'), 'rel_ht_rg'] = 8.4
    params.loc[params['src_id'].str.endswith('c1'), 'sz_rg'] = 3.907
    params.loc[params['src_id'].str.endswith('c3'), 'rel_ht_rg'] = 20
    params.loc[params['src_id'].str.endswith('c3'), 'sz_rg'] = 40.7
    print('\nnon-matching rh:')
    print(params[params['rel_ht_rg'] != params['rel_ht']])
    print('\nnon-matching sz:')
    print(params[params['sz_rg'] != params['sz']])
    return params[['facid','src_id']].drop_duplicates()

def check_temporal(run_group, work_path):
    # QA the temporal file
    temp = os.path.join(work_path, 'temporal', '%s_temporal.csv' %run_group)
    temp = pd.read_csv(temp, skipinitialspace=True, index_col=False)
    scalars = [x for x in temp.columns if x.startswith('Scalar')]
    nomonth = temp[temp['qflag'] != 'MONTH'].copy()
    if len(nomonth) > 0:
        print('\nTemporal facilities without monthly profiles')
        print(nomonth)
    temp.ix[temp['qflag'] == 'MONTH', 'sum'] = temp.ix[temp['qflag'] == 'MONTH', scalars[:12]].sum(axis=1)
    temp['sum'] = temp['sum'].round(4)
    temp.to_csv(os.path.join(work_path, 'qa', '%s_temporal_check.csv' %run_group), index=False, 
      columns=['facid','src_id','qflag','sum'])
    temp = temp[['facid','src_id']].copy()
    return temp

def check_emis(run_group, work_path):
    # Read in the point emissions helper file
    emis = os.path.join(work_path, 'emis', '%s_emis.csv' %(run_group))
    emis = pd.read_csv(emis, dtype={'region_cd': str})
    emis['region_cd'] = emis['facid'].str.split('F').str[1].str[:5]
    emis.drop('state', axis=1, inplace=True)
    return emis 

def emis_compare(inv_rg, emis, run_group, work_path):
    '''
    Compare the inventory and AERMOD helper emissions
    '''
    inv_poll = inv_rg[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
    emis = emis.groupby(['region_cd','source_group','smoke_name'], as_index=False).sum()
    emis_poll = emis[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
    inv_poll = pd.merge(inv_poll, emis_poll, on='smoke_name', how='outer', suffixes=['_inv','_aer'])
    qa_out = os.path.join(work_path, 'qa', '%s_poll_qa.csv' %run_group)
    inv_poll.to_csv(qa_out, index=False)
    inv_sccs = inv_rg.groupby(['source_group','scc'], as_index=False).sum()
    qa_out = os.path.join(work_path, 'qa', '%s_scc_qa.csv' %run_group)
    inv_sccs.to_csv(qa_out, index=False)
    inv_rg = inv_rg.groupby(['region_cd','source_group'], as_index=False).sum()
    emis = emis.groupby(['region_cd','source_group'], as_index=False).sum()
    inv_rg = pd.merge(inv_rg, emis, on=['region_cd','source_group'], how='outer', 
      suffixes=['_inv','_aer'])
    qa_out = os.path.join(work_path, 'qa', '%s_emis_qa.csv' %run_group)
    inv_rg.to_csv(qa_out, index=False)

def get_invtable(fn):
    # Read in the inventory table to get the pollutants used in the AERMOD run
    invtable = pd.read_fwf(fn, comment='#', colspecs=[(0,11), (16,32), (41,42), (43,49)], 
      names=['smoke_name','poll','keep','spec_factor'], 
      converters={'smoke_name': str.strip, 'poll': str.strip, 'keep': str.strip})
    invtable = invtable[invtable['keep'].str.upper() == 'Y'].copy()
    invtable.drop_duplicates(inplace=True)
    return invtable

def get_inv(inv_list, invtable):
    # Read in all of the inventories from the inventory list
    # Select only the sources that have the selected HAPs and are not airports
    polls = list(invtable['poll'].drop_duplicates())
    inv = pd.DataFrame()
    usecols = ['poll','ann_value','region_cd','scc'] 
    dtype = {'poll': str, 'region_cd': str, 'scc': str}
    for fn in inv_list:
        with open(fn) as f:
            cmt_cnt = 0
            for l in f:
                if l.startswith('#'):
                    cmt_cnt = cmt_cnt + 1
                else:
                    break
        df = pd.read_csv(fn, skiprows=cmt_cnt, usecols=usecols, dtype=dtype) 
        df = df[(df['poll'].isin(polls)) & (df['ann_value'] > 0)].copy()
        df = pd.merge(df, invtable, on='poll', how='left')
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        inv = pd.concat((inv, df[['region_cd','scc','smoke_name','ann_value']]))
    inv = inv.groupby(['region_cd','scc','smoke_name'], as_index=False).sum()
    return inv[inv['ann_value'] > 0].copy()

def main():
    work_path = os.environ['WORK_PATH'] 
    invtable = os.environ['INVTABLE']
    source_groups = os.environ['SOURCE_GROUPS']
    if not os.path.exists(os.path.join(work_path, 'qa')):
        os.makedirs(os.path.join(work_path, 'qa'))
    # Read in the inventory table to get the pollutants used in the AERMOD run
    invtable = get_invtable(invtable)
    inv_list = get_inv_list()
    inv = get_inv(inv_list, invtable)
    src_groups = pd.read_csv(source_groups, usecols=['source_group','run_group','scc'], 
      dtype={'scc': str})
    inv = pd.merge(inv, src_groups, on='scc', how='left')
    run_groups = list(inv['run_group'].drop_duplicates())
    for run_group in run_groups:
        print('------ BEGIN %s ------' %run_group)
        inv_rg = inv[inv['run_group'] == run_group].copy()
        locs = check_locs(run_group, work_path)
        # Sum up locations by facility and source
        fac_locs, src_locs = (len(locs[['facid','src_id']].drop_duplicates('facid')), len(locs))
        print('\n%s Grid fac locations: %s    Source locations: %s' %(run_group, fac_locs, src_locs))
        qa_counts = [str(fac_locs), str(src_locs)]
        hdr = ['%s_loc_facs' %run_group, '%s_loc_src' %run_group]
        params = check_params(run_group, work_path)
        fac_params = len(params['facid'].drop_duplicates())
        src_params = len(params)
        print('%s Total unique param facs: %s   Total unique param sources: %s' %(run_group, fac_params, src_params))
        qa_counts += [str(fac_params), str(src_params)]
        hdr += ['%s_params_facs' %run_group, '%s_params_src' %run_group]
        temp = check_temporal(run_group, work_path)
        fac_temp, src_temp = (len(temp[['facid','src_id']].drop_duplicates('facid')), len(temp))
        print('%s Grid fac temporal: %s    Source temporal: %s' %(run_group, fac_locs, src_locs))
        qa_counts += [str(fac_temp), str(src_temp)]
        hdr += ['%s_temp_facs' %run_group, '%s_temp_src' %run_group]
        emis = check_emis(run_group, work_path)
        fac_emis, src_emis = (len(emis.drop_duplicates('facid')), 
            len(emis.drop_duplicates(['facid','src_id'])))
        print('%s Emis facs: %s   Emis sources: %s' %(run_group, fac_emis, src_emis))
        qa_counts += [str(fac_emis), str(src_emis)]
        hdr += ['%s_emis_facs' %run_group, '%s_emis_src' %run_group]
        # Write the facility and source counts to the "counts"  QA file
        count_out = os.path.join(work_path, 'qa', '%s_counts.csv' %run_group)
        with open(count_out, 'w') as f:
            f.write('%s\n%s\n' %(','.join(hdr), ','.join(qa_counts)))
        # AERMOD Source crosscheck - verify that all AERMOD sources are in all helper files
        xcheck = pd.DataFrame()
        for key, df in {'locs': locs, 'params': params, 'temp': temp, 'emis': emis}.items():
            df = df[['facid','src_id']].drop_duplicates()
            df[key] = 'Y'
            if xcheck.empty:
                xcheck = df.copy()
            else:
                xcheck = pd.merge(xcheck, df, on=['facid','src_id'], how='outer')
        xcheck = xcheck.fillna('N')
        qa_out = os.path.join(work_path, 'qa', '%s_src_id_qa.csv' %run_group)
        xcheck.to_csv(qa_out, index=False)
        # Do the emissions comparisons
        emis_compare(inv_rg, emis, run_group, work_path)
        print('------ END %s ------' %run_group)

if __name__ == '__main__':
    main()













