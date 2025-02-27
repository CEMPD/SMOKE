#!/usr/bin/env python
'''
QA for AERMOD area helper files
Requires python 3.5+, pandas, and pyproj
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
    return a df of met_cell and src_id
    '''
    fn = os.path.join(work_path, 'locations', '%s_locations.csv' %run_group)
    locs = pd.read_csv(fn, skipinitialspace=True)
    print('\nWrong run group:')
    print(locs[locs['run_group'] != run_group])
    # Look at top and bottom col/rows
    locs['cols'] = locs['met_cell'].str[1:4]
    cols = list(locs['cols'].drop_duplicates().sort_values())
    locs['rows'] = locs['met_cell'].str[5:]
    rows = list(locs['rows'].drop_duplicates().sort_values())
    print('\nLow Col: %s  High Col: %s  Low Row: %s  High Row: %s' %(cols[0], cols[-1], rows[0], rows[-1]))
    # Check that all sources within a met cell are within the same UTM zone
    multi_zone = locs[['met_cell','src_id','utm_zone']].drop_duplicates(['met_cell','utm_zone'])
    print('\nCells with multiple UTM zones:')
    print(multi_zone[multi_zone.duplicated('met_cell')])
    # Check UTM locations
    utmd = locs[['met_cell','src_id','utmx','utmy','utm_zone','lat','lon']].drop_duplicates()
    utmp = UTM()
    utmd['zone_check'] = utmd['lon'].apply(lambda lon: utmp.get_zone(lon))
    print('\nNon matching UTM zones:')
    print(utmd[utmd['utm_zone'].astype(int) != utmd['zone_check'].astype(int)])
    utmd[['check_x','check_y']] = utmd[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utmp.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    print('\nNon matching UTM coords:')
    print(utmd[(utmd['utmx'].astype('f').round(0) != utmd['check_x'].astype('f').round(0)) \
      | (utmd['utmy'].astype('f').round(0) != utmd['check_y'].astype('f').round(0))])
    return locs[['met_cell','src_id','utmx','utmy']].copy()

def get_group_params(fn, grid):
    '''
    Get the rungroup parameters
    '''
    df = pd.read_csv(fn)
    df.rename(columns={'release_height': 'rel_ht', 'sigma_z': 'sz'}, inplace=True)
    if grid != '12US2':
        # Fix the run group names based on the domain. Should be able to remove later
        df.ix[df['run_group'].str[-1] == '4', 'run_group'] = \
          df.ix[df['run_group'].str[-1] == '4', 'run_group'].str[:-1] + grid[:-1]
        df.ix[df['run_group'].str[-2:] == '12', 'run_group'] = \
          df.ix[df['run_group'].str[-2:] == '12', 'run_group'].str[:-2] + grid[:-1]
    return df.drop_duplicates('run_group')

def check_params(run_group, work_path, group_params):
    '''
    Check the source parameters for correct release params and for cell alignment/res
    '''
    params = pd.DataFrame()
    fn = os.path.join(work_path, 'parameters', '%s_area_params.csv' %run_group)
    params = pd.read_csv(fn, skipinitialspace=True)
    print('\nWrong run group:')
    print(params[params['run_group'] != run_group])
    print('\nWrong verts:')
    print(params[params['verts'] != 4])
    # Check release parameters by run group
    params = pd.merge(params, group_params, on='run_group', how='left', suffixes=['','_rg'])
    print('\nnon-matching rh:')
    print(params[params['rel_ht_rg'] != params['rel_ht']])
    print('\nnon-matching sz:')
    print(params[params['sz_rg'] != params['sz']])
    params['cols'] = params['met_cell'].str[1:4].astype(int)
    params['rows'] = params['met_cell'].str[5:].astype(int)
    params['res'] = params['src_id'].str.split('_').str.get(0)
    params['src'] = params['src_id'].str.split('_').str.get(1)
    subset_cells = list(params['met_cell'])[:200]
    # Check cell alignment alignment when a source cell size == met cell size
    # Check the alignment on a subset of met cells
    align = params.loc[params['met_cell'].isin(subset_cells), 
      ['cols','rows','src','res','utmx','utmy','x4','y4']].copy()
    met_subset = align.loc[align['src'] == '1', ['cols','rows','utmx','utmy']].copy()
    subset_check = align.loc[((align['res'] != '4') & (align['src'] == '1')) | \
      ((align['res'] == '4') & (align['src'] == '3')), ['cols','rows','x4','y4']].copy() 
    subset_check['cols'] = subset_check['cols'] + 1
    met_subset = pd.merge(met_subset, subset_check, how='left', on=['cols','rows'])
    print('\nUnaligned met cells from subset:')
    print(met_subset[(met_subset['x4'].notnull()) & ((met_subset['utmx'].round(0) != met_subset['x4'].round(0)) | \
      (met_subset['utmy'].round(0) != met_subset['y4'].round(0)))])
    if len(params[params['res'] == '4']) > 1:
        # Check alignment and width of a subset of sources
        src_subset = params.loc[(params['src'].isin(('4','5'))) & (params['met_cell'].isin(subset_cells)),
          ['met_cell','src','utmx','utmy','x4','y4']].copy()
        src_subset = pd.merge(src_subset.loc[src_subset['src'] == '4', ['met_cell','x4','y4']],
          src_subset.loc[src_subset['src'] == '5', ['met_cell','utmx','utmy']], on='met_cell', how='left')
        print('Unaligned source cells from subset:')
        print(src_subset[((src_subset['x4'].notnull()) & (src_subset['utmx'].notnull())) & ((src_subset['utmx'].round(0) != src_subset['x4'].round(0)) | \
          (src_subset['utmy'].round(0) != src_subset['y4'].round(0)))])
    # Check the resolution of the source cells. This is approximate because UTM m != Lambert m
    params['vres'] = (params['y3'] - params['y4'])/1000.
    params['hres'] = (params['x3'] - params['x2'])/1000.
    print('\nIncorrect resolution source count:')
    bad_res = params[(params['vres'].round(0).astype(int) != params['res'].astype(int)) | \
      (params['hres'].round(0).astype(int) != params['res'].astype(int))].copy()
    print(len(bad_res))
    qa_out = os.path.join(work_path, 'qa', '%s_badres_qa.csv' %run_group)
    bad_res.to_csv(qa_out, index=False)
    return params[['met_cell','src_id','utmx','utmy']].drop_duplicates()

def check_temp_hr(run_group, work_path):
    '''
    Read the hourly temporal files and check the number of counties and sums
    '''
    # Read in and do a count of temporal files
    hr = pd.DataFrame()
    from glob import glob
    temp_glob = '%s/temporal/%s_[0-9]*_hourly.csv' %(work_path, run_group)
    temp = pd.DataFrame()
    if run_group[:2] in ('AG','RW'):
        fips = 'region_cd'
        usecols = [fips,'factor']
    else:
        fips = 'fips'
        usecols = [fips,'scalar']
    for fn in glob(temp_glob):
        df = pd.read_csv(fn, dtype={fips: str}, usecols=usecols, 
          skipinitialspace=True)
        if fips == 'region_cd':
            df.rename(columns={'region_cd': 'fips'}, inplace=True)
        df = df.groupby('fips', as_index=False).agg(['sum','count'])
        df.reset_index(inplace=True)
        df.columns = ['fips','sum','count']
        df['sum'] = df['sum'].round(6)
        temp = pd.concat((temp, df))
    temp.to_csv(os.path.join(work_path, 'qa', 'temporal_hour_%s_check.csv' %run_group), index=False)
    return temp[['fips','count']] 

def check_temporal(run_group, work_path):
    '''
    Read in the non-hourly temporal files and check the number of cells and sums
    '''
    temp = os.path.join(work_path, 'temporal', '%s_temporal.csv' %run_group)
    temp = pd.read_csv(temp, skipinitialspace=True, index_col=False)
    scalars = [x for x in temp.columns if x.startswith('Scalar')]
    temp.ix[temp['qflag'] == 'MONTH', 'sum'] = temp.ix[temp['qflag'] == 'MONTH', scalars[:12]].sum(axis=1)
    temp.ix[temp['qflag'] == 'HROFDY', 'sum'] = temp.ix[temp['qflag'] == 'HROFDY', scalars[:24]].sum(axis=1)
    temp.ix[temp['qflag'] == 'MHRDOW7', 'sum'] = temp.ix[temp['qflag'] == 'MHRDOW7', scalars].sum(axis=1) * \
      (8760./2016.)
    if 'MHRDOW' in list(temp['qflag']):
        wkday = [scalars[x] for x in range(288)]
        sat = [scalars[x] for x in range(288,576)]
        sun = [scalars[x] for x in range(576,864)]
        temp.ix[temp['qflag'] == 'MHRDOW', 'sum'] = (5 * temp.ix[temp['qflag'] == 'MHRDOW', wkday].sum(axis=1) + \
          temp.ix[temp['qflag'] == 'MHRDOW', sat].sum(axis=1) + temp.ix[temp['qflag'] == 'MHRDOW', sun].sum(axis=1)) \
          * (8760./2016.)
    temp['sum'] = temp['sum'].round(4)
    temp.to_csv(os.path.join(work_path, 'qa', '%s_temporal_check.csv' %run_group), index=False, 
      columns=['met_cell','src_id','qflag','sum'])
    temp = temp[['met_cell','src_id']].copy()
    return temp

def check_emis(run_group, work_path, cell_res):
    '''
    Read in the area emissions helper file for comparison to inventory values
    '''
    emis = os.path.join(work_path, 'emis', '%s_%s_emis.csv' %(cell_res, run_group))
    emis = pd.read_csv(emis, dtype={'region_cd': str})
    if 'january' in emis.columns:
        months = ['january','february','march','april','may','june','july','august','september',
          'october','november','december']
        ann = emis[['region_cd','source_group','ann_value'] + months].groupby(['region_cd','source_group'], 
          as_index=False).sum()
        emis = emis[['region_cd','met_cell','src_id','source_group','smoke_name','ann_value']].copy()
        ann['ann_check'] = ann[months].sum(axis=1)
        ann['pdiff'] = abs(ann['ann_value'] - ann['ann_check']) / ann['ann_value']
        print('\nAnnual value mismatch:')
        print(ann[ann['pdiff'] > 0.001])
    return emis 

def get_invtable(fn):
    '''
    Read in the inventory table used for the AERMOD run
    '''
    invtable = pd.read_fwf(fn, comment='#', colspecs=[(0,11), (16,32), (41,42), (43,49)], 
      names=['smoke_name','poll','keep','spec_factor'], 
      converters={'smoke_name': str.strip, 'poll': str.strip, 'keep': str.strip})
    invtable = invtable[invtable['keep'].str.upper() == 'Y'].copy()
    invtable.drop_duplicates(inplace=True)
    return invtable

def get_inv(inv_list, invtable, grid):
    '''
    Read in the inventories used for the AERMOD run
    '''
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
        df.loc[df['poll'].str[3:5] == '__', 'poll'] = df.loc[df['poll'].str[3:5] == '__', 'poll'].str[5:]
        df = df[(df['poll'].isin(polls)) & (df['ann_value'] > 0)].copy()
        if grid == '3HI1':
            df = df[df['region_cd'].str[:2] == '15'].copy()
        elif grid == '9AK2':
            df = df[df['region_cd'].str[:2] == '02'].copy()
        elif grid == '3PR1':
            df = df[df['region_cd'].str[:2].isin(('72','78'))].copy()
        else:
            print('\nNOTE: Selecting CONUS sources only. Dropping Tribal.')
            df = df[~ df['region_cd'].str[:2].isin(('02','15','72','78','88'))].copy()
        df = pd.merge(df, invtable, on='poll', how='left')
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        inv = pd.concat((inv, df[['region_cd','scc','smoke_name','ann_value']]))
    inv = inv.groupby(['region_cd','scc','smoke_name'], as_index=False).sum()
    return inv[inv['ann_value'] > 0].copy()

def check_utm(locs, params):
    '''
    Check UTM alignment between locs and params
    '''
    utm_check = pd.merge(locs, params, on=['met_cell','src_id'], how='outer', suffixes=['_l','_p'])
    utm_check = utm_check[(utm_check['utmx_l'] != utm_check['utmx_p']) \
      | (utm_check['utmy_l'] != utm_check['utmy_p'])].copy()
    if len(utm_check) > 0:
        print('\nUnaligned Location and Parameter UTM SW corners:')
        print(utm_check)

def check_hourly_temp(run_group, inv_rg, work_path):
    '''
    Hourly county temporal helper files
    '''
    xwalk = os.path.join(work_path, 'xwalk', '%s_county-to-gridcell.csv' %run_group)
    xwalk = pd.read_csv(xwalk, dtype={'region_cd': str}, usecols=['region_cd','met_cell','src_id'])
    temp = check_temp_hr(run_group, work_path)
    temp.rename(columns={'fips': 'region_cd'}, inplace=True)
    temp = pd.merge(xwalk, temp, on='region_cd', how='left')
    temp_check = temp[temp['count'].isnull()].copy()
    if len(temp_check) > 0:
        print('\nUnmatched temporal sources:')
        print(temp_check)
    '''
    Check region_cd for run groups vs inven
    '''
    fips_check = pd.merge(inv_rg[['region_cd','run_group']].drop_duplicates(),
      xwalk[['region_cd','src_id']].drop_duplicates('region_cd'), on='region_cd', how='left')
    fips_check = fips_check[(fips_check['run_group'].isnull()) | (fips_check['src_id'].isnull())].copy()
    if len(fips_check) > 0:
        print('\nUnmatched FIPS x-walk vs. inven (these might be fine depending on alignment):')
        print(fips_check)
    return temp

def fix_run_groups(source_groups, grid):
    '''
    Fix the run groups names based on the grid name
    Hopefully this step will be removed in the future as the rungroup names are improved
    '''
    source_groups.ix[source_groups['run_group'].str[-1] == '4', 'run_group'] = \
      source_groups.ix[source_groups['run_group'].str[-1] == '4', 'run_group'].str[:-1] + grid[:-1]
    source_groups.ix[source_groups['run_group'].str[-2:] == '12', 'run_group'] = \
      source_groups.ix[source_groups['run_group'].str[-2:] == '12', 'run_group'].str[:-2] + grid[:-1]
    return source_groups 

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

def main():
    grid = os.environ['REGION_IOAPI_GRIDNAME']
    work_path = os.environ['WORK_PATH'] 
    invtable = os.environ['INVTABLE']
    group_params = os.environ['GROUP_PARAMS']
    source_groups = os.environ['SOURCE_GROUPS']
    if not os.path.exists(os.path.join(work_path, 'qa')):
        os.makedirs(os.path.join(work_path, 'qa'))
    # Read in the inventory table to get the pollutants used in the AERMOD run
    invtable = get_invtable(invtable)
    inv_list = get_inv_list()
    inv = get_inv(inv_list, invtable, grid)
    src_groups = pd.read_csv(source_groups, usecols=['source_group','run_group','scc'], 
      dtype={'scc': str})
    # Get the resolution of the domain
    if grid.startswith('12'):
        domain_res = '12'
    else:
        domain_res = grid[0]
        src_groups = fix_run_groups(src_groups, grid)
    inv = pd.merge(inv, src_groups, on='scc', how='left')
    run_groups = list(inv['run_group'].drop_duplicates())
    group_params = get_group_params(group_params, grid)
    for run_group in run_groups:
        print('------ BEGIN %s ------' %run_group)
        inv_rg = inv[inv['run_group'] == run_group].copy()
        if domain_res == '12' and run_group.endswith('4'):
            cell_res = '4'
        else:
            cell_res = domain_res
        locs = check_locs(run_group, work_path)
        # Sum up locations by facility and source
        cell_locs, src_locs = (len(locs[['met_cell','src_id']].drop_duplicates('met_cell')), len(locs))
        print('\n%s Grid cell locations: %s    Source locations: %s' %(run_group, cell_locs, src_locs))
        qa_counts = [str(cell_locs), str(src_locs)]
        hdr = ['%s_loc_cells' %run_group, '%s_loc_src' %run_group]
        params = check_params(run_group, work_path, group_params)
        check_utm(locs, params)
        cell_params = len(params['met_cell'].drop_duplicates())
        src_params = len(params)
        print('%s Total unique param cells: %s   Total unique param sources: %s' %(run_group, cell_params, src_params))
        qa_counts += [str(cell_params), str(src_params)]
        hdr += ['%s_params_cells' %run_group, '%s_params_src' %run_group]
        # Temporal checks for rungroups that use hourly/county or non-hourly
        if run_group[:3] in ('RWC','LDO','HDO','HOT','AG1'):
            temp = check_hourly_temp(run_group, inv_rg[['region_cd','run_group']], work_path)
        else:
            temp = check_temporal(run_group, work_path)
        cell_temp, src_temp = (len(temp[['met_cell','src_id']].drop_duplicates('met_cell')), len(temp))
        print('%s Grid cell temporal: %s    Source temporal: %s' %(run_group, cell_locs, src_locs))
        qa_counts += [str(cell_temp), str(src_temp)]
        hdr += ['%s_temp_cells' %run_group, '%s_temp_src' %run_group]
        emis = check_emis(run_group, work_path, cell_res)
        fac_emis, src_emis = (len(emis.drop_duplicates('met_cell')), 
            len(emis.drop_duplicates(['met_cell','src_id'])))
        print('%s Emis cells: %s   Emis sources: %s' %(run_group, fac_emis, src_emis))
        qa_counts += [str(fac_emis), str(src_emis)]
        hdr += ['%s_emis_cells' %run_group, '%s_emis_src' %run_group]
        # Write the facility and source counts to the "counts"  QA file
        count_out = os.path.join(work_path, 'qa', '%s_counts.csv' %run_group)
        with open(count_out, 'w') as f:
            f.write('%s\n%s\n' %(','.join(hdr), ','.join(qa_counts)))
        # AERMOD Source crosscheck - verify that all AERMOD sources are in all helper files
        xcheck = pd.DataFrame()
        for key, df in {'locs': locs, 'params': params, 'temp': temp, 'emis': emis}.items():
            df = df[['met_cell','src_id']].drop_duplicates()
            df[key] = 'Y'
            if xcheck.empty:
                xcheck = df.copy()
            else:
                xcheck = pd.merge(xcheck, df, on=['met_cell','src_id'], how='outer')
        xcheck = xcheck.fillna('N')
        qa_out = os.path.join(work_path, 'qa', '%s_srcid_qa.csv' %run_group)
        xcheck.to_csv(qa_out, index=False)
        # Do the inventory vs helper file emissions comparisons
        emis_compare(inv_rg, emis, run_group, work_path)
        print('------ END %s ------' %run_group)


if __name__ == '__main__':
    main()













