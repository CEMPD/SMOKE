#!/usr/bin/env python

import os, sys, string
import pandas as pd
import smk2ae
from smk2ae.hem_groups import HEMGroups
from smk2ae.utm import UTM
from smk2ae.ar_ff10 import AnnualFF10
from smk2ae.temporal import Temporal, match_temporal, fill_default

state_dict = {'02': 'AK', '15': 'HI', '72': 'PR', '78': 'VI'}

def read_shapes(shapes_file, fips_list):
    '''
    Read in the shapes from the shapes file

    The shapes file is formatted like this:
    facid,lon,lat,numvert,area
    '''
    df = pd.read_csv(shapes_file, usecols=['facid','lon','lat','numvert','area'])
    df['region_cd'] = df['facid'].str[1:6].str.zfill(6)
    df['tract'] = df['facid'].str[1:12]
    df['polynumber'] = df['facid'].str[-4:]
    df['src_id'] = df['tract'].str[-10:] + df['polynumber'].str[-2:]
    df.rename(columns={'area': 'area_poly'}, inplace=True)
    tract_area = df[['tract','area_poly']].drop_duplicates().groupby('tract', as_index=False).sum()
    tract_area.rename(columns={'area_poly': 'area_tract'}, inplace=True)
    df = pd.merge(df, tract_area, on='tract', how='left')
    return df

class nonconus_surgs:
    def __init__(self, surg_path, surg_prefix):
        self.path = surg_path
        self.prefix = surg_prefix
        self.table = pd.DataFrame() 
        self._loaded_surgs = []

    def read_xref(self, srg_xref):
        '''
        Read in the nonconus surrogate x-ref
        '''
        self.xref = pd.read_csv(srg_xref, dtype={'fips': '|S6', 'scc': '|S10', 'profcode': '|S14'})
        self.xref['region_cd'] = self.xref['fips'].str.zfill(6)
        self.xref['scc'] = self.xref['scc'].str.zfill(10)
        self.xref.drop('fips', axis=1, inplace=True)
        self.xref.fillna(0, inplace=True)
        self.xref.ix[self.xref['scc'].str.endswith('000'), 'scc7'] = \
          self.xref.ix[self.xref['scc'].str.endswith('000'), 'scc'].str[:7]
        self.xref.ix[self.xref['region_cd'].str.endswith('000'), 'cost'] = \
          self.xref.ix[self.xref['region_cd'].str.endswith('000'), 'region_cd'].str[:3]
        duplicates = self.xref[self.xref.duplicated(subset=['region_cd','scc'], keep=False)]
        if not duplicates.empty:
            print duplicates
            raise ValueError, 'Duplicate surrogate cross-reference'

    def load_surgs(self, fips_list=False, prof_codes=False):
        '''
        Load the non CONUS gridding surrogates
        profcode, fips, tract, pct
        '''
        if not prof_codes:
            prof_codes = list(self.xref['profcode'].drop_duplicates())
        for surg_num in prof_codes:
            if surg_num not in self._loaded_surgs:
                fname = os.path.join(self.path, '%s_%s.csv' %(self.prefix, surg_num))
                df = pd.read_csv(fname, dtype={'profcode': '|S14', 'fips': '|S6', 'tract': '|S12'})
                df = df[df['pct'] > 0].copy()
                df['region_cd'] = df['fips'].str.zfill(6)
                if fips_list:
                    df = df[df['region_cd'].isin(fips_list)].copy()
                # Renormalize the surrogates by FIPS
                fips_sum = df[['fips','pct']].groupby('fips', as_index=False).sum()
                df = pd.merge(df, fips_sum, on='fips', how='left', suffixes=['','_sum'])
                df['pct'] = df['pct'] / df['pct_sum']
                df.drop(['fips','pct_sum'], axis=1, inplace=True)
                self.table = pd.concat((self.table, df))
                self._loaded_surgs.append(surg_num)

    def match_emis(self, emis):
        '''
        Match the surrogates to the emissions
        '''
        hierarchy = [['region_cd','scc'], ['region_cd','scc7'], ['cost','scc'], ['cost','scc7'], 
         ['scc',], ['scc7',], ['region_cd',], ['cost',]]
        short_keys = {'scc': 'scc7', 'region_cd': 'cost'}
        emis['scc7'] = emis['scc'].str[:7]
        emis['cost'] = emis['region_cd'].str[:3]
        matched_df = pd.DataFrame()
        for key_cols in hierarchy:
            xref = self.xref.copy()
            zero_cols = [col for col in xref.columns if col not in key_cols and col not in \
              ['profcode',]+short_keys.values()]
            for zero_col in zero_cols:
                if short_keys[zero_col] not in key_cols:
                    xref = xref[xref[zero_col].astype('i') == 0].copy()
            xref = xref[key_cols+['profcode',]].drop_duplicates()
            df = pd.merge(emis, xref, on=key_cols, how='left')
            matched_df = pd.concat((matched_df, df[df['profcode'].notnull()]))
            emis = df[df['profcode'].isnull()].copy()
            if emis.empty:
                break
            else:
                emis.drop('profcode', axis=1, inplace=True)
        if not emis.empty:
            print emis.drop_duplicates(['region_cd','scc'])
            raise ValueError, 'Unmatched emissions to surrogates'
        self.load_surgs(list(matched_df['region_cd']), list(matched_df['profcode']))
        matched_df = pd.merge(matched_df, self.table, on=['profcode','region_cd'], how='left')
        matched_df.drop(['profcode','scc7','cost'], axis=1, inplace=True)
        matched_df.drop_duplicates(inplace=True)
        return matched_df

def write_locations(type_df, shapes, run_group, state):
    '''
    Write the locations file
    '''
    out_cols = ['run_group','state','region_cd','tract','polynumber','facid','src_id',
      'area_poly','area_tract','utmx','utmy','utm_zone','lon','lat']
    type_df = pd.merge(type_df, shapes, on=['region_cd','tract'], how='left')
#    type_df.drop_duplicates('facid', inplace=True)
    if not type_df[type_df['utmx'].isnull()].empty:
        print 'WARNING: Missing polygon tracts'
        print type_df[type_df['utmx'].isnull()].drop_duplicates('tract')
    type_df['state'] = state_dict[state] 
    fname = os.path.join(os.environ['WORK_PATH'], 'locations', 
      'tract_%s_%s_locations.csv' %(state, run_group))
    type_df.sort().to_csv(fname, columns=out_cols, index=False)

rel_parms = [['NPHI12',10.,4.7],['NPLO12',3.9,3.6],['NONRD12',2.,1.],['RWC12',6.4,3.2],
  ['NPHI4',10.,4.7],['NPLO4',3.9,3.6],['NONRD4',2.,1.],['RWC4',6.4,3.2],
  ['LDOFF12',1.3,1.2],['LDON12',1.3,1.2],['HDOFF12',3.4,3.2],['HDON12',3.4,3.2],
  ['LDOFF4',1.3,1.2],['LDON4',1.3,1.2],['HDOFF4',3.4,3.2],['HDON4',3.4,3.2],
  ['OILGAS4',10.,4.7],['HOTEL4',3.4,0.5]]
rp_df = pd.DataFrame(rel_parms, columns=['run_group','rel_ht','sz'])
def write_parameters(type_df, shapes, run_group, state):
    '''
    Write the parameters file
    
    The maximum number of vertices for a polygon needs to be calculated first
    The coordinates for each polygon is written to a single line
    '''
    param_df = pd.merge(type_df, shapes, on=['region_cd','tract'], how='left')
    param_df = param_df[param_df['facid'].notnull()].copy()
    fname = os.path.join(os.environ['WORK_PATH'], 'parameters', 'tract_%s_%s_area_params.csv' %(state, run_group))
    maxverts = int(param_df['numvert'].max())
    coord_cols = []
    [coord_cols.extend(['utmx','utmy']) for vert_num in xrange(1, maxverts + 1)]
    ll_cols = []
    [ll_cols.extend(['lon','lat']) for vert_num in xrange(1, maxverts + 1)]
    out_cols = ['run_group','state','region_cd','tract','polynumber','facid','src_id','type',
      'area_poly','rh','numvert','sz'] + coord_cols + ll_cols
    param_df = pd.merge(param_df, rp_df, on='run_group', how='left')  
    param_df['state'] = state_dict[state]
    with open(fname,'w') as f:
        f.write('%s\n' %','.join(out_cols))
        for facid in list(param_df['facid'].drop_duplicates()):
            src_df = param_df[param_df['facid']==facid].copy()
            numverts = int(src_df['numvert'].values[0])
            out_line = [run_group, state_dict[state], src_df['region_cd'].values[0], 
              src_df['tract'].values[0], src_df['polynumber'].values[0], facid, 
              src_df['src_id'].values[0], 'AREAPOLY', str(src_df['rel_ht'].values[0]),
              str(numverts), str(src_df['sz'].values[0])]
            for ix, row in src_df.iterrows():
                out_line.extend([str(row['utmx']),str(row['utmy'])])
            if numverts < maxverts:
                for x in xrange(maxverts - numverts):
                    out_line.extend(['',''])
            for ix, row in src_df.iterrows():
                out_line.extend([str(row['lon']),str(row['lat'])])
            if numverts < maxverts:
                for x in xrange(maxverts - numverts):
                    out_line.extend(['',''])
            f.write('%s\n' %','.join(out_line))

def write_daily_prof(df, temp, run_group, state):
    '''
    Write the daily profile based hourly profiles as a normalized list by year/month/day/hour
    Daily factors should be based on sum of hourly haps over sum of annual haps
    '''
    cols = ['region_cd','scc','run_group','facid','tract','polynumber','src_id','ann_value']
    hierarchy = [['region_cd','scc'],['region_cd',]]
    df = df[cols].copy().sort('ann_value').drop_duplicates(['run_group','region_cd','facid'], 
      take_last=True)
    temp.aermod.fillna(0, inplace=True)
    value_cols = ['month','day','hour','factor']
    df = df[['run_group','region_cd','scc',
      'ann_value']].copy().drop_duplicates(['run_group','region_cd','scc'])
    df = match_temporal(df, temp.aermod, value_cols, hierarchy, temp.use_daily)
    df.drop('scc', axis=1, inplace=True)
    df['year'] = os.environ['BASE_YEAR'][2:4]
    df.drop_duplicates(['region_cd','run_group','year','month','day','hour'], inplace=True)
    keys = ['region_cd','run_group','year']
    mean_factor = df[['region_cd','factor']].copy().groupby('region_cd', as_index=False).mean()
    df = pd.merge(df, mean_factor, on='region_cd', how='left', suffixes=['','_mean'])
    df['factor'] = (df['factor']/df['factor_mean']) 
    df[['month','day','hour']] = df[['month','day','hour']].astype('i')
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','tract_%s_%s_hourly.csv' %(state, 
      run_group))
    df.to_csv(fname, columns=['run_group','region_cd','year','month','day','hour','factor'],
      index=False)

# Default temporal profile
def_prof = pd.Series(['HROFDY',0.4751525,0.4463554,0.4367563,0.4487551,0.5039496,0.59994,0.7463254,
  0.9311069,1.120688,1.267073,1.370263,1.449455,1.487851,1.514249,1.523848,1.497450,1.425457,
  1.315068,1.274273,1.221478,1.019898,0.7847215,0.6167383,0.5231477], 
  index=['qflag',]+['Scalar%s' %x for x in xrange(1,25)])
def_oil_prof = pd.Series(['MONTH',1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0], 
  index=['qflag',]+['Scalar%s' %x for x in xrange(1,13)])

def write_no_daily_prof(run_df, temp, run_group, state):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    '''
    cols = ['region_cd','scc','run_group','facid','tract','polynumber','src_id','ann_value']
    run_df = run_df[cols].copy().drop_duplicates()
    # Use main contributing region_cd and SCC for source temporalization
    run_df.sort('ann_value', inplace=True)
    run_df.drop_duplicates('facid', keep='last', inplace=True)
    scalar_cols = [s_col for s_col in temp.aermod.columns if s_col.startswith('Scalar')]
    value_cols = ['qflag',] + scalar_cols
    # Only match by fips/scc or fips; SCC only gets default
    hierarchy = [['region_cd','scc'],['region_cd',]]
    run_df = match_temporal(run_df, temp.aermod, value_cols, hierarchy)
    run_df.drop(['scc','ann_value'], axis=1, inplace=True)
    run_df.drop_duplicates(inplace=True)
    if not run_df[run_df[def_prof.index[0]].isnull()].empty:
        if run_group == 'OILGAS4':
            run_df = fill_default(run_df, def_oil_prof)
        else:
            run_df = fill_default(run_df, def_prof)
    qflag_list = list(run_df['qflag'].drop_duplicates())
    if 'MHRDOW7' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in xrange(1,2017)]
    elif 'MHRDOW' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in xrange(1,865)]
    elif 'HROFDY' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in xrange(1,25)]
    else:
        scalar_cols = ['Scalar%s' %x for x in xrange(1,13)]
    cols = ['run_group','state','region_cd','tract','polynumber','facid','src_id','qflag']
    run_df['state'] = state_dict[state]
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','tract_%s_%s_temporal.csv' %(state,
      run_group))
    run_df.to_csv(fname, index=False, columns=cols+scalar_cols)
      
def write_aermod_emis(inv, surgs, shapes, hem, run_group, state):
    '''
    Write the "Post AERMOD" emissions files
    '''
    out_cols = ['run_group','state','region_cd','tract','polynumber','facid','src_id',
      'source_group','smoke_name','ann_value']
    src_groups = hem[['scc','source_group']].drop_duplicates()
    group_sccs = list(hem.ix[hem['run_group'] == run_group, 'scc'].drop_duplicates())
    shapes = shapes[['region_cd','tract','facid','polynumber','src_id',
      'area_poly','area_tract']].copy().drop_duplicates(['region_cd','facid'])
    emis = inv.inv_emis[(inv.inv_emis['scc'].isin(group_sccs)) & \
      (inv.inv_emis['region_cd'].str[1:3] == state)].copy()
    group_df = surgs.match_emis(emis)
    group_df = pd.merge(group_df, shapes, on=['region_cd','tract'], how='left')
    group_df['pct'] = group_df['pct'] * (group_df['area_poly'] / group_df['area_tract'])
    group_df['ann_value'] = group_df['ann_value'] * group_df['pct'] 
    if inv.seasons:
        seasons = ['winter','spring','summer','fall']
        out_cols += seasons
        for season in seasons:
            group_df[season] = group_df[season] * group_df['pct']
    group_df = pd.merge(group_df, src_groups, on='scc', how='left') 
    group_df.drop(['pct','scc','area_poly','area_tract'], axis=1, inplace=True)
    group_df = group_df.groupby(['region_cd','source_group','tract','polynumber','facid',
      'src_id','smoke_name'], as_index=False).sum()
    group_df['run_group'] = run_group
    group_df['state'] = state_dict[state]
    fname = os.path.join(os.environ['WORK_PATH'],'emis','%s_%s_emis.csv' %(state, run_group))
    group_df.to_csv(fname, columns=out_cols, index=False)

def check_env_vars(vars):
    '''
    Verify that all required environment variables are set
    '''
    for env_var in vars:
        try:
            os.environ[env_var]
        except KeyError:
            raise KeyError, 'Missing environment variable: %s' %env_var

def get_inv_list():
    '''
    Get a list of all the inventory files from the environment variables
    '''
    inv_list = [] 
    for inv_var in ['EMISINV_%s' %alpha for alpha in string.ascii_uppercase]:
        if inv_var in os.environ:
            inv_list.append(os.environ[inv_var])
    return inv_list

def main():
    var_list = ('EMISINV_A','HAPS_LIST','BASE_YEAR','WORK_PATH','POLY_FILE','STATE_FIPS',
      'HEM_GROUPS','SRG_XREF','SRG_PATH','SRG_PREFIX','RUN_GROUPS')
    check_env_vars(var_list)
    inv_list = get_inv_list()
    state_fips = os.environ['STATE_FIPS']
    inv = AnnualFF10(inv_list, os.environ['HAPS_LIST'], state_fips)
    hem = HEMGroups(os.environ['HEM_GROUPS'])
    inv.emis = pd.merge(inv.emis, hem.xref, on='scc', how='left')
    # Narrow down to just the selected run groups.
    if os.environ['RUN_GROUPS']:
        run_groups = [rg.strip() for rg in os.environ['RUN_GROUPS'].split(',')]
        inv.emis = inv.emis[inv.emis['run_group'].isin(run_groups)].copy()
        inv.update_lists()
    else:
        run_groups = list(inv.emis['run_group'].drop_duplicates())
    if inv.emis.empty:
        raise ValueError, 'No emissions found for given states and run groups.'
    surgs = nonconus_surgs(os.environ['SRG_PATH'], os.environ['SRG_PREFIX'])
    surgs.read_xref(os.environ['SRG_XREF'])
    inv.emis = surgs.match_emis(inv.emis)
    if not inv.emis[inv.emis['tract'].isnull()].empty:
        print inv.emis[inv.emis['tract'].isnull()].drop_duplicates()
        raise ValueError, 'Missing tract cross-reference'
    shapes = read_shapes(os.environ['POLY_FILE'], inv.fips_list)
    utm_zones = shapes[['facid','lon']].copy().drop_duplicates('facid')
    utm = UTM()
    utm_zones['utm_zone'] = utm_zones['lon'].apply(utm.get_zone)
    shapes = pd.merge(shapes, utm_zones[['facid','utm_zone']], on='facid', how='left')
    shapes[['utmx','utmy']] = shapes[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    temp = Temporal(os.environ['ATREF'],os.environ['ATPRO_HOURLY'],os.environ['ATPRO_WEEKLY'],
        os.environ['ATPRO_MONTHLY'], inv.scc_list, inv.fips_list)
    for run_group in list(inv.emis['run_group'].drop_duplicates()):
        for state in list(inv.emis['region_cd'].str[1:3].drop_duplicates()):
            write_aermod_emis(inv, surgs, shapes, hem.xref, run_group, state)
            emis_group = inv.emis.ix[(inv.emis['run_group'] == run_group) & \
             (inv.emis['region_cd'].str[1:3] == state), ['region_cd','run_group','tract']].copy()
            emis_group.drop_duplicates(inplace=True)
            write_locations(emis_group, shapes.copy().drop_duplicates('facid'), run_group, state)
            write_parameters(emis_group, shapes, run_group, state) 
            temp_group = pd.merge(inv.emis[(inv.emis['run_group'] == run_group) & \
              (inv.emis['region_cd'].str[1:3] == state)], shapes[['region_cd','facid','tract',
              'polynumber','src_id']].drop_duplicates(), on=['region_cd','tract'], how='left')
            if temp.use_daily:
                write_daily_prof(temp_group, temp, run_group, state)
            else:
                write_no_daily_prof(temp_group, temp, run_group, state)

if __name__ == '__main__':
	main()




