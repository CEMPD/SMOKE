import os
import os.path
import pandas as pd
from smk2ae.temporal import match_temporal, fill_default 
from smk2ae.utm import UTM

def proc_area(df, grid_info, temp):
    '''
    General function for processing area sources
      Iterates over each col and row to be written.
    '''
    utm = UTM()
    df = calc_cell_corner(df, grid_info, utm)
    write_locations(df)
    write_parameters(df, grid_info, utm)
    if temp.use_daily:
        write_daily_prof(df, temp)
    else:
        write_no_daily_prof(df, temp)

def calc_cell_corner(df, grid_info, utm):
    '''
    Calc SW UTM corner
    '''
    df['utm_zone'] = df['met_centroid_lon'].apply(utm.get_zone)
    df[['lon','lat']] = grid_info.colrow_to_ll(df['col'], df['row'])
    df[['utmx','utmy']] = df[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    return df

def calc_polygons(df, grid, utm):
    '''
    Calc remainder of UTM polygon
    '''
    df[['lon2','lat2']] = grid.colrow_to_ll(df['col'], df['row']+1)
    df[['lon3','lat3']] = grid.colrow_to_ll(df['col']+1, df['row']+1)
    df[['lon4','lat4']] = grid.colrow_to_ll(df['col']+1, df['row'])
    df[['x2','y2']] = df[['lon2','lat2','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon2'],x['lat2'],x['utm_zone'])), axis=1)
    df[['x3','y3']] = df[['lon3','lat3','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon3'],x['lat3'],x['utm_zone'])), axis=1)
    df[['x4','y4']] = df[['lon4','lat4','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon4'],x['lat4'],x['utm_zone'])), axis=1)
    return df

def write_locations(df):
    '''
    Write the location information for the grid cell and any sub-cells
    This is the "locations" file
    '''
    cols = ['run_group','met_cell','src_id','utmx','utmy','utm_zone','lon','lat']
    df = df[cols].copy().drop_duplicates()
    df['type'] = 'AREAPOLY'
    # Iterate over the source groupings for the area source
    for run_group in list(df['run_group'].drop_duplicates()):
        run_df = df[df['run_group'] == run_group].copy().sort('src_id')
        fname = os.path.join(os.environ['WORK_PATH'],'locations','%s_locations.csv' %run_group)
        run_df.to_csv(fname, columns=cols, index=False)

rel_parms = [['NPHI12',10.,4.7],['NPLO12',3.9,3.6],['NONRD12',2.,1.],['RWC12',6.4,3.2],
  ['NPHI4',10.,4.7],['NPLO4',3.9,3.6],['NONRD4',2.,1.],['RWC4',6.4,3.2],
  ['LDOFF12',1.3,1.2],['LDON12',1.3,1.2],['HDOFF12',3.4,3.2],['HDON12',3.4,3.2],
  ['LDOFF4',1.3,1.2],['LDON4',1.3,1.2],['HDOFF4',3.4,3.2],['HDON4',3.4,3.2],
  ['OILGAS4',10.,4.7],['HOTEL4',3.4,0.5]]
rp_df = pd.DataFrame(rel_parms, columns=['run_group','rel_ht','sz'])
def write_parameters(df, grid_info, utm):
    '''
    Write the release parameters for the grid cell and any sub-cells
    This is the "params" file
    '''
    cols = ['run_group','met_cell','src_id','utmx','utmy','col','row','utm_zone','lon','lat']
    df = df[cols].copy().drop_duplicates()
    df = calc_polygons(df, grid_info, utm)
    df = pd.merge(df, rp_df, on='run_group', how='left')
    df['verts'] = 4
    out_cols = ['run_group','met_cell','src_id','rel_ht','verts','sz','utmx','utmy','x2','y2',
      'x3','y3','x4','y4','lon','lat','lon2','lat2','lon3','lat3','lon4','lat4']
    for run_group in list(df['run_group'].drop_duplicates()):
        run_df = df[df['run_group'] == run_group].copy().sort('src_id')
        fname = os.path.join(os.environ['WORK_PATH'],'parameters','%s_area_params.csv' %run_group)
        run_df.to_csv(fname, columns=out_cols, index=False)

def write_daily_prof(df, temp):
    '''
    Write the daily profile based hourly profiles as a normalized list by year/month/day/hour
    Daily factors should be based on sum of hourly haps over sum of annual haps
    '''
    cols = ['region_cd','scc','run_group','met_cell','ann_value','src_id']
    hierarchy = [['region_cd','scc'],['region_cd',]]
    df = df[cols].copy().sort('ann_value').drop_duplicates(['run_group','region_cd','met_cell'], 
      take_last=True)
    temp.aermod.fillna(0, inplace=True)
    value_cols = ['month','day','hour','factor']
    for run_group in list(df['run_group'].drop_duplicates()):
        run_df = df[df['run_group'] == run_group].copy()
        fips_map = get_fips_map(run_df[['run_group','region_cd','met_cell','src_id',
          'ann_value']].copy())
        write_county_xwalk(fips_map, run_group)
        run_df = run_df[['run_group','region_cd','scc',
          'ann_value']].copy().drop_duplicates(['run_group','region_cd','scc'])
        for state in list(run_df['region_cd'].str[1:3].drop_duplicates()):
            st_df = run_df.ix[run_df['region_cd'].str[1:3] == state].copy() 
            st_df = match_temporal(st_df, temp.aermod, value_cols, hierarchy, temp.use_daily)
            st_df.drop('scc', axis=1, inplace=True)
            st_df['year'] = os.environ['BASE_YEAR'][2:4]
            st_df.drop_duplicates(['region_cd','run_group','year','month','day','hour'], 
              inplace=True)
            keys = ['region_cd','run_group','year']
#            sum_df = st_df[keys+['ann_value',]].groupby(keys, as_index=False).sum()
#            sum_df.rename(columns={'ann_value': 'sum'}, inplace=True)
#            src_df = pd.merge(st_df, sum_df, on=keys, how='left')
            mean_factor = st_df[['region_cd','factor']].copy().groupby('region_cd', as_index=False).mean()
            st_df = pd.merge(st_df, mean_factor, on='region_cd', how='left', suffixes=['','_mean'])
#            st_df['factor'] = st_df['factor'] * 8760.
            st_df['factor'] = (st_df['factor']/st_df['factor_mean']) 
            st_df[['month','day','hour']] = st_df[['month','day','hour']].astype('i')
            fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_%s_hourly.csv' %(run_group, state))
            st_df.to_csv(fname, columns=['run_group','region_cd','year','month','day','hour','factor'],
              index=False)

# Default temporal profile
def_prof = pd.Series(['HROFDY',0.4751525,0.4463554,0.4367563,0.4487551,0.5039496,0.59994,0.7463254,
  0.9311069,1.120688,1.267073,1.370263,1.449455,1.487851,1.514249,1.523848,1.497450,1.425457,
  1.315068,1.274273,1.221478,1.019898,0.7847215,0.6167383,0.5231477], 
  index=['qflag',]+['Scalar%s' %x for x in xrange(1,25)])
def get_no_daily_prof(run_df, temp):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    '''
    scalar_cols = [s_col for s_col in temp.aermod.columns if s_col.startswith('Scalar')]
    value_cols = ['qflag',] + scalar_cols
    # Only match by fips/scc or fips; SCC only gets default
    hierarchy = [['region_cd','scc'],['region_cd',]]
    run_df = match_temporal(run_df, temp.aermod, value_cols, hierarchy)
    run_df.drop(['region_cd','scc'], axis=1, inplace=True)
    run_df.drop_duplicates(inplace=True)
    if not run_df[run_df[def_prof.index[0]].isnull()].empty:
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
    return run_df[['run_group','met_cell','src_id','qflag']+scalar_cols].copy()

def write_no_daily_prof(df, temp):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    '''
    cols = ['region_cd','scc','run_group','met_cell','src_id','ann_value']
    df = df[cols].copy().drop_duplicates(['region_cd','run_group','met_cell','src_id'])
    for run_group in list(df['run_group'].drop_duplicates()):
        run_df = df[df['run_group'] == run_group].copy().sort('src_id')
        states = list(run_df['region_cd'].str[1:3].drop_duplicates())
        if run_group in ('LDOFF12','LDON12','HDOFF12','HDON12','LDOFF4','LDON4','HDON4','HDOFF4',
          'HOTEL4'):
            fips_map = get_fips_map(run_df[['run_group','region_cd','met_cell','src_id',
              'ann_value']].copy())
            write_county_xwalk(fips_map, run_group)
            for state in states:
                st_df = get_no_daily_prof(run_df[run_df['region_cd'].str[1:3] == state], temp)
                write_onroad_hourly(st_df, fips_map[fips_map['region_cd'].str[1:3] == state], 
                  run_group, state)
        else:
            fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_temporal.csv' %run_group)
            with open(fname, 'w') as temp_out:
                for idx, state in enumerate(states):
                    st_df = get_no_daily_prof(run_df[run_df['region_cd'].str[1:3] == state], temp)
                    if idx == 0:
                        st_df.to_csv(temp_out, index=False)
                    else:
                        st_df.to_csv(temp_out, index=False, header=False)

def get_fips_map(df):
    df.sort('ann_value', inplace=True)
    df.drop_duplicates(['run_group','met_cell','src_id'], keep='last', inplace=True)
    df.drop('ann_value', axis=1, inplace=True)
    return df

def write_county_xwalk(fips_map, run_group):
    xwalk = os.path.join(os.environ['WORK_PATH'],'xwalk','%s_county-to-gridcell.csv' %run_group)
    fips_map.to_csv(xwalk, index=False)

def write_onroad_hourly(df, fips_map, run_group, state):
    '''
    Write the onroad temporal format
    '''
    df = pd.merge(df, fips_map, on=['run_group','met_cell','src_id'], how='left')
    df.drop(['met_cell','src_id','qflag'], axis=1, inplace=True)
    df.drop_duplicates(['run_group','region_cd'], inplace=True)
    df = pd.melt(df, id_vars=['run_group','region_cd'], var_name='tstep', value_name='scalar')
    df['tstep'] = df['tstep'].str[6:].astype('f')
    from os import environ
    df['year'] = str(environ['BASE_YEAR'])[2:]
    df = pd.concat((df, vcalc_mdh(df['tstep'])), axis=1)
#    df[['month','dow','hour']] = df['tstep'].apply(calc_mdh)
    cols = ['run_group','region_cd','year','month','dow','hour','scalar']
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_%s_hourly.csv' %(run_group, state))
    df.to_csv(fname, columns=cols, index=False)

def calc_mdh(tstep):
    '''
    Calculate month, dow, and hour from tstep
    '''
    from math import ceil, floor
    tstep = float(tstep[6:]) - 1
    dow = ceil(tstep / (24. * 12.))
    dow_start = ((dow - 1) * (24. * 12.))
    mon = ceil((dow_start - dow_start) / 24.)
    mon_start = dow_start + ((mon - 1) * 24.)
    hour = tstep - mon_start + 1
    return pd.Series([mon,dow,hour])

def vcalc_mdh(tstep):
    '''
    Calculate month, dow, and hour from tstep
    '''
    from numpy import ceil, floor
    dow = ceil(tstep / (24. * 12.))
    dow_start = ((dow - 1) * (24. * 12.))
    mon = ceil((tstep - dow_start) / 24.)
    mon_start = dow_start + ((mon - 1) * 24.)
    hour = (tstep - mon_start) 
    df = pd.concat([mon.astype('i'), dow.astype('i'), hour.astype('i')], axis=1)
    df.columns = ['month','dow','hour']
    return df


