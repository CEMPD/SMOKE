import os
import os.path
import pandas as pd
from smk2ae.temporal import match_temporal, fill_default 
from smk2ae.utm import UTM

def proc_area(df, grid_info, temp, params):
    '''
    General function for processing area sources
      Iterates over each col and row to be written.
    '''
    utm = UTM()
    df = calc_cell_corner(df, grid_info, utm)
    write_locations(df)
    write_parameters(df, grid_info, utm, params)
    if temp:
        if temp.use_daily:
            write_daily_prof(df, temp)
        else:
            write_no_daily_prof(df, temp)
    else:
        '''
        See if I need to write the temporal county xrefs
        '''
        run_group = df['run_group'].values[0]
        if run_group[:3] in ('LDO','HDO','HOT'):
            fips_map = get_fips_map(df[['run_group','region_cd','met_cell','src_id',
              'ann_value']].copy())
            write_county_xwalk(fips_map, run_group)


def calc_cell_corner(df, grid_info, utm):
    '''
    Calc SW UTM corner. This is used to identify UTM zone for the cell
    '''
    df['utm_zone'] = df['met_swcorner_lon'].apply(utm.get_zone)
    df[['lon','lat']] = grid_info.colrow_to_ll(df['col'], df['row'])
    df[['utmx','utmy']] = df[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    return df

def calc_polygons(df, grid, utm):
    '''
    Calc remainder of UTM polygon (vertices 2-4) based on lat/lon and UTM zone
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
    run_group = df['run_group'].values[0]
    df.sort_values('src_id', inplace=True)
    fname = os.path.join(os.environ['WORK_PATH'],'locations','%s_locations.csv' %run_group)
    df.to_csv(fname, columns=cols, index=False, float_format='%.12g')

def write_parameters(df, grid_info, utm, params):
    '''
    Write the release parameters for the grid cell and any sub-cells
    This is the "params" file
    '''
    cols = ['run_group','met_cell','src_id','utmx','utmy','col','row','utm_zone','lon','lat']
    df = df[cols].copy().drop_duplicates()
    # The polygons are grid cells
    df = calc_polygons(df, grid_info, utm)
    df = pd.merge(df, params, on='run_group', how='left')
    if len(df[df['rel_ht'].isnull()]) > 0:
        print('Unmatched run groups in GROUP PARAMS file')
        print(df.loc[df['rel_ht'].isnull(), 'run_group'].drop_duplicates())
        raise ValueError('Unmatched run group parameters')
    df['verts'] = 4
    out_cols = ['run_group','met_cell','src_id','rel_ht','verts','sz','utmx','utmy','x2','y2',
      'x3','y3','x4','y4','lon','lat','lon2','lat2','lon3','lat3','lon4','lat4']
    run_group = df['run_group'].values[0]
    df.sort_values('src_id', inplace=True)
    fname = os.path.join(os.environ['WORK_PATH'],'parameters','%s_area_params.csv' %run_group)
    df.to_csv(fname, columns=out_cols, index=False, float_format='%.12g')

def write_daily_prof(df, temp):
    '''
    Write the daily profile based hourly profiles as a normalized list by year/month/day/hour
    Daily factors should be based on sum of hourly haps over sum of annual haps
    '''
    cols = ['region_cd','scc','run_group','met_cell','ann_value','src_id']
    hierarchy = [['region_cd','scc'],['region_cd',]]
    df = df[cols].copy().sort_values('ann_value').drop_duplicates(['run_group','region_cd','met_cell'], 
      take_last=True)
    temp.profs.fillna(0, inplace=True)
    value_cols = ['month','day','hour','factor']
    run_group = df['run_group'].values[0]
    fips_map = get_fips_map(df[['run_group','region_cd','met_cell','src_id',
      'ann_value']].copy())
    write_county_xwalk(fips_map, run_group)
    df = df[['run_group','region_cd','scc',
      'ann_value']].copy().drop_duplicates(['run_group','region_cd','scc'])
    for state in list(df['region_cd'].str[:2].drop_duplicates()):
        st_df = df.loc[df['region_cd'].str[:2] == state].copy() 
        st_df = match_temporal(st_df, temp.profs, value_cols, hierarchy, temp.use_daily)
        st_df.drop('scc', axis=1, inplace=True)
        st_df['year'] = os.environ['BASE_YEAR'][2:4]
        st_df.drop_duplicates(['region_cd','run_group','year','month','day','hour'], 
          inplace=True)
        keys = ['region_cd','run_group','year']
        mean_factor = st_df[['region_cd','factor']].copy().groupby('region_cd', as_index=False).mean()
        st_df = pd.merge(st_df, mean_factor, on='region_cd', how='left', suffixes=['','_mean'])
        st_df['factor'] = (st_df['factor']/st_df['factor_mean']) 
        st_df[['month','day','hour']] = st_df[['month','day','hour']].astype('i')
        fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_%s_hourly.csv' %(run_group, state))
        st_df.to_csv(fname, columns=['run_group','region_cd','year','month','day','hour','factor'],
          index=False, float_format='%.12g')

# Default temporal profile
def_prof = pd.Series(['HROFDY',0.4751525,0.4463554,0.4367563,0.4487551,0.5039496,0.59994,0.7463254,
  0.9311069,1.120688,1.267073,1.370263,1.449455,1.487851,1.514249,1.523848,1.497450,1.425457,
  1.315068,1.274273,1.221478,1.019898,0.7847215,0.6167383,0.5231477], 
  index=['qflag',]+['Scalar%s' %x for x in range(1,25)])
def get_no_daily_prof(run_df, temp):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    '''
    scalar_cols = [s_col for s_col in temp.profs.columns if s_col.startswith('Scalar')]
    value_cols = ['qflag',] + scalar_cols
    # Only match by fips/scc or fips; SCC only gets default
    hierarchy = [['region_cd','scc'],['region_cd',]]
    run_df = match_temporal(run_df, temp.profs, value_cols, hierarchy)
    run_df.drop(['region_cd','scc'], axis=1, inplace=True)
    run_df.drop_duplicates(inplace=True)
    if not run_df[run_df[def_prof.index[0]].isnull()].empty:
        run_df = fill_default(run_df, def_prof)
    qflag_list = list(run_df['qflag'].drop_duplicates())
    if 'MHRDOW7' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in range(1,2017)]
    elif 'MHRDOW' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in range(1,865)]
    elif 'HROFDY' in qflag_list:
        scalar_cols = ['Scalar%s' %x for x in range(1,25)]
    else:
        scalar_cols = ['Scalar%s' %x for x in range(1,13)]
    return run_df[['run_group','met_cell','src_id','qflag']+scalar_cols].copy()

def write_no_daily_prof(df, temp):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    '''
    cols = ['region_cd','scc','run_group','met_cell','src_id','ann_value']
    df = df[cols].copy().drop_duplicates(['region_cd','run_group','met_cell','src_id'])
    run_group = df['run_group'].values[0]
    df.sort_values('src_id', inplace=True)
    states = list(df['region_cd'].str[:2].drop_duplicates())
    if run_group[:3] in ('LDO','HDO','HOT'):
        fips_map = get_fips_map(df[['run_group','region_cd','met_cell','src_id',
          'ann_value']].copy())
        write_county_xwalk(fips_map, run_group)
        print('WARNING: Temporal profiles for onroad are now run separately using specialized reports.')
        '''
        for state in states:
            st_df = get_no_daily_prof(run_df[run_df['region_cd'].str[1:3] == state], temp)
            write_onroad_hourly(st_df, fips_map[fips_map['region_cd'].str[1:3] == state], 
              run_group, state)
        '''
    else:
        fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_temporal.csv' %run_group)
        with open(fname, 'w') as temp_out:
            for idx, state in enumerate(states):
                st_df = get_no_daily_prof(df[df['region_cd'].str[:2] == state], temp)
                if idx == 0:
                    st_df.to_csv(temp_out, index=False)
                else:
                    st_df.to_csv(temp_out, index=False, header=False)

def get_fips_map(df):
    '''
    Map county and cell by sum of emissions
    '''
    df.sort_values('ann_value', inplace=True)
    df.drop_duplicates(['run_group','met_cell','src_id'], keep='last', inplace=True)
    df.drop('ann_value', axis=1, inplace=True)
    return df

def write_county_xwalk(fips_map, run_group):
    '''
    Write the county to grid cell crosswalk
    '''
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
    cols = ['run_group','region_cd','year','month','dow','hour','scalar']
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_%s_hourly.csv' %(run_group, state))
    df.to_csv(fname, columns=cols, index=False, float_format='%.12g')

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


