from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import csv
import os
import os.path
from math import pi, ceil, log10
import pandas as pd
from smk2ae.temporal import match_temporal 
from smk2ae.utm import UTM
from smk2ae.qa import QA

def proc_points(inv, temp, hrly, work_path, grid_info):
    '''
    Process the non-airports
    '''
    '''
    The parameters are calculated here by pt and fugitive types because the parameters
      are used to define unique source IDs. The pt and fug param definitions are different.
    '''
    qa = QA(work_path)
    qa.sources = inv.stk.shape[0]
    print('Total Input Sources: %s' %qa.sources)
    pt_df = get_point_params(inv.stk[inv.stk['erptype'] != 1].copy())
    fug_df = get_fug_params(inv.stk[inv.stk['erptype'] == 1].copy())
    src_map = get_src_map(pt_df, fug_df)
    inv.stk = pd.merge(inv.stk, src_map, on=['facility_id','unit_id','rel_point_id','scc'], 
      how='left')
    write_relpt_xwalk(inv.stk, work_path)
    inv.stk.drop_duplicates(['facility_id','src_id'], inplace=True)
    qa.uniq_srcs = inv.stk.shape[0]
    print('Total Unique Source IDs Across all facilities: %s' %qa.uniq_srcs)
    inv.stk = calc_predom_cell(inv, qa)
    inv.stk = get_utm_location(inv.stk, grid_info)
    write_src_xwalk(inv, src_map, work_path)
    write_locations(inv.stk, work_path)
    pt_df = pd.merge(pt_df, src_map, on=['facility_id','unit_id','rel_point_id','scc'], how='left')
    write_point_params(pt_df, work_path)
    fug_df = pd.merge(fug_df, src_map, on=['facility_id','unit_id','rel_point_id','scc'], how='left')
    qa.points = pt_df.shape[0]
    qa.fugs = fug_df.shape[0]
    write_fug_params(fug_df, work_path)
    if not hrly.egu_units.empty:
        hrly_stks = pd.merge(inv.stk[inv.stk['ALLDAY'] == -1], 
          pt_df[['facility_id','src_id','velocity','temp']].drop_duplicates(), 
          on=['facility_id','src_id'], how='left')
        write_hourly_factors(hrly_stks, hrly.aermod, os.environ['WORK_PATH'])
        qa.hourly = hrly_stks.shape[0]
    inv.stk = inv.stk[inv.stk['ALLDAY'] != -1]
    qa.n_temp = inv.stk.shape[0]
    if not inv.stk.empty:
        write_temp_factors(inv.stk, temp.aermod, os.environ['WORK_PATH'])
    qa.total_benzene = sum(inv.emis.loc[inv.emis['smoke_name'] == 'BENZENE','ann_value'])
    qa.write_totals(inv.state)

def get_grid_location(df, grid_info):
    '''
    Generate the point locations file
    Fill with grid projection information and UTM information
    '''
    df[['grid_x','grid_y']] = df[['longitude','latitude']].apply(lambda x: \
      pd.Series(grid_info.get_coords(x['longitude'],x['latitude'])), axis=1)
    df[['col','row']] = df[['grid_x','grid_y']].apply(lambda x: \
      pd.Series(grid_info.get_cell(x['grid_x'],x['grid_y'])), axis=1)
    return df

def get_utm_location(df, grid_info):
    '''
    Select UTM zone based on centroid of met grid cell. Calc UTM X and Y.
    '''
    utm = UTM()
    df['centroid_lon'] = grid_info.colrow_to_ll(df['col'].astype('f') + 0.5, 
      df['row'].astype('f') + 0.5)['lon'].astype('f')
    utm = UTM()
    df['utm_zone'] = df['centroid_lon'].apply(utm.get_zone)
    df[['utm_x','utm_y']] = df[['longitude','latitude','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['longitude'],x['latitude'],x['utm_zone'])), axis=1)
    return df

def write_locations(df, work_path):
    '''
    Write the unique locations of each source ID by facility
    '''
    cols = ['state','facility_id','facility_name','src_id','grid_x','grid_y','longitude',
        'latitude','utm_x','utm_y','utm_zone','col','row']
    df = df[cols].copy()
    fname = os.path.join(work_path, 'locations', 'point_locations.csv')
    df.to_csv(fname, columns=cols, index=False, quotechar=' ')

def calc_predom_cell(inv, qa):
    '''
    Calculate the predominate grid cell for a facility based on HAPs emissions
    '''
    for fac in list(inv.stk['facility_id'].drop_duplicates()):
        fac_df = inv.stk.loc[inv.stk['facility_id'] == fac, ['col','row','ann_value']].copy()
        fac_df = fac_df.groupby(['col','row'], as_index=False).sum()
        if len(fac_df[['col','row']]) > 1:
            fac_df.sort_values('ann_value', ascending=False, inplace=True)
            fac_df.drop_duplicates(['col','row'],inplace=True)
            qa.moved_facs.append(str(fac))
            inv.stk.loc[inv.stk['facility_id'] == fac, 'col'] = fac_df['col'].values[0]
            inv.stk.loc[inv.stk['facility_id'] == fac, 'row'] = fac_df['row'].values[0]
    qa.write_moved_facs(inv.state)
    return inv.stk

def get_point_params(df):
    '''
    Calculate into metric units
    '''
    df.drop(['fug_height','fug_width_ydim','fug_length_xdim'],axis=1,inplace=True)
    df['height'] = 0.3048 * df['stkhgt']
    df['diameter'] = 0.3048 * df['stkdiam']
    df['velocity'] = 0.3048 * df['stkvel']
    df['temp'] = (df['stktemp'] + 459.67) * (old_div(5.,9.))
    src_types = {'2': 'POINT', '3': 'POINTHOR', '4': 'POINT', '5': 'POINTCAP', '6': 'POINTHOR', 
      '99': 'POINT'}
    df['aermod_src_type'] = df['erptype'].apply(lambda x: src_types[str(int(x))])
    return df

def write_point_params(df, work_path):
    '''
    Write the unique point parameters of each source ID by facility
    '''
    cols = ['facility_id','facility_name','src_id','aermod_src_type','height',
        'temp','velocity','diameter']
    df = df[cols].copy().drop_duplicates(['facility_id','src_id'])
    fname = os.path.join(work_path, 'parameters', 'point_srcparam.csv')
    df.to_csv(fname, columns=cols, index=False, quotechar=' ')

def get_fug_params(df):
    '''
    Fill the fugitive source parameters; with defaults as needed 
    '''
    df.drop(['stkhgt','stktemp'],axis=1,inplace=True)
    df['aermod_src_type'] = 'AREA'
    df['rel_ht'] = df['fug_height'] * 0.3048
    df['x_length'] = df['fug_width_ydim'] * 0.3048
    df['y_length'] = df['fug_length_xdim'] * 0.3048
    df.loc[df['rel_ht'] > 10, 'szinit'] = old_div(df.loc[df['rel_ht'] > 10, 'rel_ht'], 4.3)
    df.loc[df['rel_ht'] <= 10, 'szinit'] = 0
    df.rename(columns={'fug_angle': 'angle'}, inplace=True)
    df.loc[df['angle'] == -99., 'angle'] = 0
    return df

def write_fug_params(df, work_path):
    '''
    Write the unique fugitive parameters of each source ID by facility
    '''
    cols = ['facility_id','facility_name','src_id','aermod_src_type','rel_ht','x_length',
        'y_length','angle','szinit']
    df = df[cols].copy().drop_duplicates(['facility_id','src_id'])
    fname = os.path.join(work_path, 'parameters', 'point_fug_srcparam.csv')
    df.to_csv(fname, columns=cols, index=False, quotechar=' ')

def get_temp_codes(df, xref):
    '''
    Get the temporal profile codes for each row. This is used for checking unique temporal profile.

    *  Assume that only diurnal profiles only vary by weekday/weekend or
        are constant across all days of the week. 
       This assumption is good for all US point sources in 2014v2.
    '''
    if len(xref[xref['ALLDAY'].isnull()]) > 0:
        if 'TUESDAY' in list(xref.columns):
            try:
                xref.loc[xref['ALLDAY'].isnull(), 'ALLDAY'] = xref.loc[xref['ALLDAY'].isnull(), 
                  'TUESDAY']
            except KeyError:
                xref.loc[xref['ALLDAY'].isnull(), 'ALLDAY'] = xref.loc[xref['ALLDAY'].isnull(), 
                  'WEEKDAY']
    type_cols = ['ALLDAY','MONTHLY','WEEKLY']
    xref.drop_duplicates(['facility_id','region_cd','scc'], inplace=True)
    xref['region_cd'] = xref['region_cd'].apply(fix_fips)
    xref['scc'] = xref['scc'].apply(lambda x: str(x).strip())
    xref['facility_id'] = xref['facility_id'].apply(lambda x: str(x).strip())
    hierarchy = [['region_cd','scc','facility_id'],['scc','facility_id'],['region_cd','facility_id'],
        ['facility_id',], ['region_cd','scc'], ['scc',], ['region_cd',]]
    return match_temporal(df, xref, type_cols, hierarchy)

def write_temp_factors(df, temp, work_path):
    '''
    Write the temporal factors by source ID for each facility
    Adjust the output columns based on the max number of scalars used within the facility
    '''
    df = df[['region_cd','scc','facility_id','facility_name','src_id','state']].copy()
    scalar_cols = [col for col in temp.columns if col.startswith('Scalar')]
    hierarchy = [['region_cd','scc','facility_id'],['scc','facility_id'],['region_cd','facility_id'],
        ['facility_id',], ['region_cd','scc'], ['scc',], ['region_cd',]]
    df = match_temporal(df, temp, ['qflag',]+scalar_cols, hierarchy)
    cols = ['facility_id','facility_name','src_id','qflag']
    df = df[cols+scalar_cols].copy().drop_duplicates(['facility_id','src_id'])
    qflag_list = list(df['qflag'].drop_duplicates())
    if 'MHRDOW7' in qflag_list:
        col_names = cols + ['Scalar%s' %x for x in range(1,2017)]
    elif 'MHRDOW' in qflag_list:
        col_names = cols + ['Scalar%s' %x for x in range(1,865)]
    elif 'HROFDAY' in qflag_list:
        col_names = cols + ['Scalar%s' %x for x in range(1,25)]
    else:
        col_names = cols + ['Scalar%s' %x for x in range(1,13)]
    fname = os.path.join(work_path, 'temporal', 'point_temporal.csv')
    df.to_csv(fname, columns=col_names, index=False, quotechar=' ')

def get_src_map(pt_df, fug_df):
    '''
    Get a map from the original release point and ID to the unique source ID
    '''
    uniq_cols = ['facility_id','latitude','longitude','erptype','diameter','velocity','temp',
        'height','angle','x_length','y_length']
    if 'MONTHLY' in pt_df.columns:
        uniq_cols.extend(['MONTHLY','WEEKLY','ALLDAY'])
    df = pd.concat((pt_df, fug_df))
    df = df[['rel_point_id','scc','unit_id','process_id']+uniq_cols].copy()
    uniq_stks = df[uniq_cols].copy().drop_duplicates()
    for fac in list(uniq_stks['facility_id'].drop_duplicates()):
        n_srcs = len(uniq_stks[uniq_stks['facility_id'] == fac])
        srcs = ['S%s' %str(int(x)).zfill(int(ceil(log10(n_srcs)))) for x in range(1,n_srcs+1)]
        uniq_stks.loc[uniq_stks['facility_id'] == fac, 'src_id'] = srcs
    df = pd.merge(df, uniq_stks, on=uniq_cols, how='left')
    return df[['facility_id','unit_id','rel_point_id','scc','src_id']].copy().drop_duplicates()

def write_relpt_xwalk(df, work_path):
    '''
    Write release point to src_id x-walk
    '''
    fname = os.path.join(work_path,'xwalk','point_srcid_xwalk.csv')
    cols = ['state','facility_id','facility_name','fac_source_type','unit_id','process_id',
      'rel_point_id','src_id']
    df = df[cols].copy().drop_duplicates()
    df.to_csv(fname, columns=cols, index=False, quotechar=' ')

def write_src_xwalk(inv, src_map, work_path):
    '''
    Write the SRC_ID x-walk
    '''
    df = pd.merge(inv.emis, src_map, on=['facility_id','unit_id','rel_point_id','scc'], how='left')
    df = pd.merge(df, inv.stk[['facility_id','facility_name','fac_source_type','state']].drop_duplicates(), 
      on='facility_id', how='left')
    df = df.groupby(['state','facility_id','facility_name','fac_source_type','src_id','smoke_name'], 
      as_index=False).sum()
    fname = os.path.join(work_path,'xwalk','point_srcid_emis.csv')
    cols = ['state','facility_id','facility_name','fac_source_type','src_id','smoke_name','ann_value']
    df.to_csv(fname, columns=cols, index=False, quotechar=' ')

def write_hourly_factors(df, hr_facs, work_path):
    '''
    Write the hourly factors, typically for EGU
    '''
    cols = ['facility_id','src_id','year','month','day','hour','hour_factor','temp','velocity']
    df = df[['facility_id','unit_id','src_id','state','temp','velocity']].copy().drop_duplicates()
    df = pd.merge(df, hr_facs, on=['facility_id','unit_id'], how='left')
    df['year'] = df['year'].apply(lambda x: str(x)[2:])
    fac_list = sorted(list(df['facility_id'].drop_duplicates()))
    for fac in fac_list:
        fac_df = df[df['facility_id'] == fac].copy()
        state = fac_df['state'].values[0]
        fname = os.path.join(work_path,'temporal','%s_%s_hourly.csv' %(fac,state))
        fac_df.sort_values(['year','month','day','hour'], inplace=True)
        fac_df.to_csv(fname, index=False, columns=cols)

def fix_fips(fips):
    '''
    Force FIPS to a 5 character zero-padded string
    '''
    if fips:
        try:
            fips = '%0.5d' %int(fips)
        except:
            fips = str(fips)
        return fips
    else:
        return str(fips)

