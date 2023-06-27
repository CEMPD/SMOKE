#!/usr/bin/env python
'''
QA for AERMOD point helper files
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

def check_fug_params(aer, inv):
    '''
    Check the fugitive parameters against the inventory
    '''
    inv = inv[inv['erptype'] == 1].copy()
    cols = ['facility_id','src_id','fug_height','fug_length_xdim','fug_width_ydim','fug_angle']
    # Convert meters to feet
    aer['fug_height'] = aer['rel_ht'] * 3.2808
    aer['fug_length_xdim'] = aer['x_length'] * 3.2808
    aer['fug_width_ydim'] = aer['y_length'] * 3.2808
    aer['fug_angle'] = aer['angle']
    keys = ['facility_id','src_id']
    df = pd.merge(aer[cols].drop_duplicates(keys), inv[cols].drop_duplicates(keys), on=keys, 
      how='outer', suffixes=['_aer','_inv'])
    print('\nFug len/wid <= 0:')
    print(aer[(aer['fug_length_xdim'] <= 0) | (aer['fug_width_ydim'] <= 0)])
    return df

def check_point_params(aer, inv):
    '''
    Check the stack parameters against the inventory
    '''
    inv = inv[inv['erptype'] != 1].copy()
    cols = ['facility_id','src_id','stkhgt','stktemp','stkvel','stkdiam']
    # Convert meters (used in AERMOD helpers) to feet (using in inventory)
    aer['stkhgt'] = aer['height'] * 3.2808
    aer['stkvel'] = aer['velocity'] * 3.2808
    aer['stkdiam'] = aer['diameter'] * 3.2808
    aer['stktemp'] = aer['temp'] * (9/5.) - 459.67  # K to degF
    keys = ['facility_id','src_id']
    df = pd.merge(aer[cols].drop_duplicates(keys), inv[cols].drop_duplicates(keys), on=keys, 
      how='outer', suffixes=['_aer','_inv'])
    print('\nStack parameters <= 0:')
    print(aer[(aer['stkhgt'] <= 0) | (aer['stkvel'] <= 0) | (aer['stkdiam'] <= 0)])
    return df

def loc_check(sector, work_path):
    # Location file QA
    locs = os.path.join(work_path, 'locations', '%s_location.csv' %sector)
    usecols=['facility_id','src_id','grid_x','grid_y','longitude','latitude','utm_x','utm_y','utm_zone','col','row']
    locs = pd.read_csv(locs, dtype={'facility_id': str}, usecols=usecols, skipinitialspace=True)
    # Check for duplicate sources
    print('\nDuplicate sources in locations:')
    print(locs[locs[['facility_id','src_id']].duplicated()])
    # Check to see if a facility is in more than one column and row
    fcr = locs[['facility_id','col','row']].drop_duplicates()
    print('\nFacilities in multiple col/row:')
    print(fcr[fcr.duplicated('facility_id')])
    # Check to see if a facility is in more than UTM zone
    fcr = locs[['facility_id','utm_zone']].drop_duplicates()
    print('\nFacilities in multiple UTM zones:')
    print(fcr[fcr.duplicated('facility_id')])
    # Check UTM locations
    utmd = locs[['facility_id','utm_x','utm_y','utm_zone','latitude','longitude']].drop_duplicates('facility_id')
    utmp = UTM()
    utmd['zone_check'] = utmd['longitude'].apply(lambda lon: utmp.get_zone(lon))
    print('\nNon matching UTM zones:')
    print(utmd[utmd['utm_zone'].astype(int) != utmd['zone_check'].astype(int)])
    utmd[['check_x','check_y']] = utmd[['longitude','latitude','utm_zone']].apply(lambda x: \
      pd.Series(utmp.get_coords(x['longitude'],x['latitude'],x['utm_zone'])), axis=1)
    print('\nNon matching UTM coords:')
    print(utmd[(utmd['utm_x'].astype('f').round(0) != utmd['check_x'].astype('f').round(0)) \
      | (utmd['utm_y'].astype('f').round(0) != utmd['check_y'].astype('f').round(0))])
    return locs 

def fug_params_check(sector, work_path):
    fug_params = pd.DataFrame()
    # Aggregate and count parameters
    fn = os.path.join(work_path, 'parameters', '%s_fug_srcparam.csv' %sector)
    if os.path.exists(fn):
        fug_params = pd.read_csv(fn, dtype={'facility_id': str}, 
          usecols=['facility_id','src_id','aermod_src_type','rel_ht','x_length','y_length','angle'], 
          skipinitialspace=True)
        if len(fug_params[fug_params['aermod_src_type'] != 'AREA']) > 0:
            print('WARNING: Source type other than AREA found in fugitive parameter file.')
    return fug_params

def point_params_check(sector, work_path):
    point_params = pd.DataFrame()
    fn = os.path.join(work_path, 'parameters', '%s_point_srcparam.csv' %sector)
    if os.path.exists(fn):
        point_params = pd.read_csv(fn, dtype={'facility_id': str}, 
          usecols=['facility_id','src_id','aermod_src_type','height','temp','velocity','diameter'], 
          skipinitialspace=True)
        if len(point_params[point_params['aermod_src_type'] == 'AREA']) > 0:
            print('WARNING: Source type AREA found in point parameter file.')
    return point_params

def check_temp_hr(sector, work_path):
    # Read in and do a count of temporal files
    hr = pd.DataFrame()
    if sector.startswith('ptegu') or 'combined' in sector:
        from glob import glob
        temp_glob = '%s/temporal/[0-9]*_hourly.csv' %work_path
        temp = pd.DataFrame()
        for fn in glob(temp_glob):
            df = pd.read_csv(fn, dtype={'facility_id': str}, usecols=['facility_id','src_id','hour_factor'], 
              skipinitialspace=True)
            df = df.groupby(['facility_id','src_id'], as_index=False).agg(['sum','count'])
            df.reset_index(inplace=True)
            df.columns = ['facility_id','src_id','sum','count']
            df['sum'] = df['sum'].round(2)
            temp = pd.concat((temp, df))
            temp.to_csv(os.path.join(work_path, 'qa', 'temporal_hour_%s_check.csv' %sector), index=False)
        hr = temp[['facility_id','src_id']].copy()
    return hr

def check_temporal(sector, work_path):
    temp = pd.DataFrame()
    if 'ptegu' not in sector:
        # QA the temporal file
        temp = os.path.join(work_path, 'temporal', '%s_temporal.csv' %sector)
        temp = pd.read_csv(temp, dtype={'facility_id': str}, skipinitialspace=True, index_col=False)
        scalars = [x for x in temp.columns if x.startswith('Scalar')]
        temp.loc[temp['qflag'] == 'MONTH', 'sum'] = temp.loc[temp['qflag'] == 'MONTH', scalars[:12]].sum(axis=1)
        temp.loc[temp['qflag'] == 'HROFDY', 'sum'] = temp.loc[temp['qflag'] == 'HROFDY', scalars[:24]].sum(axis=1)
        temp.loc[temp['qflag'] == 'MHRDOW7', 'sum'] = temp.loc[temp['qflag'] == 'MHRDOW7', scalars].sum(axis=1) * \
          (8760./2016.)
        if 'MHRDOW' in list(temp['qflag']):
            wkday = [scalars[x] for x in range(288)]
            sat = [scalars[x] for x in range(288,576)]
            sun = [scalars[x] for x in range(576,864)]
            temp.loc[temp['qflag'] == 'MHRDOW', 'sum'] = (5 * temp.loc[temp['qflag'] == 'MHRDOW', wkday].sum(axis=1) + \
              temp.loc[temp['qflag'] == 'MHRDOW', sat].sum(axis=1) + temp.loc[temp['qflag'] == 'MHRDOW', sun].sum(axis=1)) \
              * (8760./2016.)
        temp['sum'] = temp['sum'].round(4)
        temp.to_csv(os.path.join(work_path, 'qa', 'temporal_%s_check.csv' %sector), index=False, 
          columns=['facility_id','src_id','qflag','sum'])
        temp = temp[['facility_id','src_id']].copy()
    return temp

def check_emis(sector, work_path):
    # Read in the point emissions helper file
    emis = os.path.join(work_path, 'xwalk', '%s_srcid_emis.csv' %sector)
    emis = pd.read_csv(emis, dtype={'facility_id': str}, usecols=['facility_id','src_id','smoke_name','ann_value'], 
      skipinitialspace=True) 
    emis = emis.groupby(['facility_id','src_id','smoke_name'], as_index=False).sum()
    return emis

def get_invtable(fn):
    # Read in the inventory table to get the pollutants used in the AERMOD run
    invtable = pd.read_fwf(fn, comment='#', colspecs=[(0,11), (16,32), (41,42), (43,49)], 
      names=['smoke_name','poll','keep','spec_factor'], 
      converters={'smoke_name': str.strip, 'poll': str.strip, 'keep': str.strip})
    invtable = invtable[invtable['keep'].str.upper() == 'Y'].copy()
    invtable.drop_duplicates(inplace=True)
    return invtable

def get_inv(sector, inv_list, invtable):
    # Read in all of the inventories from the inventory list
    # Select only the sources that have the selected HAPs and are not airports
    polls = list(invtable['poll'].drop_duplicates())
    inv = pd.DataFrame()
    usecols = ['facility_id','poll','ann_value','fac_source_type','rel_point_id','unit_id',
      'process_id','region_cd','erptype','stkhgt','stkvel','stkdiam','stktemp','fug_height',
      'fug_angle','fug_width_ydim','fug_length_xdim']
    dtype = {'poll': str, 'facility_id': str, 'fac_source_type': str, 'rel_point_id': str, 
      'unit_id': str, 'process_id': str, 'region_cd': str}
    for fn in inv_list:
        with open(fn) as f:
            cmt_cnt = 0
            for l in f:
                if l.startswith('#'):
                    cmt_cnt = cmt_cnt + 1
                else:
                    break
        df = pd.read_csv(fn, skiprows=cmt_cnt, usecols=usecols, dtype=dtype) 
        df = df[(df['poll'].isin(polls)) & (df['fac_source_type'] != '100') & (df['ann_value'] > 0)].copy()
        if sector in ('point','ptegu','point_conus'):
            # Remove inv values outside of the CONUS domain
            df = df[~ df['region_cd'].astype(int).astype(str).str.zfill(5).str[:2].isin(('02','15','72','78'))].copy()
        elif 'nonconus' in sector:
            df = df[df['region_cd'].astype(int).astype(str).str.zfill(5).str[:2].isin(('02','15','72','78'))].copy()
        df = pd.merge(df, invtable, on='poll', how='left')
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        inv = pd.concat((inv, df))
    return inv


def main():
    sector = os.environ['SECTOR']
    work_path = os.environ['WORK_PATH'] 
    invtable = os.environ['INVTABLE']
    # Determine which sector type to use
    if sector.startswith('ptegu'):
        sector = 'ptegu'
    elif 'combined' in sector:
        sector = 'point_combined'
    elif '_nonconus' in sector:
        sector = 'point_nonconus'
    elif '_conus' in sector:
        sector = 'point_conus'
    else:
        sector = 'point'
    if not os.path.exists(os.path.join(work_path, 'qa')):
        os.makedirs(os.path.join(work_path, 'qa'))
    qa_counts = []
    # Check the locations file
    locs = loc_check(sector, work_path)
    # Sum up locations by facility and source
    fac_locs, src_locs = (len(locs[['facility_id','src_id']].drop_duplicates('facility_id')), len(locs))
    print('\nFacility locations: %s    Source locations: %s' %(fac_locs, src_locs))
    qa_counts += [str(fac_locs), str(src_locs)]
    # Check the fugitive parameters file
    fug_params = fug_params_check(sector, work_path)
    if len(fug_params) > 0:
        fac_params, src_params = (len(fug_params.drop_duplicates('facility_id')), len(fug_params)) 
        print('Fugitive facilities: %s   Fugitive sources: %s' %(fac_params, src_params))
    else:
        fac_params, src_params = [0,0]
    qa_counts += [str(fac_params), str(src_params)]
    #Check the point parameters file
    point_params = point_params_check(sector, work_path)
    if len(point_params) > 0:
        fac_params, src_params = (len(point_params.drop_duplicates('facility_id')), len(point_params))
        print('Point facilities: %s   Point sources: %s' %(fac_params, src_params))
    else:
        fac_params, src_params = [0,0]
    qa_counts += [str(fac_params), str(src_params)]
    params = pd.concat((fug_params[['facility_id','src_id']], point_params[['facility_id','src_id']]))
    fac_params, src_params = (len(params.drop_duplicates('facility_id')), len(params))
    print('Total unique param facilities: %s   Total unique param sources: %s' %(fac_params, src_params))
    qa_counts += [str(fac_params), str(src_params)]
    # Check the temporal files
    hr = check_temp_hr(sector, work_path)
    temp = check_temporal(sector, work_path)
    temp = pd.concat((temp, hr))
    fac_temp, src_temp = (len(temp[['facility_id','src_id']].drop_duplicates('facility_id')), len(temp))
    print('Temporal facilities: %s   Temporal sources: %s' %(fac_temp, src_temp))
    qa_counts += [str(fac_temp), str(src_temp)]
    # Check if the sources in the locations file are also in the temporal file and vice-versa
    locs['z'] = 'Y'
    temp['s'] = 'Y'
    temp = pd.merge(temp, locs, on=['facility_id','src_id'], how='outer')
    temp[(temp['s'].isnull()) | (temp['z'].isnull())].to_csv(os.path.join(work_path, 'qa', 
      'missing_%s_sources.csv' %sector), index=False)
    # Get the emis information
    emis = check_emis(sector, work_path)
    fac_emis, src_emis = (len(emis.drop_duplicates('facility_id')), 
        len(emis.drop_duplicates(['facility_id','src_id'])))
    print('Xwalk facilities: %s   Xwalk sources: %s' %(fac_emis, src_emis))
    qa_counts += [str(fac_emis), str(src_emis)]
    # Get the inventory emissions
    invtable = get_invtable(invtable)
    inv_list = get_inv_list()
    inv = get_inv(sector, inv_list, invtable)
    # Do a pollutant emissions comparison
    inv_polls = inv[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
    emis_polls = emis[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
    polls = pd.merge(emis_polls, inv_polls, on='smoke_name', how='outer', suffixes=['_aer','_inv'])
    qa_out = os.path.join(work_path, 'qa', '%s_poll_qa.csv' %sector)
    polls.to_csv(qa_out, index=False)
    # Count the number of unique facilities in the inventories
    inv_facs = len(inv.drop_duplicates('facility_id'))
    print('Inv facilities: %s ' %inv_facs)
    qa_counts.append(str(inv_facs))
    # Write the facility and source counts to the "counts"  QA file
    count_out = os.path.join(work_path, 'qa', '%s_counts.csv' %sector)
    with open(count_out, 'w') as f:
        f.write('fac_locs,src_locs,fug_fac_parms,src_fac_parms,pt_fac_parms,pt_src_parms,tot_fac_parms,tot_src_parms,fac_temp,src_temp,fac_emis,src_emis,fac_inv\n')
        f.write('%s\n' %','.join(qa_counts))
    # Read in the point source xwalk helper file
    xwalk = os.path.join(work_path, 'xwalk', '%s_srcid_xwalk.csv' %sector)
    xwalk = pd.read_csv(xwalk, usecols=['facility_id','unit_id','process_id','rel_point_id','src_id'],
      dtype={'facility_id': str, 'rel_point_id': str, 'unit_id': str, 'process_id': str})
    inv = pd.merge(inv, xwalk, on=['facility_id','unit_id','process_id','rel_point_id'], how='left')
    inv['src_id'] = inv['src_id'].fillna('')
    # Use the xwalk file to go from inventory sources to AERMOD sources
    inv_emis = inv[['facility_id','src_id','ann_value']].groupby(['facility_id','src_id'], as_index=False).sum()
    emis = emis[['facility_id','src_id','ann_value']].groupby(['facility_id','src_id'], as_index=False).sum()
    # Write an emissions comparison for the sources between the inventory and helper files
    df = pd.merge(inv_emis, emis, on=['facility_id','src_id'], how='outer', suffixes=['_inv','_aer'])
    print('\nMissing facility:')
    print(df[df['ann_value_aer'].isnull()])
    qa_out = os.path.join(work_path, 'qa', '%s_emis_qa.csv' %sector)
    df.to_csv(qa_out, index=False, columns=['facility_id','src_id','ann_value_inv','ann_value_aer'])
    # Stack/fugitive parameter validation and QA
    inv.drop_duplicates(['facility_id','src_id'], inplace=True)
    if len(fug_params) > 0:
        df = check_fug_params(fug_params, inv)
        qa_out = os.path.join(work_path, 'qa', '%s_fug_param_qa.csv' %sector)
        df.to_csv(qa_out, index=False)
    if len(point_params) > 0:
        df = check_point_params(point_params, inv)
        qa_out = os.path.join(work_path, 'qa', '%s_point_param_qa.csv' %sector)
        df.to_csv(qa_out, index=False)
    # AERMOD Source crosscheck - verify that all AERMOD sources are in all helper files
    xcheck = pd.DataFrame()
    for key,df in {'locs': locs, 'params': params, 'temp': temp, 'emis': emis, 'xwalk': xwalk}.items():
        df = df[['facility_id','src_id']].drop_duplicates()
        df[key] = 'Y'
        if xcheck.empty:
            xcheck = df.copy()
        else:
            xcheck = pd.merge(xcheck, df, on=['facility_id','src_id'], how='outer')
    xcheck = xcheck.fillna('N')
    qa_out = os.path.join(work_path, 'qa', '%s_srcid_qa.csv' %sector)
    xcheck.to_csv(qa_out, index=False)

if __name__ == '__main__':
    main()
