#!/usr/bin/env python
'''
QA for AERMOD airport files
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

work_path = os.environ['WORK_PATH'] 
invtable = os.environ['INVTABLE']

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
    usecols = ['facility_id','poll','ann_value','fac_source_type','rel_point_id','unit_id',
      'process_id','region_cd','scc']
    dtype = {'poll': str, 'facility_id': str, 'fac_source_type': str, 'rel_point_id': str, 
      'unit_id': str, 'process_id': str, 'region_cd': str, 'scc': str}
    for fn in inv_list:
        with open(fn) as f:
            cmt_cnt = 0
            for l in f:
                if l.startswith('#'):
                    cmt_cnt = cmt_cnt + 1
                else:
                    break
        df = pd.read_csv(fn, skiprows=cmt_cnt, usecols=usecols, dtype=dtype) 
        df = df[(df['poll'].isin(polls)) & (df['fac_source_type'] == '100') & (df['ann_value'] > 0)].copy()
        if len(states[states['state'].isin(('02','15','72','78'))]) == 0:
            # Remove inv values outside of the CONUS domain
            df = df[~ df['region_cd'].astype(int).astype(str).str.zfill(5).str[:2].isin(('02','15','72','78'))].copy()
        # Cut lead in half for these two SCCs
        df.ix[(df['scc'].isin(('2275050011','2275060011'))) & (df['poll'] == '7439921'), 'ann_value'] = \
          df.ix[df['scc'].isin(('2275050011','2275060011')) & (df['poll'] == '7439921'), 'ann_value'] * 0.5
        df = pd.merge(df, invtable, on='poll', how='left')
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        inv = pd.concat((inv, df))
    return inv

if not os.path.exists(os.path.join(work_path, 'qa')):
    os.makedirs(os.path.join(work_path, 'qa'))

runways = pd.read_csv(os.environ['RUNWAYS'], usecols=['facility_id','siteno'],
  dtype={'facility_id': str})

# Location file QA
qa_counts = []
line_locs = os.path.join(work_path, 'locations', 'airport_combined_line_locations.csv')
usecols = ['facility_id','src_id','col','row','utm_zone']
line_locs = pd.read_csv(line_locs, dtype={'facility_id': str}, usecols=usecols, skipinitialspace=True)
# Check runways
rw_check = pd.merge(runways, line_locs[['facility_id','src_id']], on='facility_id', how='outer')
print('Runway x-check:')
print(rw_check[(rw_check['src_id'].isnull()) | (rw_check['siteno'].isnull())])
# Sum up line locations by facility and source
fac_locs, src_locs = (len(line_locs[['facility_id','src_id']].drop_duplicates('facility_id')), len(line_locs))
print('Line locations: %s    Line source locations: %s' %(fac_locs, src_locs))
qa_counts += [str(fac_locs), str(src_locs)]
locs = os.path.join(work_path, 'locations', 'airport_combined_nonrunway_locations.csv')
usecols = ['facility_id','src_id','utm_zone','col','row','utm_x','utm_y','latitude','longitude']
locs = pd.read_csv(locs, dtype={'facility_id': str}, usecols=usecols, skipinitialspace=True)
# Check nonrunways
rw_check = pd.merge(runways, locs[['facility_id','src_id']], on='facility_id', how='inner')
print('Nonrunway x-check:')
print(rw_check)
# Sum up nonrunway locations by facility and source
fac_locs, src_locs = (len(locs[['facility_id','src_id']].drop_duplicates('facility_id')), len(locs))
print('Nonrunway locations: %s    Nonrunway source locations: %s' %(fac_locs, src_locs))
qa_counts += [str(fac_locs), str(src_locs)]
# Check UTM locations for nonrunways
utmd = locs[['facility_id','utm_x','utm_y','utm_zone','latitude','longitude']].drop_duplicates('facility_id')
utmp = UTM()
utmd['zone_check'] = utmd['longitude'].apply(lambda lon: utmp.get_zone(lon))
print('Non matching UTM zones:')
print(utmd[utmd['utm_zone'].astype(int) != utmd['zone_check'].astype(int)])
utmd[['check_x','check_y']] = utmd[['longitude','latitude','utm_zone']].apply(lambda x: \
  pd.Series(utmp.get_coords(x['longitude'],x['latitude'],x['utm_zone'])), axis=1)
print('Non matching UTM coords:')
print(utmd[(utmd['utm_x'].astype('f').round(0) != utmd['check_x'].astype('f').round(0)) \
  | (utmd['utm_y'].astype('f').round(0) != utmd['check_y'].astype('f').round(0))])
# Check to see if a facility is in more than column and row
locs = pd.concat((line_locs, locs[['facility_id','src_id','col','row','utm_zone']]))
fcr = locs[['facility_id','col','row']].drop_duplicates()
print('Facilities in multiple col/row:')
print(fcr[fcr.duplicated('facility_id')])
# Check to see if a facility is in more than UTM zone
fcr = locs[['facility_id','utm_zone']].drop_duplicates()
print('Facilities in multiple UTM zones:')
print(fcr[fcr.duplicated('facility_id')])
# Sum up locations by facility and source
fac_locs, src_locs = (len(locs[['facility_id','src_id']].drop_duplicates('facility_id')), len(locs))
print('Facility locations: %s    Source locations: %s' %(fac_locs, src_locs))
qa_counts += [str(fac_locs), str(src_locs)]

# Aggregate and count parameters
fug_params = os.path.join(work_path, 'parameters', 'airport_combined_nonrunway_params.csv')
if os.path.exists(fug_params):
    fug_params = pd.read_csv(fug_params, dtype={'facility_id': str}, 
      usecols=['facility_id','src_id','relhgt','lengthx','lengthy'], 
      skipinitialspace=True)
    if len(fug_params[(fug_params['lengthx'] != 10) | (fug_params['lengthy'] != 10)]) > 0:
        print('WARNING: Nonrunway x and y != 10')
    if len(fug_params[(fug_params['relhgt'] != 3)]) > 0:
        print('WARNING: Nonrunway relhgt != 3')
    fug_count = fug_params.drop_duplicates(['facility_id','src_id'])
    fac_params, src_params = (len(fug_count.drop_duplicates('facility_id')), len(fug_count)) 
    print('Nonrunway facilities: %s   Nonrunway sources: %s' %(fac_params, src_params))
else:
    fac_params, src_params = [0,0]
    fug_params = pd.DataFrame()
qa_counts += [str(fac_params), str(src_params)]
line_params = os.path.join(work_path, 'parameters', 'airport_combined_line_params.csv')
if os.path.exists(line_params):
    line_params = pd.read_csv(line_params, dtype={'facility_id': str}, 
      usecols=['facility_id','src_id','src_type','relhgt','width'], 
      skipinitialspace=True)
    if len(line_params[line_params['src_type'] != 'LINE']) > 0:
        print('WARNING: Source type other than LINE found in line parameter file.')
    if len(line_params[(line_params['relhgt'] != 3)]) > 0:
        print('WARNING: Line relhgt != 3')
    if len(line_params[(line_params['width'] != 25) & (line_params['width'] != 50)]) > 0:
        print('WARNING: Invalid line width')
    line_count = line_params.drop_duplicates(['facility_id','src_id'])
    fac_params, src_params = (len(line_count.drop_duplicates('facility_id')), len(line_count))
    print('Point facilities: %s   Point sources: %s' %(fac_params, src_params))
else:
    fac_params, src_params = [0,0]
    line_params = pd.DataFrame()
qa_counts += [str(fac_params), str(src_params)]
params = pd.concat((fug_count, line_count))
fac_params, src_params = (len(params.drop_duplicates('facility_id')), len(params))
print('Total unique param facilities: %s   Total unique param sources: %s' %(fac_params, src_params))
qa_counts += [str(fac_params), str(src_params)]

# Read in and do a count of temporal files
# QA the temporal file
temp = os.path.join(work_path, 'temporal', 'airport_combined_line_temporal.csv')
temp = pd.read_csv(temp, dtype={'facility_id': str}, skipinitialspace=True, index_col=False)
nr = os.path.join(work_path, 'temporal', 'airport_combined_nonrunway_temporal.csv')
nr = pd.read_csv(nr, dtype={'facility_id': str}, skipinitialspace=True, index_col=False)
temp = pd.concat((temp, nr))
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
temp.to_csv(os.path.join(work_path, 'qa', 'temporal_airport_check.csv'), index=False, 
  columns=['facility_id','src_id','qflag','sum'])
fac_temp, src_temp = (len(temp[['facility_id','src_id']].drop_duplicates('facility_id')), len(temp))
print('Temporal facilities: %s   Temporal sources: %s' %(fac_temp, src_temp))
qa_counts += [str(fac_temp), str(src_temp)]

# Read in the point emissions helper file
emis = os.path.join(work_path, 'xwalk', 'airport_combined_srcid_emis.csv')
emis = pd.read_csv(emis, dtype={'facility_id': str, 'state': str}, skipinitialspace=True, 
  usecols=['state','facility_id','smoke_name','ann_value']) 
emis_polls = emis[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
states = emis[['facility_id','state']].drop_duplicates('state')
emis = emis.groupby('facility_id', as_index=False).sum()
fac_emis = len(emis.drop_duplicates('facility_id')) 
print('Emis facilities: %s  ' %fac_emis)
qa_counts += [str(fac_emis),]

invtable = get_invtable(invtable)
polls = list(invtable['poll'].drop_duplicates())

# Read in all of the inventories from the inventory list
# Select only the sources that have the selected HAPs and are airports
inv_list = get_inv_list()
inv = get_inv(inv_list, invtable)

# Do a pollutant emissions comparison
inv_polls = inv[['smoke_name','ann_value']].groupby('smoke_name', as_index=False).sum()
polls = pd.merge(emis_polls, inv_polls, on='smoke_name', how='outer', suffixes=['_aer','_inv'])
qa_out = os.path.join(work_path, 'qa', 'airport_poll_qa.csv')
polls.to_csv(qa_out, index=False)

# Count the number of unique facilities in the inventories
inv_facs = len(inv.drop_duplicates('facility_id'))
print('Inv facilities: %s ' %inv_facs)
qa_counts.append(str(inv_facs))
# Write the facility and source counts to the "counts"  QA file
count_out = os.path.join(work_path, 'qa', 'airport_counts.csv')
with open(count_out, 'w') as f:
    h = ['lineloc_fac','lineloc_src','nrloc_fac','nrloc_src','loc_fac','loc_src','nrparam_fac',
      'nrparam_src','lineparam_fac','lineparam_src','param_fac','param_src','temp_fac','temp_src','emis_fac','emis_inv']
    f.write('%s\n%s\n' %(','.join(h), ','.join(qa_counts)))

# Use the xwalk file to go from inventory sources to AERMOD sources
inv_emis = inv[['facility_id','ann_value']].groupby('facility_id', as_index=False).sum()
# Write an emissions comparison for the sources between the inventory and helper files
df = pd.merge(inv_emis, emis, on='facility_id', how='outer', suffixes=['_inv','_aer'])
qa_out = os.path.join(work_path, 'qa', 'airport_emis_qa.csv')
df.to_csv(qa_out, index=False, columns=['facility_id','ann_value_inv','ann_value_aer'])

# AERMOD Source crosscheck - verify that all AERMOD sources are in all helper files
xcheck = pd.DataFrame()
for key,df in {'locs': locs, 'params': params, 'temp': temp}.items():
    df = df[['facility_id','src_id']].drop_duplicates()
    df[key] = 'Y'
    if xcheck.empty:
        xcheck = df.copy()
    else:
        xcheck = pd.merge(xcheck, df, on=['facility_id','src_id'], how='outer')
xcheck = xcheck.fillna('N')
qa_out = os.path.join(work_path, 'qa', 'airport_srcid_qa.csv')
xcheck.to_csv(qa_out, index=False)

