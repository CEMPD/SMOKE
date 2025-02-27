#!/usr/bin/env python
'''
QA for AERMOD nonrunway helper files
Requires python 3.5, pyproj, and pandas
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

def loc_check(sector, work_path):
    # Location file QA
    locs = os.path.join(work_path, 'locations', '%s_nonrunway_locations.csv' %sector)
    usecols=['facility_id','src_id','longitude','latitude','utm_x','utm_y','utm_zone']
    locs = pd.read_csv(locs, dtype={'facility_id': str}, usecols=usecols, skipinitialspace=True)
    return locs 

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

def get_inv(sector, inv_list, invtable):
    # Read in all of the inventories from the inventory list
    # Select only the sources that have the selected HAPs and are not airports
    polls = list(invtable['poll'].drop_duplicates())
    inv = pd.DataFrame()
    usecols = ['facility_id','latitude','longitude','poll','ann_value']
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
        df = df[(df['poll'].isin(polls)) & (df['ann_value'] > 0)].copy()
        df.drop_duplicates('facility_id', inplace=True)
        inv = pd.concat((inv, df))
    return inv


def main():
    sector = os.environ['SECTOR']
    work_path = os.environ['WORK_PATH'] 
    invtable = os.environ['INVTABLE']
    if not os.path.exists(os.path.join(work_path, 'qa')):
        os.makedirs(os.path.join(work_path, 'qa'))
    # Check the locations file
    locs = loc_check(sector, work_path)
    invtable = get_invtable(invtable)
    # Get the inventory emissions
    inv_list = get_inv_list()
    inv = get_inv(sector, inv_list, invtable)
    inv = pd.merge(locs, inv, on='facility_id', how='left', suffixes=['_aer','_inv'])
    utm = UTM()
    inv['utm_zone_inv'] = inv['longitude_inv'].apply(lambda lon: utm.get_zone(lon))
    inv[['utm_x_inv','utm_y_inv']] = inv[['longitude_inv','latitude_inv','utm_zone_inv']].apply(lambda x: \
      pd.Series(utm.get_coords(x['longitude_inv'],x['latitude_inv'],x['utm_zone_inv'])), axis=1)
    qa_out = os.path.join(work_path, 'qa', '%s_locs_check.csv' %sector)
    inv.to_csv(qa_out, index=False)

if __name__ == '__main__':
    main()
