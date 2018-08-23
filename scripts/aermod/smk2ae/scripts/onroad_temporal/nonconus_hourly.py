#!/usr/bin/env python
# Calculate the onroad hourly profiles by county for AERMOD
# This script calculates the non-CONUS run group values based on the US values

import os.path
import pandas as pd

def calc_profile(run_group, inpath):
    '''
    Calculate the profiles for the rungroups
    '''
    fn = os.path.join(inpath,'%s_localtime.csv' %run_group)
    df = pd.read_csv(fn, usecols=['date','hour','emis'])
    df = df.groupby(['date','hour'], as_index=False).sum()
    df['year'] = df['date'].str[-2:]
    df['month'] = df['date'].str[:2].astype(int)
    df['day'] = df['date'].str[3:5].astype(int)
    df = df[['year','month','day','hour','emis']].copy()
    df_sum = sum(df['emis'])
    df['emis'] = df['emis'] / df_sum
    df.insert(0, 'rungroup', run_group)
    df.rename(columns={'emis': 'scalar'}, inplace=True)
    df = df[~ ((df['year'] == '13') & (df['hour'].isin((17,18,19))))].copy()
    df.loc[df['year'] == '13', 'year'] = '14'
    return df[['rungroup','year','month','day','hour','scalar']].sort_values(['year','month','day','hour'])

def write_state_profile(state_df, outpath):
    '''
    Write the profile by state
    '''
    run_group = state_df['rungroup'].values[0]
    state = state_df['fips'].values[0][:2]
    state_df['scalar'] = state_df['scalar'].fillna(0)
    fn = os.path.join(outpath,'%s_%s_hourly.csv' %(run_group, state))
    columns = ['rungroup','fips','year','month','day','hour','scalar']
    state_df.to_csv(fn, index=False, float_format='%.12g', columns=columns)

def main():
    inpath = '/work/EMIS/users/bte/WO150.8_2014v2/aermod_helper/2014_onroad_diurnal'
    outpath = '/work/EMIS/users/bte/WO150.8_2014v2/aermod_helper/2014_onroad_diurnal/out'
    run_groups = ['LDOFF12','HOTEL4','HDON4','LDON4','HDOFF12']
    # Define dictionary of states to generate and all counties contained in those states
    states = {'15': [1,3,5,7,9], '02': [13,16,20,50,60,68,70,90,100,105,110,122,130,150,164,170,180,
      185,188,195,198,220,230,240,261,270,275,282,290], '78': [10,20,30], '72': [1,3,5,7,9,11,13,15,
      17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,54,55,57,59,61,63,65,67,69,71,73,75,
      77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127,
      129,131,133,135,137,139,141,143,145,147,149,151,153]}
    domain_groups = {'02': '9AK', '15': '3HI', '72': '3PR', '78': '3PR'}
    for run_group in run_groups:
        print(run_group)
        prof = calc_profile(run_group, inpath)
        for state, counties in states.items():
            st_df = pd.DataFrame(['%s%0.3d' %(state, county) for county in counties], 
              columns=['fips',])
            st_df['rungroup'] = run_group
            st_df = pd.merge(st_df, prof, on='rungroup', how='left')
            if run_group.endswith('4'):
                domain_group = run_group[:-1] + domain_groups[state]
            else:
                domain_group = run_group[:-2] + domain_groups[state]
            st_df['rungroup'] = domain_group
            write_state_profile(st_df, outpath)

if __name__ == '__main__':
    main()
