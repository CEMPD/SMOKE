#!/usr/bin/env python
# Calculate the onroad hourly profiles by county for AERMOD

import os.path
import pandas as pd

def calc_profile(run_group, inpath):
    '''
    Calculate the profiles for the rungroups
    '''
    fn = os.path.join(inpath,'%s_localtime.csv' %run_group)
    df = pd.read_csv(fn, usecols=['fips','date','hour','emis'])
    df['year'] = df['date'].str[-2:]
    df['month'] = df['date'].str[:2].astype(int)
    df['day'] = df['date'].str[3:5].astype(int)
    df = df[['fips','year','month','day','hour','emis']].copy()
    df_sum = df[['fips','emis']].groupby('fips', as_index=False).sum()
    df = pd.merge(df, df_sum, on='fips', how='left', suffixes=['','_sum'])
    df['emis'] = df['emis'] / df['emis_sum']
    df.insert(0, 'rungroup', run_group)
    df['fips'] = df['fips'].astype(str).str.zfill(5)
    df.rename(columns={'emis': 'scalar'}, inplace=True)
    return df[['rungroup','fips','year','month','day','hour','scalar']]

def write_state_profile(state_df, outpath):
    '''
    Write the profile by state
    '''
    run_group = state_df['rungroup'].values[0]
    state = state_df['fips'].values[0][:2]
    state_df['scalar'] = state_df['scalar'].fillna(0)
    fn = os.path.join(outpath,'%s_%s_hourly.csv' %(run_group, state))
    state_df.to_csv(fn, index=False, float_format='%.12g')

def main():
    # Select the input path, output path, and run groups for the onroad temporal profiles
    inpath = '/work/EMIS/users/bte/WO150.8_2014v2/aermod_helper/2014_onroad_diurnal'
    outpath = '/work/EMIS/users/bte/WO150.8_2014v2/aermod_helper/2014_onroad_diurnal/out'
    run_groups = ['LDOFF12','HOTEL4','HDON4','LDON4','HDOFF12']
    for run_group in run_groups:
        print(run_group)
        prof = calc_profile(run_group, inpath)
        states = list(prof['fips'].str[:2].drop_duplicates())
        for state in states:
            state_df = prof[prof['fips'].str[:2] == state].copy()
            write_state_profile(state_df, outpath)

if __name__ == '__main__':
    main()
