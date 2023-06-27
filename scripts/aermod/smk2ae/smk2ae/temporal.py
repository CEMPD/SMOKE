from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import numpy as np
import pandas as pd

class Temporal(object):
    '''
    Import the temporal crossreference file and merge against the temporal profiles
    '''
    def __init__(self, tref=None, hour=None, week=None, month=None, scc_list=[], 
      fips_list=[], fac_list=[]):
        self.use_daily = False
        if tref:
            self.xref = self._get_tref(tref, scc_list, fips_list, fac_list)
            self._hour = self._get_hour_pro(hour)
            self._week = self._get_weekly_pro(week)
            self._month = self._get_monthly_pro(month)
            if self.use_daily:
                self._daily = self._get_daily_pro()
            self.profs = self.gen_temp_prof()

    def _get_tref(self, tref, scc_list, fips_list, fac_list):
        '''
        Read in the temporal x-ref

        * Assume that profiles only match on SCC and/or FIPS; Does not match on pollutant
         or facility or unit or release point or process ID.
        This is a simplification that currently works for ptnonipm and ptegu.
        '''
        names = ['scc','region_cd','facility_id','unit_id','rel_point_id','process_id','poll',
            'type','code','comment']
        usecols = ['facility_id','scc','region_cd','type','code']
        dtype={'region_cd': str, 'scc': str, 'code': str,'facility_id': str}
        xref = pd.read_csv(tref, comment='#', names=names, usecols=usecols, dtype=dtype)
        xref.fillna('', inplace=True)
        key_cols = ['scc','region_cd','facility_id']
        # Replace 0s with null strings
        for col in key_cols:
            xref.loc[xref[col].str.zfill(15) == '000000000000000', col] = ''
        # Fill region_cd to 5 characters using 0s
        xref.loc[xref['region_cd'] != '', 'region_cd'] = xref.loc[xref['region_cd'] != '', 
          'region_cd'].str[-5:].str.zfill(5)
        # Narrow down the x-refs to the ones that are in the inventory. This saves processing time.
        xref.set_index(key_cols, inplace=True)
        xref.sort_index(inplace=True)
        out_xref = pd.DataFrame()
        # Iterate over a matching hierarchy for the temporal x-ref. Hierarchy goes from most specific
        #  and highest rank to least specific and lowest rank.
        for idx_keys in ((scc_list, fips_list, fac_list), (scc_list, fips_list, ''), 
          (scc_list, '', fac_list), ('', fips_list, fac_list), (scc_list, '', ''),
          ('', fips_list, ''), ('', '', fac_list)):
            try:
                xref_idx = xref.loc[idx_keys]
            except KeyError:
                pass
            else:
                out_xref = pd.concat((out_xref, xref_idx))
        out_xref = out_xref.reset_index().drop_duplicates(key_cols.extend('type'))
        # Weekly profiles are not used if there are daily profiles
        if 'DAILY' in list(out_xref['type'].drop_duplicates()):
            self.use_daily = True
            out_xref = out_xref[out_xref['type'] != 'WEEKLY'].copy()
        else:
            self.use_daily = False
        return out_xref

    def _get_hour_pro(self, hour):
        '''
        Get the hourly temporal profile and cross-reference it to fips/scc/facility
        Not worried about pollutants in these x-refs because no temp profiles are HAP only
        '''
        hours = ['hr%0.2d' %hr for hr in range(1,25)]
        names = ['code',]+hours+['comment',]
        usecols = ['code',]+hours
        dtype = {'code': str}
        for hr in hours:
            dtype[hr] = float 
        # The hourly profile sometimes has column names in the file. If it does, skip them.
        with open(hour) as f:
            skip_lines = 0
            for l in f:
                if l.startswith('#') or l.lower().startswith('profile_id'):
                    skip_lines += 1
                else:
                    break
        df = pd.read_csv(hour, skiprows=skip_lines, names=names, usecols=usecols, dtype=dtype,
          index_col=False)
        # There are multiple types of diurnal profiles that go from general (ALLDAY) to specific (MONDAY...SUNDAY)
        types = ['ALLDAY','WEEKDAY','WEEKEND','MONDAY','TUESDAY','WEDNESDAY','THURSDAY','FRIDAY',
            'SATURDAY','SUNDAY']
        df = pd.merge(self.xref[self.xref['type'].isin(types)], df, on='code', how='left')
        return df[usecols].copy().drop_duplicates()

    def _get_daily_pro(self):
        '''
        Get the daily temporal profile and cross-reference it to fips/scc/facility
        Not worried about pollutants in these x-refs because no temp profiles are HAP only
        '''
        from os import environ
        names = ['code','month']+['d%s' %day for day in range(1,32)]
        # Pull the TPRO_DAILY from the environment variable. Not a great way to do this..
        df = pd.read_csv(environ['ATPRO_DAILY'], comment='#', names=names, dtype={'code': str},
          index_col=False)
        df = pd.merge(self.xref[self.xref['type'] == 'DAILY'], df, on='code', how='left').drop_duplicates()
        return df[names].copy().drop_duplicates()

    def _get_weekly_pro(self, week):
        '''
        Get the weekly temporal profile and cross-reference it to fips/scc/facility
        Not worried about pollutants in these x-refs because no temp profiles are HAP only
        '''
        names = ['code','mon','tue','wed','thu','fri','sat','sun','comments']
        usecols = ['code','mon','tue','wed','thu','fri','sat','sun']
        df = pd.read_csv(week, comment='#', names=names, usecols=usecols, dtype={'code': str})
        df = pd.merge(self.xref[self.xref['type']=='WEEKLY'], df, on='code', how='left')
        return df[usecols].copy().drop_duplicates()

    def _get_monthly_pro(self, month):
        '''
        Get the monthly temporal profile and cross-reference it to fips/scc/facility
        Not worried about pollutants in these x-refs because no temp profiles are HAP only
        '''
        names = ['code','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','cmt']
        usecols = ['code','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        df = pd.read_csv(month, comment='#', names=names, usecols=usecols, dtype={'code': str})
        df = pd.merge(self.xref[self.xref['type']=='MONTHLY'], df, on='code', how='left')
        return df[usecols].copy().drop_duplicates()

    def _calc_daily_prof(self, row):
        '''
        Calculate the hourly temporal profiles using the month->day daily profiles
        This is for sectors that use daily profiles rather than monthly/weekly.
        '''
        from calendar import monthrange
        month_prof = self._month[self._month['code'] == row['MONTHLY']]
        daily_prof = self._daily[self._daily['code'] == row['DAILY']]
        hour_prof = self._hour[self._hour['code'] == row['ALLDAY']]
        month_cols = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
        hr_table = []
        for mon_num, mon in enumerate(month_cols):
            mon_val = float(month_prof[mon])
            day_month_df = daily_prof[daily_prof['month'] == (mon_num + 1)].copy()
            if day_month_df.empty:
                print('WARNING: Missing daily profile')
            else:
                # Assume not a leap year (2014). Change to year later.
                for day in range(1,monthrange(2014, mon_num + 1)[1]+1):
                    try:
                        day_val = float(day_month_df['d%s' %day])
                    except ValueError:
                        pass
                    else:
                        for hr in range(1,25):
                            factor = mon_val*day_val*float(hour_prof['hr%0.2d' %hr])
                            hr_table.append(list(row) + [mon_num+1, day, hr, factor])
        if hr_table:
            return pd.DataFrame(hr_table, columns=list(row.index) + ['month','day','hour','factor'])
        else:
            return pd.DataFrame() 

    def fill_missing(self, keys):
        '''
        Fill missing profile codes by moving to SCC for entries with scc/region combos
        Looks for unmatched xrefs that have scc xrefs
        '''
        for key in keys:
            idx = (self.xref['region_cd'] != '') & (self.xref[key].isnull())
            miss_len = len(self.xref[idx])
            if miss_len > 0:
                print('NOTE: Filling %s %s county profiles with SCC default profiles' %(miss_len, key))
                fill_profs = self.xref.loc[self.xref['region_cd'] == '', ['scc',key]].drop_duplicates('scc')
                self.xref = pd.merge(self.xref, fill_profs, on='scc', how='left', suffixes=['','_def'])
                self.xref.loc[idx, key] = self.xref.loc[idx, '%s_def' %key]
                self.xref.drop('%s_def' %key, axis=1, inplace=True)

    def gen_temp_prof(self):
        '''
        Calculate the AERMOD ready temporal profiles from the temporal XREF
        This iterates through all selected cross references to retrieve a complete set of scalars
          for that profile.
        '''
        keys = ['facility_id','scc','region_cd']
        self.xref = pd.pivot_table(self.xref, index=keys,
          columns='type', values='code', aggfunc=lambda x: ' '.join(x)).reset_index()
        if self.xref.empty:
            raise ValueError('No valid temporal cross references found')
        vals = [col for col in self.xref.columns if col not in keys]
        self.fill_missing(vals)
        uniq_profs = self.xref[vals].copy().drop_duplicates()
        aermod_profs = pd.DataFrame()
        # Iterate over all of the unique profiles as a combination of x-refs
        for idx, row in uniq_profs.iterrows():
            if self.use_daily:
                prof = self._calc_daily_prof(row)
            else:
                prof = self._calc_temp_prof(row)
            if not prof.empty:
                aermod_profs = pd.concat((aermod_profs, prof))
                if len(prof.columns) > len(aermod_profs.columns):
                    aermod_profs.columns = prof.columns
        if self.use_daily:
            keys = [col for col in aermod_profs.columns if col != 'factor']
        else:
            keys = [col for col in aermod_profs.columns if not col.startswith('Scalar')]
        vals = list(uniq_profs.columns)
        if aermod_profs.empty:
            print('WARNING: No temporal profiles found')
        else:
            # Merge the calculated profiles onto the xrefs
            aermod_profs.drop_duplicates(keys, inplace=True)
            if self.use_daily:
                aermod_profs[['month','day','hour']] = \
                  aermod_profs[['month','day','hour']].astype('i')
                aermod_profs = pd.merge(self.xref, aermod_profs, on=vals, how='left')
            else:
                aermod_profs = pd.merge(self.xref, aermod_profs, on=vals, how='left')
        # Keep track of columns because not all of the profiles have the same number of scalars
        columns = [col for col in aermod_profs.columns if col not in vals]
        return aermod_profs[columns]

    def _calc_month(self, month_prof, cols, col_names, sum_mon):
        '''
        Profile flat across day of week and hour of day, only varies by month
        '''
        mons = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        col_names += ['Scalar%s' %x for x in range(1,13)]
        cols.append('MONTH')
        mon_vals = [(float(month_prof[mon])/sum_mon) for mon in mons]
        cols.extend(mon_vals)
        return cols, col_names

    def _calc_hrofdy(self, hour_prof, cols, col_names, hours):
        '''
        Profile flat across month of year and day of week, only varies by hour of day
        '''
        col_names += ['Scalar%s' %x for x in range(1,25)]
        cols.append('HROFDY')
        sum_hour = float(hour_prof[hours].sum(axis=1))
        hr_vals = [(float(hour_prof['hr%0.2d' %hr])/sum_hour) for hr in range(1,25)]
        cols.extend(hr_vals)
        return cols, col_names

    def _calc_temp_prof(self, row):
        '''
        Calculate the AERMOD ready temporal profile for each region/SCC in temporal XREF
        This is the calculation of a unique profile using each unique combination of xref keys.
        '''
        row = row[row.notnull()].copy()
        month_prof = self._month[self._month['code'] == row['MONTHLY']].drop_duplicates()
        mons = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        sum_mon = float(month_prof[mons].sum(axis=1))
        week_prof = self._week[self._week['code'] == row['WEEKLY']].drop_duplicates()
        try:
            hour_prof = self._hour[self._hour['code'] == row['ALLDAY']].drop_duplicates()
        except KeyError:
            hour_prof = self._hour[self._hour['code'] == row['WEEKDAY']].drop_duplicates()
            row['ALLDAY'] = ''
        hours = ['hr%0.2d' %hr for hr in range(1,25)]
        col_names = list(row.index) + ['qflag',]
        cols = list(row.values)
        # If a profile is flat over a week and diurnal, then it only varies by MONTH so it is monthly
        if str(row['WEEKLY']) == '7' and str(row['ALLDAY']) == '24':
            cols, col_names = self._calc_month(month_prof, cols, col_names, sum_mon)
        # If a profile is flat over the months and a week, then it only varies by hour so it is hourly
        elif str(row['MONTHLY']) == '262' and str(row['WEEKLY']) == '7':
            cols, col_names = self._calc_hrofdy(hour_prof, cols, col_names, hours)
        # Otherwise it varies across multiple profile types
        else:
            '''
            Otherwise assume variation by month of year, day of week, and hour of day
            '''
            # Try to get the weekday average.
            try:
                wkday_avg = sum(week_prof[['mon','tue','wed','thu','fri']].values[0])/5.
            except KeyError:
                raise KeyError(week_prof)
            if (week_prof['tue'].values[0]/wkday_avg) == 1:
                '''
                 Assume that day of week profiles only vary between weekdays, Saturday, and Sunday
                   rather than every individual day of the week.
                '''
                day_dict = {'tue': 'TUESDAY', 'sat': 'SATURDAY', 'sun': 'SUNDAY'}
                day_list = ['tue','sat','sun']
                # This profile will end up with 864 scalars
                col_names += ['Scalar%s' %x for x in range(1,865)]
                cols.append('MHRDOW')
            else:
                '''
                Assume variation by weekday
                '''
                day_dict = {'mon': 'MONDAY', 'tue': 'TUESDAY', 'wed': 'WEDNESDAY', 
                  'thu': 'THURSDAY', 'fri': 'FRIDAY', 'sat': 'SATURDAY', 'sun': 'SUNDAY'}
                day_list = ['mon','tue','wed','thu','fri','sat','sun']
                col_names += ['Scalar%s' %x for x in range(1,2017)]
                cols.append('MHRDOW7')
            sum_day = float(week_prof[['mon','tue','wed','thu','fri','sat','sun']].sum(axis=1))
            # Check for diurnal profiles by day of weekday, weekday, weekend, otherwise go with
            #  the all day profile
            for day in day_list:
                day_key = day_dict[day]
                day_val = float(week_prof[day])/sum_day
                try:
                    code = row[day_key]
                except KeyError:
                    if day in ('mon','tue','wed','thu','fri'):
                        day_key = 'WEEKDAY'
                    else:
                        day_key = 'WEEKEND'
                    try:
                        code = row[day_key]
                    except KeyError:
                        hour_day = hour_prof.copy()
                    else:
                        hour_day = self._hour[self._hour['code'] == row[day_key]].drop_duplicates()
                else:
                    hour_day = self._hour[self._hour['code'] == row[day_key]].drop_duplicates()
                sum_hour = float(hour_day[hours].sum(axis=1))
                for mon in mons:
                    mon_val = float(month_prof[mon])/sum_mon
                    hr_vals = [(float(hour_day['hr%0.2d' %hr])/sum_hour)*mon_val*day_val for hr in range(1,25)]
                    cols.extend(hr_vals)
        return pd.DataFrame([cols], columns=col_names)

def match_temporal(df, xref, val_cols, hierarchy, use_daily=False):
    '''
    Match data to temporal based on xref keys using a hierarchy of keys
    * This is an incomplete hierarchy based on region, scc, and facility only.
        This hierarchy should be fine for HAPS sources in the 2014v2 NEI.
    '''
    key_cols = [col for col in df.columns if col not in val_cols]
    matched_df = pd.DataFrame(columns=[val_cols[0],])
    # Iterate over the index groups in the passed hierarchy
    for match_cols in hierarchy:
        df.set_index(match_cols, inplace=True)
        h_xref = xref.copy()
        # Try to matc those columns, otherwise move on to the next level
        for col in key_cols:
            if col in h_xref.columns:
                if col in match_cols:
                    h_xref = h_xref[h_xref[col] != ''].copy()
                else:
                    h_xref = h_xref[h_xref[col] == ''].copy()
        if h_xref.empty:
            df.reset_index(inplace=True)
        else:
            h_xref = h_xref[match_cols+val_cols].copy().set_index(match_cols)
            df = pd.merge(df, h_xref, left_index=True, right_index=True, how='left')
            if use_daily:
                uniq_cols = key_cols + ['month','day','hour']
            else:
                uniq_cols = key_cols
            df = df.reset_index().drop_duplicates(uniq_cols)
            matched_df = pd.concat((matched_df, df[df[val_cols[0]].notnull()]))
            df = df[df[val_cols[0]].isnull()].copy()
            df.drop(val_cols, axis=1, inplace=True)
            if df.empty:
                break
    df = pd.concat((matched_df, df))
    # Check for any source that weren't matched to the tref
    unmatched = df[df['qflag'].isnull()].copy()
    if len(unmatched) > 0:
        print(unmatched[:20])
        raise ValueError('Unmatched sources to temporal profiles')
    return df

def fill_default(df, default_prof):
    '''
    Fill missing temporal profiles with a default profile
    '''
    null_df = df[df[default_prof.index[0]].isnull()].copy()
    df = df[~ df[default_prof.index[0]].isnull()].copy()
    id_cols = [col for col in null_df.columns if col not in default_prof.index]
    null_df = pd.merge(null_df[id_cols], pd.DataFrame(data=[default_prof.values] * len(null_df),
      columns=default_prof.index, index=null_df.index), left_index=True, right_index=True)
    null_df.drop_duplicates(id_cols, inplace=True)
    return pd.concat((df, null_df))



