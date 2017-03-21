import numpy as np
import pandas as pd

class Temporal(object):
    '''
    Import the temporal crossreference file and merge against the temporal profiles
    '''
    def __init__(self, tref, hour, week, month, scc_list, fips_list, fac_list=[]):
        self.xref = self._get_tref(tref, scc_list, fips_list, fac_list)
        self._hour = self._get_hour_pro(hour)
        self._week = self._get_weekly_pro(week)
        self._month = self._get_monthly_pro(month)
        if self.use_daily:
            self._daily = self._get_daily_pro()
        self.aermod = self.gen_temp_prof()

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
        dtype={'region_cd': '|S6', 'scc': '|S10', 'code': '|S8','facility_id': '|S15'}
        xref = pd.read_csv(tref, comment='#', names=names, usecols=usecols, dtype=dtype)
        key_cols = ['scc','region_cd','facility_id']
        for col in key_cols:
            xref.ix[xref[col].str.zfill(15) == '000000000000000', col] = ''
        xref.ix[xref['region_cd'] != '', 'region_cd'] = xref.ix[xref['region_cd'] != '', 
          'region_cd'].str.zfill(6)
        # Narrow down the x-refs to the ones that are in the inventory. This saves processing time.
        xref.set_index(key_cols, inplace=True)
        xref.sort_index(inplace=True)
        out_xref = pd.DataFrame()
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
        hours = ['hr%0.2d' %hr for hr in xrange(1,25)]
        names = ['code',]+hours+['comment',]
        usecols = ['code',]+hours
        dtype = {'code': '|S8'}
        for hr in hours:
            dtype[hr] = 'f'
        try:
            df = pd.read_csv(hour, comment='#', names=names, usecols=usecols, dtype=dtype,
              index_col=False)
        except AttributeError:
            raise AttributeError, 'Please comment column headers in %s.\n \
              If error persists check hourly factors for non-floats.' %hour
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
        df = pd.read_csv(environ['ATPRO_DAILY'], comment='#', names=names, dtype={'code': '|S8'},
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
        df = pd.read_csv(week, comment='#', names=names, usecols=usecols, dtype={'code': '|S8'})
        df = pd.merge(self.xref[self.xref['type']=='WEEKLY'], df, on='code', how='left')
        return df[usecols].copy().drop_duplicates()

    def _get_monthly_pro(self, month):
        '''
        Get the monthly temporal profile and cross-reference it to fips/scc/facility
        Not worried about pollutants in these x-refs because no temp profiles are HAP only
        '''
        names = ['code','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','cmt']
        usecols = ['code','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        df = pd.read_csv(month, comment='#', names=names, usecols=usecols, dtype={'code': '|S8'})
        df = pd.merge(self.xref[self.xref['type']=='MONTHLY'], df, on='code', how='left')
        return df[usecols].copy().drop_duplicates()

    def _calc_daily_prof(self, row):
        '''
        Calculate the hourly temporal profiles using the month->day daily profiles
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
                print 'WARNING: Missing daily profile'
            else:
                # Assume not a leap year (2014). Change to year later.
                for day in xrange(1,monthrange(2014, mon_num + 1)[1]+1):
                    try:
                        day_val = float(day_month_df['d%s' %day])
                    except ValueError:
                        pass
                    else:
                        for hr in xrange(1,25):
                            factor = mon_val*day_val*float(hour_prof['hr%0.2d' %hr])
                            hr_table.append(list(row) + [mon_num+1, day, hr, factor])
        if hr_table:
            return pd.DataFrame(hr_table, columns=list(row.index) + ['month','day','hour','factor'])
        else:
            return pd.DataFrame() 

    def gen_temp_prof(self):
        '''
        Calculate the AERMOD ready temporal profiles from the temporal XREF
        '''
        keys = ['facility_id','scc','region_cd']
        self.xref = pd.pivot_table(self.xref, index=keys,
          columns='type', values='code', aggfunc=lambda x: ' '.join(x)).reset_index()
        vals = [col for col in self.xref.columns if col not in keys]
        uniq_profs = self.xref[vals].copy().drop_duplicates()
        aermod = pd.DataFrame()
        for idx, row in uniq_profs.iterrows():
            if self.use_daily:
                prof = self._calc_daily_prof(row)
            else:
                prof = self._calc_temp_prof(row)
            if not prof.empty:
                aermod = pd.concat((aermod, prof))
                if len(prof.columns) > len(aermod.columns):
                    aermod.columns = prof.columns
        if self.use_daily:
            keys = [col for col in aermod.columns if col != 'factor']
        else:
            keys = [col for col in aermod.columns if not col.startswith('Scalar')]
        aermod.drop_duplicates(keys, inplace=True)
        vals = list(uniq_profs.columns)
        if self.use_daily:
            aermod[['month','day','hour']] = aermod[['month','day','hour']].astype('i')
            aermod = pd.merge(self.xref, aermod, on=vals, how='left')
        else:
            aermod = pd.merge(self.xref, aermod, on=vals, how='left')
        columns = [col for col in aermod.columns if col not in vals]
        return aermod[columns]

    def _calc_month(self, month_prof, cols, col_names, avg_mon):
        '''
        Profile flat across day of week and hour of day, only varies by month
        '''
        mons = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        col_names += ['Scalar%s' %x for x in xrange(1,13)]
        cols.append('MONTH')
        mon_vals = [(float(month_prof[mon])/avg_mon) for mon in mons]
        cols.extend(mon_vals)
        return cols, col_names

    def _calc_hrofdy(self, hour_prof, cols, col_names, hours):
        '''
        Profile flat across month of year and day of week, only varies by hour of day
        '''
        col_names += ['Scalar%s' %x for x in xrange(1,25)]
        cols.append('HROFDY')
        avg_hour = float(hour_prof[hours].mean(axis=1))
        hr_vals = [(float(hour_prof['hr%0.2d' %hr]) / avg_hour) for hr in xrange(1,25)]
        cols.extend(hr_vals)
        return cols, col_names

    def _calc_temp_prof(self, row):
        '''
        Calculate the AERMOD ready temporal profile for each region/SCC in temporal XREF
        '''
        row = row[row.notnull()].copy()
        month_prof = self._month[self._month['code'] == row['MONTHLY']].drop_duplicates()
        mons = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
        avg_mon = float(month_prof[mons].mean(axis=1))
        week_prof = self._week[self._week['code'] == row['WEEKLY']].drop_duplicates()
        try:
            hour_prof = self._hour[self._hour['code'] == row['ALLDAY']].drop_duplicates()
        except KeyError:
            hour_prof = self._hour[self._hour['code'] == row['WEEKDAY']].drop_duplicates()
            row['ALLDAY'] = ''
        hours = ['hr%0.2d' %hr for hr in xrange(1,25)]
        col_names = list(row.index) + ['qflag',]
        cols = list(row.values)
        if str(row['WEEKLY']) == '7' and str(row['ALLDAY']) == '24':
            cols, col_names = self._calc_month(month_prof, cols, col_names, avg_mon)
        elif str(row['MONTHLY']) == '262' and str(row['WEEKLY']) == '7':
            cols, col_names = self._calc_hrofdy(hour_prof, cols, col_names, hours)
        else:
            '''
            Otherwise assume variation by month of year, day of week, and hour of day
            '''
            try:
                wkday_avg = sum(week_prof[['mon','tue','wed','thu','fri']].values[0]) / 5.
            except KeyError:
                raise KeyError, week_prof
            if (week_prof['tue'].values[0] / wkday_avg) == 1:
                '''
                 Assume that day of week profiles only vary between weekdays, Saturday, and Sunday
                   rather than every individual day of the week.
                '''
                day_dict = {'tue': 'TUESDAY', 'sat': 'SATURDAY', 'sun': 'SUNDAY'}
                day_list = ['tue','sat','sun']
                col_names += ['Scalar%s' %x for x in xrange(1,865)]
                cols.append('MHRDOW')
            else:
                '''
                Assume variation by weekday
                '''
                day_dict = {'mon': 'MONDAY', 'tue': 'TUESDAY', 'wed': 'WEDNESDAY', 
                  'thu': 'THURSDAY', 'fri': 'FRIDAY', 'sat': 'SATURDAY', 'sun': 'SUNDAY'}
                day_list = ['mon','tue','wed','thu','fri','sat','sun']
                col_names += ['Scalar%s' %x for x in xrange(1,2017)]
                cols.append('MHRDOW7')
            avg_day = float(week_prof[['mon','tue','wed','thu','fri','sat','sun']].mean(axis=1))
            for day in day_list:
                day_key = day_dict[day]
                day_val = float(week_prof[day])/avg_day
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
                avg_hour = float(hour_day[hours].mean(axis=1))
                for mon in mons:
                    mon_val = float(month_prof[mon])/avg_mon
                    hr_vals = [(float(hour_day['hr%0.2d' %hr])/avg_hour)*mon_val*day_val for hr in xrange(1,25)]
                    cols.extend(hr_vals)
        return pd.DataFrame([cols], columns=col_names)

def match_temporal(df, xref, val_cols, hierarchy, use_daily=False):
    '''
    Match data to temporal based on xref keys using a hierarchy of keys
    * This is an incomplete hierarchy based on region, scc, and facility only.
        This hierarchy should be fine for HAPS sources in the current NEI.
    '''
    key_cols = [col for col in df.columns if col not in val_cols]
    matched_df = pd.DataFrame(columns=[val_cols[0],])
    for match_cols in hierarchy:
        df.set_index(match_cols, inplace=True)
        h_xref = xref.copy()
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
    return pd.concat((matched_df, df))

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



