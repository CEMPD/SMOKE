from __future__ import division
from __future__ import print_function
from builtins import object
from past.utils import old_div
import os.path
from datetime import datetime, timedelta
import pandas as pd

class HourlyReports(object):
    '''
    Import the hourly reports typically for EGU sources
    '''
    def __init__(self, sector, case, year, report_path, st_fips=''):
        self.aermod = pd.DataFrame()
        self.header = []
        today = datetime.strptime('0101%s' %year, '%d%m%Y')
        while today.year == int(year):
            fpath = self._name_report(sector, case, report_path, today)
            if os.path.exists(fpath):
                day_df = self._load_report(fpath, today, st_fips)
                self.aermod = pd.concat((self.aermod, day_df))
            else:
                print('WARNING: No report file: %s' %fpath)
            today += timedelta(days=1)
        self.egu_units = pd.DataFrame()
        if not self.aermod.empty:
            self.aermod = self._calc_scalar(self.aermod)
            self.egu_units = self.aermod[['year','unit_id']].copy().drop_duplicates()
            self.egu_units['uniq'] = 'Y'
            self.egu_units = self.egu_units[['unit_id','uniq']].copy()

    def _name_report(self, sector, case, report_path, today):
        fname = 'rep_ptegu_%s_%s_hour.txt' %(case, today.strftime('%Y%m%d'))
        return os.path.join(report_path, fname)
    
    def strip(self, x):
        return x.strip()

    def _parse_header(self, fpath):
        '''
        Read in the EGU hourly report. Format may change.
        '''
        col_dict = {'plant id': 'facility_id', 'char 1': 'unit_id'} #, 'stk tmp': 'stktemp', 'stk vel': 'stkvel'}
        with open(fpath) as f:
            ln = 0
            for line in f:
                if line.startswith('#'):
                    ln += 1
                    if 'date' in line[:7].lower():
                        hdr = [cell.strip().lower() for cell in line[1:].strip().split(';') if cell.strip()]
                        for (i, col) in enumerate(hdr):
                            if col in list(col_dict.keys()):
                                hdr[i] = col_dict[col]
                        self.header = hdr
                        self.header_len = ln
                        break

    def _load_report(self, fpath, today, st_fips):
        '''
        Load hourly SMOKE report
        '''
        if not self.header:
            self._parse_header(fpath)
        df = pd.read_csv(fpath, sep=';', skiprows=self.header_len, skipinitialspace=True, names=self.header,
          usecols=['hour','facility_id','unit_id','co','voc_inv','pm2_5','region'], 
          dtype={'facility_id': str, 'region': str, 'hour': str, 
          'unit_id': str}, converters={'facility_id': self.strip, 'unit_id': self.strip})
        df = df.groupby(['facility_id','unit_id','region','hour'], as_index=False).sum()
        if 'voc' in df.columns and 'voc_inv' not in df.columns:
            df.rename(columns={'voc': 'voc_inv'}, inplace=True)
        if st_fips:
            df['region'] = df['region'].str.zfill(5)
            df = df[df['region'].str.startswith(st_fips)]
        df['year'] = today.year
        df['month'] = today.month
        df['day'] = today.day
        try:
            df['emis'] = df['co'] + df['pm2_5'] + df['voc_inv']
        except KeyError:
            raise KeyError("Missing CO, PM2_5, or VOC_INV in hourly PTEGU report")
        return df[['facility_id','unit_id','year','month','day','hour','emis']].copy()

    def _calc_scalar(self, prof):
        '''
        Calculate the hourly scalar from the SMOKE reports
        '''
        sum_df = prof[['facility_id','unit_id','emis']].copy().groupby(['facility_id','unit_id'],
          as_index=False).sum()
        sum_df.rename(columns={'emis': 'sum'}, inplace=True)
        prof = pd.merge(prof, sum_df, on=['facility_id','unit_id'], how='left')
        prof['hour_factor'] = (old_div(prof['emis'],prof['sum'])) * 8760.
        return prof[['facility_id','unit_id','year','month','day','hour','hour_factor']].copy()

