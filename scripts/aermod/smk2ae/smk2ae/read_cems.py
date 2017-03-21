import pandas as pd
cems = '/sol/work/EMIS/users/bte/WO150.1_2014plat/AERMOD/qa/temporal/cems_2014/HOUR_UNIT_2014_revised.txt'
inv = '/sol/work/EMIS/em_v7/2014_NATA/2014_NATA_v1/inputs/ptegu/ptegu_2014NEIv1_final_POINT_20160825_09sep2016_v0.csv'
xwalk = '/sol/work/EMIS/users/bte/WO150.1_2014plat/AERMOD/xwalk/point_srcid_xwalk.csv'

cems = pd.read_csv(cems, names=['oris','boiler','date','hour','nox','so2','noxr','optime','gload','sload','heat','heatm','so2m','noxm','noxrm','flow'],
  usecols=['oris','boiler','date','hour','heat'], dtype={'oris': '|S6', 'boiler': '|S6', 'date': '|S6'}, na_values=['-9',-9])
cems.fillna(0, inplace=True)

class CEMS:
    '''
    Import the CEMS as
    '''
/sol/work/EMIS/em_v7/2014platform/work/CEMS/revised/SMOKE/H UR_UNIT_2013_12_31dec.txt
class HourlyReports:
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
                if self.aermod.empty:
                    self.aermod = day_df.copy()
                else:
                    self.aermod = pd.concat((self.aermod, day_df))
            else:
                print 'WARNING: No report file: %s' %fpath
            today += timedelta(days=1)
        if self.aermod.empty:
            self.rel_point_list = []
        else:
            self.aermod = self._calc_scalar(self.aermod)
            self.rel_point_list = list(self.aermod['rel_point_id'].drop_duplicates())


    def _calc_scalar(self, prof):
        '''
        Calculate the hourly scalar from the SMOKE reports
        '''
        sum_df = prof[['facility_id','rel_point_id','emis']].copy().groupby(['facility_id',
            'rel_point_id'],as_index=False).sum() 
        sum_df.rename(columns={'emis': 'sum'}, inplace=True)
        prof = pd.merge(prof, sum_df, on=['facility_id','rel_point_id'], how='left')
        prof['hour_factor'] = (prof['emis']/prof['sum']) * 8760.
        return prof[['facility_id','rel_point_id','year','month','day','hour','hour_factor']].copy()

