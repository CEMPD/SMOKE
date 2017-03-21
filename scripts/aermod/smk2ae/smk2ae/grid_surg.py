import os.path
import numpy as np
import pandas as pd

class GridSurg(object):
    '''
    Import the temporal crossreference file and merge against the temporal profiles
    '''
    def __init__(self, gref, srg_path, srg_desc, scc_list):
        self.xref = self._get_gref(gref, scc_list)
        if self.xref.empty:
            raise ValueError, 'No surrogates found for inventory SCCs in xref'
        srg_path = self._parse_srg_path(srg_path)
        self.desc = self._get_srgdesc(srg_desc, srg_path)
        self.surgs = dict()
        self._load_surgs()

    def _parse_srg_path(self, srg_path):
        '''
        Get the directory name of the surrogate profile from the SRGPRO environment variable
        '''
        return os.path.dirname(srg_path)

    def _fix_fips(self, fips):
        '''
        Convert the FIPS into a 5 character zero-padded string
        Assume that all sources are US only
        '''
        try:
            fips = int(fips)
        except ValueError:
            return fips
        else:
            if fips < 100000 and fips > 0:
                return '%0.6d' %fips
            else:
                return ''

    def _get_gref(self, gref, scc_list):
        '''
        Read in the gridding x-ref file (AGREF)
        '''
        names = ['region_cd','scc','code']
        df = pd.read_csv(gref, names=names, comment='#', sep=';', 
            dtype={'region_cd': '|S6', 'scc': '|S10', 'code': '|S8'})
        df['region_cd'] = df['region_cd'].apply(self._fix_fips) 
        df.ix[df['scc'].str.startswith('00'), 'scc'] = df.ix[df['scc'].str.startswith('00'), 
          'scc'].str[2:]
        df['code'] = df['code'].str.split('!').str[0].str.strip()
        return df[df['scc'].isin(scc_list)].copy()

    def _get_srgdesc(self, srg_desc, srg_path):
        '''
        Read in the surrogate description file to get the names of the surrogates
         by surrogate code
        '''
        names = ['country','code','desc','fname']
        try:
            df = pd.read_csv(srg_desc, names=names, comment='#', sep=';', usecols=['code','fname'],
              dtype={'code': '|S3'})
        except:
            df = pd.read_csv(srg_desc, names=names, comment='#', sep=',', usecols=['code','fname'],
              dtype={'code': '|S3'})
        df['fname'] = '%s/' %srg_path + df['fname']
        return df

    def _load_surgs(self):
        '''
        Load the surrogates by code and file name into a surrogate object
         stored in a srg code indexed dictionary
        '''
        xref = pd.merge(self.xref, self.desc, on='code', how='left')
        xref = xref[['code','fname']].copy().drop_duplicates()
        codes = list(xref['code'].drop_duplicates())
        for code in codes:
            fname = xref.ix[xref['code'] == code, 'fname'].values[0]
            self.surgs[code] = Surrogate(code, fname)
            
class Surrogate:
    def __init__(self, srg_id, srg_path):
        '''
        Read a gridding surrogate into a dataframe
        '''
        print 'Loading surrogate: %s' %srg_path
        self.id = srg_id
        self.table = self._read_srg(srg_path)

    def _read_srg(self, srg_path):
        '''
        Parse in the gridding surrogate from the mixed formatted (is it tab? is it space?)
         surrogate files
        '''
        srg = []
        with open(srg_path) as f:
            for line in f:
                if line.startswith('#'):
                    if line[:6].startswith('#GRID') and 'GRIDDESC' not in line:
                        line = line.replace('\t', ' ')
                        line = [cell.strip() for cell in line.split(' ') if cell.strip()]
                        self.xorig = float(line[2])
                        self.yorig = float(line[3])
                        self.cell = float(line[4])
                        self.cols = int(line[6])
                        self.rows=int(line[7])
                else:
                    line = line.split('!')[0]
                    line = line.replace('\t', ' ')
                    line = [cell.strip() for cell in line.split(' ') if cell.strip()]
                    if len(line) != 5:
                        raise ValueError, 'Invalid surrogate format: %s' %' '.join(line)
                    srg.append(line)
        df = pd.DataFrame(srg, columns=['code','region_cd','col','row','frac'])
        df['frac'] = df['frac'].astype('f')
        df['region_cd'] = df['region_cd'].str.zfill(6)
        return df[['region_cd','col','row','frac']].copy()

def match_surrogate(df, xref):
    '''
    Match data to temporal based on xref keys using a hierarchy of keys
    * This is an incomplete hierarchy based on region, scc, and facility only.
        This hierarchy should be fine for HAPS sources in the current NEI.
    '''
    hierarchy = [['region_cd','scc'],['scc',],['region_cd',]]
    matched_df = pd.DataFrame()
    for match_cols in hierarchy:
        df = pd.merge(df, xref[match_cols+['code',]], on=match_cols, how='left')
        if matched_df.empty:
            matched_df = df[df['code'].notnull()].copy()
        else:
            matched_df = pd.concat((matched_df, df[df['code'].notnull()]))
        df = df[df['code'].isnull()].copy()
        df.drop('code', axis=1, inplace=True)
        if df.empty:
            break
    return matched_df

def grid_sources(emis, surg, grid_info):
    '''
    Apply the gridding matrix for the surrogate code group
    '''
    df = pd.DataFrame()
    emis = emis.copy().groupby(['region_cd','scc','code','run_group'], 
      as_index=False, sort=False).sum()
    code_list = list(emis['code'].drop_duplicates())
    for code in code_list:
        code_emis = emis[emis['code'] == code].copy()
        try:
            src_surg = surg.surgs[code]
        except KeyError:
            raise KeyError, 'Missing surrogate for code %s' %code
        if src_surg.cell != grid_info.xcell:
            raise ValueError, 'Gridding surrogate cell size does not match grid cell size'
        code_emis = pd.merge(code_emis, src_surg.table, on='region_cd', how='left')
        code_emis['ann_value'] = code_emis['ann_value'] * code_emis['frac']
        if not code_emis[(code_emis['col'].isnull()) | (code_emis['row'].isnull())].empty:
            bad_cell = code_emis[(code_emis['col'].isnull()) | (code_emis['row'].isnull())]
            print 'WARNING: Dropping %s records with null column or row numbers.' %len(bad_cell)
            print bad_cell[:20] 
            code_emis = code_emis[(code_emis['col'].notnull()) & (code_emis['row'].notnull())].copy()
        code_emis['col'] = code_emis['col'].astype('i') + calc_offset(src_surg.xorig, grid_info.xorig, grid_info.xcell)
        code_emis['row'] = code_emis['row'].astype('i') + calc_offset(src_surg.yorig, grid_info.yorig, grid_info.xcell)
        code_emis = code_emis[(code_emis['col'] > 0) & (code_emis['col'] <= grid_info.ncols) & \
          (code_emis['row'] > 0) & (code_emis['row'] <= grid_info.nrows)].copy()
        df = pd.concat((df, code_emis))
    return df

def calc_offset(surg_orig, grid_orig, cell_size):
    '''
    Offset number of grid cells from surrogate grid to emissions grid
    '''
    return int((surg_orig - grid_orig)/cell_size)









 
