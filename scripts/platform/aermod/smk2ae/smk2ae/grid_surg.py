from __future__ import division
from __future__ import print_function
from builtins import object
import os.path
import numpy as np
import pandas as pd

class GridSurg(object):
    '''
    Import the temporal crossreference file and merge against the temporal profiles
    For more information about these files see the SMOKE documentation.
    '''
    def __init__(self, gref, srg_path, srg_desc, scc_list):
        self.xref = self._get_gref(gref, scc_list)
        if self.xref.empty:
            raise ValueError('No surrogates found for inventory SCCs in xref')
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
                return '%0.5d' %fips
            else:
                return ''

    def _get_gref(self, gref, scc_list):
        '''
        Read in the gridding x-ref file (AGREF)
        '''
        names = ['region_cd','scc','code']
        df = pd.read_csv(gref, names=names, comment='#', sep=';', 
            dtype={'region_cd': str, 'scc': str, 'code': str})
        df['region_cd'] = df['region_cd'].apply(self._fix_fips) 
        df.loc[df['scc'].str.startswith('00'), 'scc'] = df.loc[df['scc'].str.startswith('00'), 
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
              dtype={'code': str})
        except:
            df = pd.read_csv(srg_desc, names=names, comment='#', sep=',', usecols=['code','fname'],
              dtype={'code': str})
        df['fname'] = '%s/' %srg_path + df['fname']
        return df

    def _load_surgs(self):
        '''
        Load the surrogates by code and file name into a surrogate object
         stored in a srg code indexed dictionary
        '''
        # Only load the codes that are needed based on the sources selected in the xref
        xref = pd.merge(self.xref, self.desc, on='code', how='left')
        xref = xref[['code','fname']].copy().drop_duplicates()
        codes = list(xref['code'].drop_duplicates())
        if len(xref.loc[xref['fname'].isnull()]) > 0:
            print(xref.loc[xref['fname'].isnull()])
            raise ValueError('Missing filenames from grid description')
        for code in codes:
            fname = xref.loc[xref['code'] == code, 'fname'].values[0]
            self.surgs[code] = Surrogate(code, fname)
            
class Surrogate(object):
    def __init__(self, srg_id, srg_path):
        '''
        Read a gridding surrogate into a dataframe
        '''
        print('Loading surrogate: %s' %srg_path)
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
                        raise ValueError('Invalid surrogate format: %s' %' '.join(line))
                    srg.append(line)
        df = pd.DataFrame(srg, columns=['code','region_cd','col','row','frac'])
        df['frac'] = df['frac'].astype('f')
        df['region_cd'] = df['region_cd'].str.zfill(5)
        return df[['region_cd','col','row','frac']].copy()

def match_surrogate(df, xref):
    '''
    Match data to temporal based on xref keys using a hierarchy of keys
    * This is an incomplete hierarchy based on region, scc, and facility only.
        This hierarchy should be fine for HAPS sources in the 2014v2 NEI.
    '''
    hierarchy = [['region_cd','scc'],['scc',],['region_cd',]]
    matched_df = pd.DataFrame()
    # Iterate over the passed hierarchies until a successful match is made
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
    # Notify that there was an unmatched source
    if not df.empty:
        print(df.head())
        raise ValueError('Unmatched sources for gridding')
    return matched_df

def grid_sources(emis, surg, grid_info):
    '''
    Apply the gridding matrix for the surrogate code group
    '''
    df = pd.DataFrame()
    emis = emis.copy().groupby(['region_cd','scc','code','run_group'], 
      as_index=False, sort=False).sum()
    code_list = list(emis['code'].drop_duplicates())
    # Iterate over the surrogate codes associated with the emissions
    for code in code_list:
        code_emis = emis[emis['code'] == code].copy()
        # Look for the code, hopefully it was read and subsets
        try:
            src_surg = surg.surgs[code]
        except KeyError:
            raise KeyError('Missing surrogate for code %s' %code)
        if src_surg.cell != grid_info.XCELL:
            raise ValueError('Gridding surrogate cell size does not match grid cell size')
        code_emis = pd.merge(code_emis, src_surg.table, on='region_cd', how='left')
        # Identify sources the cross reference successfully to a surrogate but are either missing a fraction
        #   or have a fraction that equals zero. These sources will ened up with dropped emissions
        zero_fracs = code_emis[['region_cd','scc','code','ann_value','frac']].groupby(['region_cd',
          'scc','code'], as_index=False).sum()
        zero_fracs['frac'] = zero_fracs['frac'].fillna(0)
        zero_fracs = zero_fracs[zero_fracs['frac'] == 0].copy()
        if len(zero_fracs) > 0:
            print('WARNING: Gridding the following sources will result in dropped emissions.')
            print(zero_fracs)
        # Multiply the emissions that use this surrogate code by the county->cell gridding fractions
        code_emis['ann_value'] = code_emis['ann_value'] * code_emis['frac']
        # Only keep emissions with a value over 0.
        code_emis = code_emis[code_emis['ann_value'] > 0].copy()
        if not code_emis[(code_emis['col'].isnull()) | (code_emis['row'].isnull())].empty:
            bad_cell = code_emis[(code_emis['col'].isnull()) | (code_emis['row'].isnull())]
            print('WARNING: Dropping %s records with null column or row numbers.' %len(bad_cell))
            print(bad_cell[:20]) 
            code_emis = code_emis[(code_emis['col'].notnull()) & (code_emis['row'].notnull())].copy()
        # Fill cells with emissions for this surrogate code
        code_emis['col'] = code_emis['col'].astype('i') + calc_offset(src_surg.xorig, grid_info.XORIG, grid_info.XCELL)
        code_emis['row'] = code_emis['row'].astype('i') + calc_offset(src_surg.yorig, grid_info.YORIG, grid_info.XCELL)
        code_emis = code_emis[(code_emis['col'] > 0) & (code_emis['col'] <= grid_info.NCOLS) & \
          (code_emis['row'] > 0) & (code_emis['row'] <= grid_info.NROWS)].copy()
        df = pd.concat((df, code_emis))
    return df

def calc_offset(surg_orig, grid_orig, cell_size):
    '''
    Offset number of grid cells from surrogate grid to emissions grid
    The emissions grid should nest and subset under the surrogate grid
    '''
    return int((surg_orig - grid_orig)/cell_size)









 
