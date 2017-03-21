#!/usr/bin/env python

import os, sys, string
import os.path
import numpy as np
import pandas as pd
import smk2ae

from smk2ae.grid import Grid

########


def check_env_vars(vars):
    '''
    Verify that all required environment variables are set
    '''
    for env_var in vars:
        try:
            os.environ[env_var]
        except KeyError:
            raise KeyError, 'Missing environment variable: %s' %env_var

def get_inv_list():
    '''
    Get a list of all the inventory files from the environment variables
    '''
    inv_list = [] 
    for inv_var in ['EMISINV_%s' %alpha for alpha in string.ascii_uppercase]:
        if inv_var in os.environ:
            inv_list.append(os.environ[inv_var])
    return inv_list

def parse_costcy(costcy):
    '''
    Parse in the COSTCY file and get a postal code to state fips x-ref
    '''
    with open(costcy) as infile:
        st_dict = {}
        cty_dict = {}
        read_state = False
        read_county = False
        for line in infile:
            if line.startswith('/STATE/'):
                read_state = True
            elif line.startswith('/COUNTY/'):
                read_state = False
                read_county = True
            elif read_state:
                if line.startswith('0'):
                    stfips = '%0.2d' %int(line[1:3])
                    st_dict[stfips] = line[3:5].strip()
    return pd.DataFrame(st_dict.items(), columns=['stfips','state'])

def insert_states(df, costcy_fname):
    '''
    Insert the state postal codes into the dataframe
    '''
    costcy = parse_costcy(costcy_fname)
    df['stfips'] = df['region_cd'].str[:2].str.zfill(2)
    return pd.merge(df, costcy, on='stfips', how='left')

def init_paths(work_path, sub_dirs):
    for sub_dir in sub_dirs:
        sub_path = os.path.join(work_path, sub_dir)
        if not os.path.exists(sub_path):
            os.makedirs(sub_path) 

def proc_point_sector(inv_list, grid_name, grid_desc, state_fips):
    '''
    Process a point sector
    '''
    from smk2ae.point_io import get_temp_codes, proc_points, get_grid_location
    from smk2ae.pt_ff10 import AnnualFF10
    grid_info = Grid(grid_name, grid_desc)
    var_list = ('PTPRO_HOURLY','PTPRO_WEEKLY','PTPRO_MONTHLY','PTREF','REPORT_PATH','CASE')
    check_env_vars(var_list)
    init_paths(os.environ['WORK_PATH'], ['xwalk','qa','locations','parameters','temporal'])
    inv = AnnualFF10(os.environ['SECTOR'], inv_list, os.environ['HAPS_LIST'], state_fips)
    if inv.stk.empty:
        raise ValueError, 'No facilities found containing HAPS. Check state list and inventories.'
    inv.stk['facility_name'] = '"' + inv.stk['facility_name'] + '"'
    inv.stk = insert_states(inv.stk, os.environ['COSTCY'])
    if state_fips:
        inv.state = inv.stk.ix[inv.stk['stfips'] == inv.st_fips, 'state'].drop_duplicates().values[0]
        print 'Running for state: %s' %inv.state
    from smk2ae.hourly_rep import HourlyReports 
    hrly = HourlyReports(os.environ['SECTOR'], os.environ['CASE'], os.environ['BASE_YEAR'], 
        os.environ['REPORT_PATH'], state_fips)
    from smk2ae.temporal import Temporal
    temp = Temporal(os.environ['PTREF'],os.environ['PTPRO_HOURLY'],os.environ['PTPRO_WEEKLY'],
        os.environ['PTPRO_MONTHLY'], inv.scc_list, inv.fips_list, inv.fac_list)
    inv.stk = get_temp_codes(inv.stk, temp.xref.copy())
    if not hrly.egu_units.empty:
        print inv.stk[:5]
        print hrly.egu_units[:5]
        inv.stk = pd.merge(inv.stk, hrly.egu_units, on='unit_id', how='left')
        # This is in here to keep EGU units unique temporally
        inv.stk.ix[inv.stk['uniq'] == 'Y', 'ALLDAY'] = -1
        inv.stk.ix[inv.stk['uniq'] == 'Y', 'WEEKLY'] =  inv.stk.ix[inv.stk['uniq'] == 'Y', 
          'unit_id'].astype('i')
        inv.stk.ix[inv.stk['uniq'] == 'Y', 'MONTHLY'] =  inv.stk.ix[inv.stk['uniq'] == 'Y', 
          'rel_point_id'].astype('i')
    # Process airports
    try:
        os.environ['AIRPORT_LOCS']
    except KeyError:
        print 'WARNING: No Airport Locations File defined'
    else:
        from smk2ae.airports import Airports
        air = Airports(os.environ['AIRPORT_LOCS'])
        air.match_airport_inv(inv)
        inv.drop_facilities(air.fac_list)
        if air.fac_list:
            air.calc_cell_xref(grid_info)
            air.write_runways(grid_info, temp.aermod)
            air.write_no_runways(grid_info, temp.aermod)
    # Process all other point sources
    if inv.fac_list:
        inv.stk = get_grid_location(inv.stk, grid_info)
        proc_points(inv, temp, hrly, os.environ['WORK_PATH'], grid_info)

def proc_area_sector(inv_list, grid_name, grid_desc, state_fips):
    '''
    Process an area sector
    '''
    from smk2ae.ar_ff10 import AnnualFF10
    from smk2ae.grid_surg import GridSurg, match_surrogate, grid_sources
    from smk2ae.hem_groups import HEMGroups
    var_list = ('ATPRO_HOURLY','ATPRO_WEEKLY','ATPRO_MONTHLY','ATREF','SRGPRO','AGREF','SRGDESC',
      'HEM_GROUPS','COSTCY','RUN_GROUPS','CELL_GRID')
    check_env_vars(var_list)
    init_paths(os.environ['WORK_PATH'], ['emis','qa','locations','parameters','temporal','xwalk'])
    inv = AnnualFF10(inv_list, os.environ['HAPS_LIST'], state_fips)
    surg = GridSurg(os.environ['AGREF'], os.environ['SRGPRO'], os.environ['SRGDESC'], inv.scc_list)
    hem = HEMGroups(os.environ['HEM_GROUPS'])
    inv.emis = pd.merge(inv.emis, hem.xref, on='scc', how='left')
    if not inv.emis[(inv.emis['run_group'].isnull()) | (inv.emis['run_group'] == '')].empty:
        print 'WARNING: Unmatched SCCs to HEM groups'
        print inv.emis.ix[inv.emis['run_group'].isnull(),'scc'].drop_duplicates()
    # Narrow down to just the selected run groups.
    if os.environ['RUN_GROUPS']:
        run_groups = [rg.strip() for rg in os.environ['RUN_GROUPS'].split(',')]
        inv.emis = inv.emis[inv.emis['run_group'].isin(run_groups)].copy()
        inv.update_lists()
    else:
        run_groups = list(inv.emis['run_group'].drop_duplicates())
    inv.emis = match_surrogate(inv.emis, surg.xref)
    met_grid = Grid(grid_name, grid_desc)
    cell_grid = Grid(os.environ['CELL_GRID'], grid_desc)
    grid_emis = grid_sources(inv.emis, surg, cell_grid)
    if (met_grid.gdnam == cell_grid.gdnam) and (met_grid.xcell == 12000.):
        grid_emis['met_col'] = grid_emis['col']
        grid_emis['met_row'] = grid_emis['row']
        grid_emis['src_id'] = '12_1'
    elif (cell_grid.xcell == 4000.) and (met_grid.xorig == cell_grid.xorig) and (met_grid.xcell == 12000.):
        grid_emis = define_met_cell(grid_emis, cell_grid, met_grid)
    else:
        raise ValueError, 'Unsupported grid cell size'
    centroids = grid_emis[['met_col','met_row']].copy().drop_duplicates()
    centroids['met_centroid_lon'] = met_grid.colrow_to_ll(centroids['met_col'].astype('f') + 0.5, 
      centroids['met_row'].astype('f') + 0.5)['lon'].astype('f')
    centroids['met_cell'] = 'G' + centroids['met_col'].astype(str).str.zfill(3) + 'R' + centroids['met_row'].astype(str).str.zfill(3)
    grid_emis = pd.merge(grid_emis, centroids, on=['met_col','met_row'], how='left')
    grid_keys = ['region_cd','scc','src_id','met_cell']
    write_aermod_emis(inv, cell_grid.xcell, run_groups, 
      grid_emis[grid_keys + ['frac',]].drop_duplicates(grid_keys), hem.xref)
    # Temporalization
    from smk2ae.temporal import Temporal
    temp = Temporal(os.environ['ATREF'],os.environ['ATPRO_HOURLY'],os.environ['ATPRO_WEEKLY'],
        os.environ['ATPRO_MONTHLY'], inv.scc_list, inv.fips_list)
    from smk2ae.area_io import proc_area
    emis_keys = ['region_cd','scc','run_group','col','row','met_cell','src_id']
    grid_emis = grid_emis[emis_keys+['ann_value','met_centroid_lon']].copy().drop_duplicates(emis_keys)
    proc_area(grid_emis, cell_grid, temp)

def define_met_cell(df, grid_info, met_grid):
    '''
    Define met col and row plus any subcell numbering
    '''
    df['met_col'] = calc_met_cell(df['col'], grid_info.xcell, met_grid.xorig, grid_info.xorig, met_grid.xcell).astype('i')
    df['met_row'] = calc_met_cell(df['row'], grid_info.xcell, met_grid.yorig, grid_info.yorig, met_grid.xcell).astype('i')
    df['subcell'] = calc_sub_cell(df['col'], df['row'], df['met_col'], df['met_row'], grid_info.xcell, met_grid.xcell)
    grid_size = int(grid_info.xcell / 1000.)
    df['src_id'] = '%s_' %grid_size + df['subcell'].astype(str)
    df.drop('subcell', axis=1, inplace=True)
    return df.sort(['met_col','met_row','src_id'])

def calc_met_cell(cell_dim, cell_size, met_orig, grid_orig, met_size):
    '''
    Calculate the met cell for a sub cell
    This is essentially a grid->grid conversion. Only tested on the same projection.
    '''
    return np.floor((((cell_dim.astype('i') - 1) * cell_size) + (met_orig - grid_orig)) / met_size) + 1

def calc_sub_cell(col, row, met_col, met_row, grid_xcell, met_xcell):
    '''
    Number the subcells if needed
    '''
    mult = int((met_xcell/grid_xcell))
    col = (col - 1) - ((met_col - 1) * mult) 
    row = (row - 1) - ((met_row - 1) * mult)
    return (row * mult) + col + 1

def write_aermod_emis(inv, cell_size, run_groups, grid_emis, hem):
    '''
    Write the "Post AERMOD" emissions files
    '''
    key_cols = ['run_group','region_cd','met_cell','src_id','source_group','smoke_name']
    src_groups = hem[['scc','source_group']].drop_duplicates()
    for run_group in run_groups:
        out_cols = key_cols + ['ann_value',]
        group_sccs = list(hem.ix[hem['run_group'] == run_group, 'scc'].drop_duplicates())
        group_df = pd.merge(inv.inv_emis[inv.inv_emis['scc'].isin(group_sccs)], 
          grid_emis[grid_emis['scc'].isin(group_sccs)], on=['region_cd','scc'], how='left')
        group_df['ann_value'] = group_df['ann_value'] * group_df['frac']
        if inv.seasons:
            seasons = ['winter','spring','summer','fall']
            out_cols += seasons
            for season in seasons:
                group_df[season] = group_df[season] * group_df['frac']
        group_df.drop('frac', axis=1, inplace=True)
        group_df = pd.merge(group_df, src_groups, on='scc', how='left') 
        group_df.drop('scc', axis=1, inplace=True)
        group_df = group_df.groupby(['region_cd','source_group','met_cell','src_id','smoke_name'],
          as_index=False, sort=False).sum()
        group_df['run_group'] = run_group
#        group_df.sort(key_cols, inplace=True)
        fname = os.path.join(os.environ['WORK_PATH'],'emis','%s_%s_emis.csv' %(int(cell_size/1000.),
          run_group))
        group_df.to_csv(fname, columns=out_cols, index=False)

def main():
    var_list = ('GRIDDESC','REGION_IOAPI_GRIDNAME','EMISINV_A','HAPS_LIST','BASE_YEAR',
      'SECTOR','WORK_PATH','STATE_FIPS')
    check_env_vars(var_list)
    state_fips = os.environ['STATE_FIPS']
    inv_list = get_inv_list()
    grid_name = os.environ['REGION_IOAPI_GRIDNAME']
    grid_desc = os.environ['GRIDDESC']
    if os.environ['SECTOR'].lower() in ('ptnonipm','pt_oilgas','ptegu','ptegu_pk','ptipm','point'):
        proc_point_sector(inv_list, grid_name, grid_desc, state_fips)
    else:
        proc_area_sector(inv_list, grid_name, grid_desc, state_fips)

if __name__ == '__main__':
	main()


