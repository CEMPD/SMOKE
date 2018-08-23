#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
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
            raise KeyError('Missing environment variable: %s' %env_var)

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
    See SMOKE and IOAPI documentation for more information
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
                # Pad state FIPS to 2 characters
                if line.startswith('0'):
                    stfips = '%0.2d' %int(line[1:3])
                    st_dict[stfips] = line[3:5].strip()
    return pd.DataFrame(list(st_dict.items()), columns=['stfips','state'])

def insert_states(df, costcy_fname):
    '''
    Insert the state postal codes into the dataframe
    '''
    costcy = parse_costcy(costcy_fname)
    df['stfips'] = df['region_cd'].str[:2].str.zfill(2)
    return pd.merge(df, costcy, on='stfips', how='left')

def init_paths(work_path, sub_dirs):
    '''
    Create all of the subdirectories needed for this run
    '''
    for sub_dir in sub_dirs:
        sub_path = os.path.join(work_path, sub_dir)
        if not os.path.exists(sub_path):
            os.makedirs(sub_path) 

def proc_point_sector(inv_list, grid_name, grid_desc, state_fips):
    '''
    Process a point sector that may have mutiple run groups
    '''
    from smk2ae.point_io import get_temp_codes, proc_points, get_grid_location
    from smk2ae.pt_ff10 import AnnualFF10
    grid_info = Grid(grid_name, grid_desc)
    # Check for point sector specific variables
    var_list = ('PTPRO_HOURLY','PTPRO_WEEKLY','PTPRO_MONTHLY','PTREF','REPORT_PATH','CASE')
    check_env_vars(var_list)
    init_paths(os.environ['WORK_PATH'], ['xwalk','qa','locations','parameters','temporal'])
    invtable = get_invtable(os.environ['INVTABLE'])
    # Read in the point FF10
    inv = AnnualFF10(os.environ['SECTOR'], inv_list, invtable, state_fips)
    if inv.stk.empty:
        raise ValueError('No facilities found containing HAPS. Check state list and inventories.')
    inv.stk['facility_name'] = '"' + inv.stk['facility_name'] + '"'
    # Get the state names and filter down if specific states were specified
    inv.stk = insert_states(inv.stk, os.environ['COSTCY'])
    if state_fips:
        inv.state = inv.stk.loc[inv.stk['stfips'] == inv.st_fips, 'state'].drop_duplicates().values[0]
        print('Running for state: %s' %inv.state)
    from smk2ae.hourly_rep import HourlyReports 
    hrly = HourlyReports(os.environ['SECTOR'], os.environ['CASE'], os.environ['BASE_YEAR'], 
        os.environ['REPORT_PATH'], state_fips)
    from smk2ae.temporal import Temporal
    temp = Temporal(os.environ['PTREF'],os.environ['PTPRO_HOURLY'],os.environ['PTPRO_WEEKLY'],
        os.environ['PTPRO_MONTHLY'], inv.scc_list, inv.fips_list, inv.fac_list)
    inv.stk = get_temp_codes(inv.stk, temp.xref.copy())
    if not hrly.egu_units.empty:
        inv.stk = pd.merge(inv.stk, hrly.egu_units, on='unit_id', how='left')
        # This is in here to keep EGU units unique temporally
        inv.stk.loc[inv.stk['uniq'] == 'Y', 'ALLDAY'] = -1
        inv.stk.loc[inv.stk['uniq'] == 'Y', 'WEEKLY'] =  inv.stk.loc[inv.stk['uniq'] == 'Y', 
          'unit_id'].astype('i')
        inv.stk.loc[inv.stk['uniq'] == 'Y', 'MONTHLY'] =  inv.stk.loc[inv.stk['uniq'] == 'Y', 
          'rel_point_id'].astype('i')
    # Process airports if the airport locations file was defined in the driver script
    try:
        os.environ['AIRPORT_LOCS']
    except KeyError:
        print('WARNING: No Airport Locations File defined')
    else:
        from smk2ae.airports import Airports
        air = Airports(os.environ['AIRPORT_LOCS'])
        air.match_airport_inv(inv)
        inv.drop_facilities(air.fac_list)
        if air.fac_list:
            air.calc_cell_xref(grid_info)
            air.write_runways(grid_info, temp.profs)
            air.write_no_runways(grid_info, temp.profs)
    # Process all other point sources
    if inv.fac_list:
        inv.stk = get_grid_location(inv.stk, grid_info)
        proc_points(inv, temp, hrly, os.environ['WORK_PATH'], grid_info)

def get_invtable(fn):
    '''
    Read in the inventory table for the kept pollutants
    See SMOKE documentation for more information about the INVTABLE
    '''
    invtable = pd.read_fwf(fn, comment='#', colspecs=[(0,11), (16,32), (41,42), (43,49)], 
      names=['smoke_name','poll','keep','spec_factor'], 
      converters={'smoke_name': str.strip, 'poll': str.strip, 'keep': str.strip})
    invtable = invtable[invtable['keep'].str.upper() == 'Y'].copy()
    invtable.drop_duplicates(inplace=True)
    return invtable

def release_params(fn, rg_suffix):
    '''
    Get the release parameters file used to define release height and sigma z
     for area sources
    '''
    df = pd.read_csv(fn)
    df.rename(columns={'release_height': 'rel_ht', 'sigma_z': 'sz'}, inplace=True)
    df['run_group'] = df['run_group'] + rg_suffix
    return df

def get_sw_corner(df, met_grid):
    '''
    Find the SW corner of the grid cell associated with the lat-lon. This is used for 
      gridding and determining the UTM zone
    '''
    df['met_swcorner_lon'] = met_grid.colrow_to_ll(df['met_col'].astype('f'), 
      df['met_row'].astype('f'))['lon'].astype('f')
    df['met_cell'] = 'G' + df['met_col'].astype(str).str.zfill(3) + 'R' + df['met_row'].astype(str).str.zfill(3)
    return df

def proc_area_sector(inv_list, grid_name, grid_desc, state_fips):
    '''
    Process an area sector or run group
    '''
    from smk2ae.ar_ff10 import AnnualFF10
    from smk2ae.grid_surg import GridSurg, match_surrogate, grid_sources
    from smk2ae.source_groups import SourceGroups
    # Check for area specific environment variables
    var_list = ('ATPRO_HOURLY','ATPRO_WEEKLY','ATPRO_MONTHLY','ATREF','SRGPRO','AGREF','SRGDESC',
      'SOURCE_GROUPS','COSTCY','RUN_GROUPS','GROUP_PARAMS','RUN_GROUP_SUFFIX')
    check_env_vars(var_list)
    # Setup the output paths
    init_paths(os.environ['WORK_PATH'], ['emis','qa','locations','parameters','temporal','xwalk'])
    invtable = get_invtable(os.environ['INVTABLE'])
    # Read the nonpoint FF10s
    inv = AnnualFF10(inv_list, invtable, state_fips)
    scc_list = list(inv.emis['scc'].drop_duplicates()) 
    # Load all of the appropriate gridding surrogates
    surg = GridSurg(os.environ['AGREF'], os.environ['SRGPRO'], os.environ['SRGDESC'], scc_list)
    src_groups = SourceGroups(os.environ['SOURCE_GROUPS'])
    src_groups.xref['run_group'] = src_groups.xref['run_group'] + os.environ['RUN_GROUP_SUFFIX']
    inv.merge_rungroups(src_groups.xref)
    # Narrow down to just the selected run groups.
    if os.environ['RUN_GROUPS']:
        run_groups = [rg.strip() for rg in os.environ['RUN_GROUPS'].split(',')]
    else:
        run_groups = list(inv.emis['run_group'].drop_duplicates())
    # Setup the met grid. This is the IOAPI grid specified, but may not be the same as the cell grid
    #  if this is a sector with 4 km cells such as HDON, OILGAS, etc.
    met_grid = Grid(grid_name, grid_desc)
    for run_group in run_groups:
        print('NOTE: Running for %s' %run_group)
        # Gridding
        run_emis = match_surrogate(inv.emis[inv.emis['run_group'] == run_group], surg.xref)
        if grid_name.startswith('12'):
            # Grid at 4 km for HDON LDON HDOFF LDOFF HOTEL and OILGAS
            if run_group[:3] in ('HDO','LDO','HOT','OIL'):
                cell_grid = Grid('4%s' %grid_name[2:], grid_desc)
            else:
                cell_grid = met_grid
        else:
            cell_grid = met_grid
        run_emis = grid_sources(run_emis, surg, cell_grid)
        # If the cell and met grid are the same, then set up the parameters to be equivalent
        if cell_grid.GDNAM == met_grid.GDNAM:
            run_emis['met_col'] = run_emis['col']
            run_emis['met_row'] = run_emis['row']
            run_emis['src_id'] = '%s_1' %int(met_grid.XCELL/1000.)
        # Otherwise find the met cell based on the cell grid
        elif (cell_grid.XCELL == 4000.) and (met_grid.XORIG == cell_grid.XORIG) and (met_grid.XCELL == 12000.):
            run_emis = define_met_cell(run_emis, cell_grid, met_grid)
        else:
            raise ValueError('Unsupported grid cell size')
        swcorners = get_sw_corner(run_emis[['met_col','met_row']].copy().drop_duplicates(), met_grid)
        run_emis = pd.merge(run_emis, swcorners, on=['met_col','met_row'], how='left')
        grid_keys = ['region_cd','scc','src_id','met_cell']
        write_aermod_emis(inv, cell_grid.XCELL, run_emis[grid_keys+['frac',]].drop_duplicates(grid_keys), src_groups.xref)
        # Onroad run groups have temporal helper files calculated outside of this processor
        if run_group[:3] in ('HDO','LDO','HOT'):
            temp = None 
        else:
            from smk2ae.temporal import Temporal
            # Temporalization
            scc_list = list(run_emis['scc'].drop_duplicates()) 
            fips_list = list(run_emis['region_cd'].drop_duplicates()) 
            temp = Temporal(os.environ['ATREF'],os.environ['ATPRO_HOURLY'],os.environ['ATPRO_WEEKLY'],
                os.environ['ATPRO_MONTHLY'], scc_list, fips_list)
        from smk2ae.area_io import proc_area
        emis_keys = ['region_cd','scc','run_group','col','row','met_cell','src_id']
        run_emis = run_emis[emis_keys+['ann_value','met_swcorner_lon']].copy().drop_duplicates(emis_keys)
        params = release_params(os.environ['GROUP_PARAMS'], os.environ['RUN_GROUP_SUFFIX'])
        proc_area(run_emis, cell_grid, temp, params)

def define_met_cell(df, grid_info, met_grid):
    '''
    Define met col and row plus any subcell numbering
    '''
    df['met_col'] = calc_met_cell(df['col'], grid_info.XCELL, met_grid.XORIG, grid_info.XORIG, met_grid.XCELL).astype('i')
    df['met_row'] = calc_met_cell(df['row'], grid_info.XCELL, met_grid.YORIG, grid_info.YORIG, met_grid.XCELL).astype('i')
    df['subcell'] = calc_sub_cell(df['col'], df['row'], df['met_col'], df['met_row'], grid_info.XCELL, met_grid.XCELL)
    grid_size = int(grid_info.XCELL/1000.)
    df['src_id'] = '%s_' %grid_size + df['subcell'].astype(str)
    df.drop('subcell', axis=1, inplace=True)
    return df

def calc_met_cell(cell_dim, cell_size, met_orig, grid_orig, met_size):
    '''
    Calculate the met cell for a sub cell
    This is essentially a grid->grid conversion. Only tested on the same projection.
    '''
    return np.floor((((cell_dim.astype('i') - 1) * cell_size) + (met_orig - grid_orig))/met_size) + 1

def calc_sub_cell(col, row, met_col, met_row, grid_xcell, met_xcell):
    '''
    Number the subcells if needed
    '''
    mult = int(met_xcell/grid_xcell)
    col = (col - 1) - ((met_col - 1) * mult) 
    row = (row - 1) - ((met_row - 1) * mult)
    return (row * mult) + col + 1

def write_aermod_emis(inv, cell_size, run_emis, xref):
    '''
    Write the "Post AERMOD" emissions files
    '''
    key_cols = ['run_group','region_cd','met_cell','src_id','source_group','smoke_name']
    out_cols = key_cols + ['ann_value',]
    group_sccs = list(run_emis['scc'].drop_duplicates())
    # Merge in the pollutants
    group_df = pd.merge(run_emis, inv.emis, on=['region_cd','scc'], how='left')
    group_df['ann_value'] = group_df['ann_value'] * group_df['frac']
    # Fill in the monthly emissions if this is a sector with montly emissions
    if inv.monthly:
        mons = ['january','february','march','april','may','june','july','august','september',
          'october','november','december']
        out_cols += mons
        for mon in mons:
            group_df[mon] = group_df[mon] * group_df['frac']
    group_df = pd.merge(group_df, xref[['scc','source_group']], on='scc', how='left')
    group_df.drop(['frac','scc'], axis=1, inplace=True)
    # Aggregate up emissions by source for writing
    group_df = group_df.groupby(['run_group','region_cd','source_group','met_cell','src_id','smoke_name'],
      as_index=False, sort=False).sum()
    run_group = group_df['run_group'].values[0]
    fname = os.path.join(os.environ['WORK_PATH'],'emis','%s_%s_emis.csv' %(int(cell_size/1000.),
      run_group))
    group_df.to_csv(fname, columns=out_cols, index=False, float_format='%.12g')

def main():
    var_list = ('GRIDDESC','REGION_IOAPI_GRIDNAME','EMISINV_A','INVTABLE','BASE_YEAR',
      'SECTOR','WORK_PATH','STATE_FIPS')
    check_env_vars(var_list)
    state_fips = os.environ['STATE_FIPS']
    inv_list = get_inv_list()
    grid_name = os.environ['REGION_IOAPI_GRIDNAME']
    grid_desc = os.environ['GRIDDESC']
    # Select point processing if the sector name falls under the list below
    if os.environ['SECTOR'].lower() in ('ptnonipm','pt_oilgas','ptegu','ptegu_pk','ptipm','point'):
        proc_point_sector(inv_list, grid_name, grid_desc, state_fips)
    # Otherwise process as area/nonpoint
    else:
        proc_area_sector(inv_list, grid_name, grid_desc, state_fips)

if __name__ == '__main__':
	main()


