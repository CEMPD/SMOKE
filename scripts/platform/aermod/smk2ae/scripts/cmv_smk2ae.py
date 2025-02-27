#!/usr/bin/env python

import os, sys, string
from collections import Counter
from osgeo import ogr
import pandas as pd
import smk2ae
from smk2ae.grid import Grid
from smk2ae.utm import UTM
from smk2ae.ar_ff10 import AnnualFF10
from smk2ae.source_groups import SourceGroups
from smk2ae.temporal import Temporal, match_temporal

def label_emis(inv, shape_facs):
    '''
    Label the FF10 emissions with facid and src_id based on class and activity
    '''
    inv.loc[inv['scc'].str[6] == '2', 'class'] = 'c1'
    inv.loc[inv['scc'].str[6] == '3', 'class'] = 'c3'
    inv.loc[inv['scc'].str.endswith('100'), 'shape_id'] = inv.loc[inv['scc'].str.endswith('100'), 
      'shape_id'].astype(str).str.zfill(5)
    inv.loc[inv['scc'].str.endswith('200'), 'shape_id'] = inv.loc[inv['scc'].str.endswith('200'), 
      'shape_id'].astype(str).str.zfill(4)
    emis = inv.groupby(['region_cd','shape_id','scc','class','smoke_name'], as_index=False).sum()
    emis = pd.merge(emis, shape_facs, on='shape_id', how='left')
    emis['ann_value'] = emis['ann_value'] * emis['area_frac']
    ports = emis[emis['facid'].str.startswith('P')].copy()
    ports['src_id'] = ports['facid'].str[:6] + ports['class']
    uw = emis[emis['facid'].str.startswith('U')].copy()
    uw['src_id'] = uw['facid'].str[:5] + uw['facid'].str[-4:] + uw['class']
    emis = pd.concat((ports, uw))
    return emis[['region_cd','facid','src_id','scc','smoke_name','ann_value']]

def read_shapes(shapes_file):
    '''
    Read in the shapes from the shapes file

    The shapes file is formatted like this:
    facid,sourceid,loctype,lon,lat,zone,numvert,area
    U0612F10001P001,TRACT,L,-75.103862,38.995806,,45,586792801
    U0612F10001P001,TRACT,L,-75.311905,38.945821,,45,586792801

    Port facility IDs follow P[5-character-shapeid]F[5-character-fips]
    Underway facility IDs follow U[4-character-shapeid]F[5-character-fips]P[3-character-polygonid]
    '''
    df = pd.read_csv(shapes_file, usecols=['facid','lon','lat','numvert','area'], dtype={'facid': str})
    # Calc polygon area fraction for splitting a shape into polygons
    df['shape_id'] = df['facid'].str.split('F').str[0].str[1:].astype(int).astype(str)
    df.loc[df['facid'].str.startswith('P'), 'shape_id'] = df.loc[df['facid'].str.startswith('P'), 
      'shape_id'].str.zfill(5)
    df.loc[df['facid'].str.startswith('U'), 'shape_id'] = df.loc[df['facid'].str.startswith('U'), 
      'shape_id'].str.zfill(4)
    area = df[['facid','shape_id','area']].drop_duplicates(['facid','shape_id'])
    area_sum = area[['shape_id','area']].groupby('shape_id', as_index=False).sum()
    area = pd.merge(area, area_sum, on='shape_id', how='left', suffixes=['','_sum'])
    area['area_frac'] = area['area'] / area['area_sum']
    df = pd.merge(df, area[['facid','area_frac']], on='facid', how='left')
    return df

def get_met_cell(shapes_df, met_grid):
    '''
    Get the dominant met/CMAQ grid cell for each polygon ID.

    This function largely does the looping over the facility shapes and determines if
      the facility falls into more than one grid cell.
    This function calls get_dominant_cell if the facility shape falls in more than one
      grid cell.
    This function returns a new "shapes" dataframe where the dominant grid cell is selected
      for the entire facility.
    '''
    shapes_df[['x','y']] = shapes_df[['lon','lat']].apply(lambda x: \
      pd.Series(met_grid.get_coords(x['lon'],x['lat'])), axis=1)
    shapes_df[['col','row']] = shapes_df[['x','y']].apply(lambda x: \
      pd.Series(met_grid.get_cell(x['x'],x['y'])), axis=1)
    shapes_df['cell'] = shapes_df['col'].astype(str) + ',' + shapes_df['row'].astype(str)
    out_df = pd.DataFrame()
    for facid in list(shapes_df['facid'].drop_duplicates()):
        shape_df = shapes_df[shapes_df['facid'] == facid].copy()
        for src_id in list(shape_df['src_id'].drop_duplicates()):
            poly_df = shape_df[shape_df['src_id'] == src_id].copy()
            if len(list(poly_df['cell'].drop_duplicates())) > 1:
                dom_cell = get_dominant_cell(poly_df[['src_id','cell','x','y','col','row']], 
                  met_grid, facid)
                poly_df['col'] = dom_cell.split(',')[0]
                poly_df['row'] = dom_cell.split(',')[1]
            poly_df.drop('cell', inplace=True, axis=1)
            out_df = pd.concat((out_df, poly_df[['facid','src_id','col','row']].drop_duplicates()))
    return out_df

def get_dominant_cell(poly_df, met_grid, facid):
    '''
    Calculate the dominant cell based on max shape area in all possible cells
    
    The polygons within the facility and src_id shape are drawn on the CMAQ/met Lambert grid.
    All grid cells that contain the Lambert coordinates of the vertices are selected and drawn.
    The grid cell for the facility is selected based on the area of greatest intersection
      with the facility shape (all polygons).
    This has a couple of loops to identify features. There may be a better way to do this.
    '''
    max_area = 0
    src_id = poly_df['src_id'].drop_duplicates().values[0]
    cell_polys = gen_cell_polys(poly_df[['cell','col','row']].copy().drop_duplicates(), met_grid)
    # Iterate over the points to identify individual rings
    init_vert = False
    ring_list = []
    for idx, vert in poly_df.iterrows():
        cur_vert = list(vert[['x','y']])
        if init_vert:
            # An individual ring is put in the ring list when the current vertex equals the
            #  initial vertex
            if cur_vert == init_vert: 
                ring.AddPoint(*vert[['x','y']])
                ring_list.append(ring)
                init_vert = False
            else:
                ring.AddPoint(*vert[['x','y']])
        else:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(*vert[['x','y']])
            init_vert = cur_vert
    # Iterate over the rings to identify inner vs. outer rings
    poly_list = []
    for x, ring in enumerate(ring_list):
        if x == 0:
            uniq_poly = ogr.Geometry(ogr.wkbPolygon)
            uniq_poly.AddGeometry(ring)
            last_outer = 0
        elif x > 0:
            # Need to create temporary polygons to use the geometry "Within" function
            cur_ring_poly = ogr.Geometry(ogr.wkbPolygon)
            cur_ring_poly.AddGeometry(ring)
            last_outer_poly = ogr.Geometry(ogr.wkbPolygon)
            last_outer_poly.AddGeometry(ring_list[last_outer])
            # A ring is inner to the outer if it is within the outer. This will be part of 
            #  a single polygon whereas a new outer gets a new polygon.
            if cur_ring_poly.Within(last_outer_poly):
                uniq_poly.AddGeometry(ring)
            else:
                poly_list.append(uniq_poly)
                uniq_poly = ogr.Geometry(ogr.wkbPolygon)
                uniq_poly.AddGeometry(ring)
                last_outer = x
    poly_list.append(uniq_poly)
    if len(poly_list) > 1:
        src_poly = ogr.Geometry(ogr.wkbMultiPolygon)
        for poly in poly_list:
            src_poly.AddGeometry(poly)
    else:
        src_poly = poly_list[0]
    if src_poly.IsValid():
        for cell, cell_poly in cell_polys.items():
            overlap_area = src_poly.Intersection(cell_poly).GetArea()
            if  overlap_area > max_area:
                max_area = overlap_area
                dom_cell = cell
    else:
        print('WARNING: Shape %s %s contains an invalid polygon' %(facid, src_id))
        # For an invalid polygon use the grid cell that contains the most source points 
        cell_count = Counter(list(poly_df['cell']))
        dom_cell = cell_count.most_common()[0][0]
    return dom_cell

def gen_cell_polys(cell_df, met_grid):
    '''
    This function calculates the Lambert grid cell polygons as needed based
      on the CMAQ/met grid column and row
    This function returns a dictionary of grid cell polygons by cell ID 
    '''
    cell_polys = {}
    for idx, cell in cell_df.iterrows():
        cell_poly = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(*met_grid.colrow_to_coords(cell['col'],cell['row']))
        ring.AddPoint(*met_grid.colrow_to_coords(cell['col'],cell['row']+1))
        ring.AddPoint(*met_grid.colrow_to_coords(cell['col']+1,cell['row']+1))
        ring.AddPoint(*met_grid.colrow_to_coords(cell['col']+1,cell['row']))
        ring.AddPoint(*met_grid.colrow_to_coords(cell['col'],cell['row']))
        cell_poly.AddGeometry(ring)
        cell_polys[cell['cell']] = cell_poly
    return cell_polys

def write_locations(df):
    '''
    Write the locations file
    '''
    df.drop_duplicates(['region_cd','facid','src_id','type','col','row','utm_zone'], inplace=True)
    out_cols = ['state','region_cd','facid','src_id','type','col','row','utmx','utmy','utm_zone',
      'lon','lat']
    fname = os.path.join(os.environ['WORK_PATH'], 'locations', 'CMV_locations.csv')
    df.to_csv(fname, columns=out_cols, index=False)

def write_parameters(df):
    '''
    Write the parameters file
    
    The maximum number of vertices for a polygon needs to be calculated first
    The coordinates for each polygon is written to a single line
    '''
    # Set the release height and sigma z based on source group names
    df.loc[df['src_id'].str.endswith('c1'), 'sz'] = 3.907
    df.loc[df['src_id'].str.endswith('c1'), 'rel_ht'] = 8.4
    df.loc[df['src_id'].str.endswith('c3'), 'sz'] = 40.7
    df.loc[df['src_id'].str.endswith('c3'), 'rel_ht'] = 20
    fname = os.path.join(os.environ['WORK_PATH'], 'parameters', 'CMV_area_params.csv')
    maxverts = df['numvert'].max()
    coord_cols = []
    [coord_cols.extend(['utmx','utmy']) for vert_num in range(1, maxverts + 1)]
    ll_cols = []
    [ll_cols.extend(['lon','lat']) for vert_num in range(1, maxverts + 1)]
    out_cols = ['state','region_cd','facid','src_id','type','area','rel_ht','numvert','sz'] + coord_cols\
      + ll_cols
    with open(fname,'w') as f:
        f.write('%s\n' %','.join(out_cols))
        for facid in list(df['facid'].drop_duplicates()):
            # Iterate over vertices in a polygon. This isn't consistent across polygons
            for src_id in list(df.ix[df['facid'] == facid, 'src_id'].drop_duplicates()):
                src_df = df[(df['facid']==facid) & (df['src_id']==src_id)].copy()
                numverts = src_df['numvert'].values[0]
                out_line = [src_df['state'].values[0], src_df['region_cd'].values[0], facid, 
                  src_id, src_df['type'].values[0], str(src_df['area'].values[0]),
                  str(src_df['rel_ht'].values[0]), str(numverts),str(src_df['sz'].values[0])]
                for ix, row in src_df.iterrows():
                    out_line.extend([str(row['utmx']),str(row['utmy'])])
                if numverts < maxverts:
                    for x in range(maxverts - numverts):
                        out_line.extend(['',''])
                for ix, row in src_df.iterrows():
                    out_line.extend([str(row['lon']),str(row['lat'])])
                if numverts < maxverts:
                    for x in range(maxverts - numverts):
                        out_line.extend(['',''])
                f.write('%s\n' %','.join(out_line))
        
def calc_monthly_temp(df, temp):
    '''
    Get monthly only temporlization. The CMV run group only varies by month.
    '''
    scalar_cols = [s_col for s_col in temp.profs.columns if s_col.startswith('Scalar')]
    value_cols = ['qflag',] + scalar_cols
    # Only match by fips/scc or fips; SCC only gets default
    hierarchy = [['region_cd','scc'],['region_cd',],['scc',]]
    df = match_temporal(df, temp.profs, value_cols, hierarchy)
    df['state'] = df['region_cd'].str[:2]
    df.drop(['region_cd','scc'], axis=1, inplace=True)
    df.drop_duplicates(inplace=True)
    if len(df[df['qflag'] != 'MONTH']) > 0:
        print('WARNING: Non-MONTH profiles found. This should not be for CMV.')
    scalar_cols = ['Scalar%s' %x for x in range(1,13)]
    return df[['state','facid','src_id','qflag']+scalar_cols]

def write_temporal(temp):
    '''
    Write the standard monthly temporal profiles "MONTH"
    '''
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','CMV_temporal.csv')
    temp.to_csv(fname, index=False)

def write_emis(emis, source_groups):
    '''
    Write the "AFTER-AERMOD" Emissions files

    These are the HAP emissions by facility and source ID
    '''
    cols = ['state','run_group','facid','src_id','source_group','smoke_name','ann_value']
    emis = pd.merge(emis, source_groups.xref[['scc','source_group']], on='scc', how='left')
    emis['state'] = emis['region_cd'].str[:2]
    emis['run_group'] = 'CMV'
    emis = emis[cols].groupby(['run_group','state','facid','src_id','source_group','smoke_name'], 
      as_index=False).sum()
    fname = os.path.join(os.environ['WORK_PATH'], 'emis', 'CMV_emis.csv')
    emis.to_csv(fname, columns=cols, index=False)

def check_env_vars(evars):
    '''
    Verify that all required environment variables are set
    '''
    for env_var in evars:
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

def get_invtable(fn):
    '''
    Read in the inventory table for the kept pollutants
    See the SMOKE documentation for more information
    '''
    invtable = pd.read_fwf(fn, comment='#', colspecs=[(0,11), (16,32), (41,42), (43,49)], 
      names=['smoke_name','poll','keep','spec_factor'], 
      converters={'smoke_name': str.strip, 'poll': str.strip, 'keep': str.strip})
    invtable = invtable[invtable['keep'].str.upper() == 'Y'].copy()
    invtable.drop_duplicates(inplace=True)
    return invtable

def proc_shapes(shapes, met_grid):
    '''
    Process the shapes into the correct cells and reproject to UTM
    '''
    cell_poly_df = get_met_cell(shapes, met_grid)
    shapes.drop(['col','row','cell'], axis=1, inplace=True)
    shapes = pd.merge(shapes, cell_poly_df, on=['facid','src_id'], how='left')
    shapes['centroid_lon'] = met_grid.colrow_to_ll(shapes['col'].astype('f'), 
      shapes['row'].astype('f'))['lon'].astype('f')
    utm = UTM()
    shapes['utm_zone'] = shapes['centroid_lon'].apply(utm.get_zone)
    shapes[['utmx','utmy']] = shapes[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    shapes['type'] = 'AREAPOLY'
    shapes['state'] = shapes['region_cd'].str[:2]
    return shapes

def main():
    var_list = ('ATPRO_HOURLY','ATPRO_WEEKLY','ATPRO_MONTHLY','ATREF','EMISINV_A','SOURCE_GROUPS',
      'COSTCY','POLY_FILE','INVTABLE','WORK_PATH','GRIDDESC','REGION_IOAPI_GRIDNAME','STATE_FIPS')
    check_env_vars(var_list)
    inv_list = get_inv_list()
    state_fips = os.environ['STATE_FIPS']
    invtable = get_invtable(os.environ['INVTABLE'])
    inv = AnnualFF10(inv_list, invtable, state_fips, use_shapes=True)
    src_groups = SourceGroups(os.environ['SOURCE_GROUPS'])
    # Read in the shapes files and drop ones that aren't needed
    shapes = read_shapes(os.environ['POLY_FILE'])
    emis = label_emis(inv.emis, shapes[['facid','shape_id','area_frac']].drop_duplicates('facid'))
    emis_id = emis[['facid','src_id','region_cd']].drop_duplicates(['facid','src_id'])
    shapes = pd.merge(shapes, emis_id, on='facid', how='left')
    shapes = shapes[shapes['src_id'].notnull()].copy()
    met_grid = Grid(os.environ['REGION_IOAPI_GRIDNAME'], os.environ['GRIDDESC'])
    shapes = proc_shapes(shapes, met_grid)
    write_emis(emis.copy(), src_groups)
    write_locations(shapes.copy())
    write_parameters(shapes) 
    # Temporalization
    scc_list = list(emis['scc'].drop_duplicates())
    fips_list = list(emis['region_cd'].drop_duplicates())
    temp = Temporal(os.environ['ATREF'],os.environ['ATPRO_HOURLY'],os.environ['ATPRO_WEEKLY'],
        os.environ['ATPRO_MONTHLY'], scc_list, fips_list)
    temp_profs = calc_monthly_temp(emis.drop_duplicates(['facid','src_id']), temp)
    write_temporal(temp_profs)

if __name__ == '__main__':
	main()




