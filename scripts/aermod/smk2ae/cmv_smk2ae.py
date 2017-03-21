#!/usr/bin/env python

import os, sys, string
from collections import Counter
from osgeo import ogr
import pandas as pd
import smk2ae
from smk2ae.grid import Grid
from smk2ae.utm import UTM

def read_emissions(emis_list, haps_list, state_fips=False):
    '''
    Read in the ACCESS exported emissions allocation input
    These files should be converted into a standard condensed format where
     region_cd, poll, scc, facid, and ann_value are all present as columns.
    '''
    emis = pd.DataFrame()
    haps = pd.read_csv(haps_list,usecols=['poll','smoke_name','spec_factor','ure'],
      dtype={'poll':'|S16'})
    for emis_name in emis_list:
        df = pd.read_csv(emis_name, usecols=['facid','scc','poll','ann_value'],
          dtype={'poll': '|S16', 'scc': '|S10', 'ann_value': 'f'})
        df['region_cd'] = df['facid'].apply(lambda x: x.split('F')[1][:5])
        if state_fips:
            df = df[df['region_cd'].str.startswith('%0.2d' %int(state_fips))].copy()
        df['region_cd'] = df['region_cd'].str.zfill(5)
        df = pd.merge(df, haps, on='poll', how='left')
        df = df[df['spec_factor'].notnull()].copy()
        df['ann_value'] = df['ann_value'] * df['spec_factor']
        df.drop(['poll','spec_factor'], axis=1, inplace=True)
        df.ix[df['scc'].isin(('2280002100','2280002200')), 'scc'] = 'c1c2'
        df.ix[df['scc'].isin(('2280003100','2280003200')), 'scc'] = 'c3'
        df = df.groupby(['region_cd','facid','scc','smoke_name','ure'], as_index=False).sum()
        emis = pd.concat((emis, df))
    emis['src_id'] = emis['facid'].apply(name_polygon)
    emis['facid'] = emis[['facid','region_cd']].apply(lambda row: '%sF%s' %(row['facid'].split('F')[0],
      row['region_cd']), axis=1)
    return emis

def name_polygon(facid):
    '''
    Format the source ID name based on the polygon number
    '''
    if facid[0] == 'U':
        return facid[:5] + 'P' + facid.split('P')[1]
    else:
        return facid

def read_shapes(shapes_file, st_fips):
    '''
    Read in the shapes from the shapes file

    The shapes file is formatted like this:
    facid,sourceid,loctype,lon,lat,zone,numvert,area
    U0612F10001P001,TRACT,L,-75.103862,38.995806,,45,586792801
    U0612F10001P001,TRACT,L,-75.311905,38.945821,,45,586792801
    '''
    df = pd.read_csv(shapes_file, usecols=['facid','lon','lat','numvert','area'])
    df['region_cd'] = df['facid'].apply(lambda x: x.split('F')[1][:5])
    if st_fips:
        df = df[df['region_cd'].str.startswith(st_fips)].copy()
    df['region_cd'] = df['region_cd'].str.zfill(5)
    df['src_id'] = df['facid'].apply(name_polygon)
    df['facid'] = df[['facid','region_cd']].apply(lambda row: '%sF%s' %(row['facid'].split('F')[0],
      row['region_cd']), axis=1)
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
        for cell, cell_poly in cell_polys.iteritems():
            overlap_area = src_poly.Intersection(cell_poly).GetArea()
            if  overlap_area > max_area:
                max_area = overlap_area
                dom_cell = cell
    else:
        #raise AttributeError,'Shape %s %s contains an invalid polygon' %(facid, src_id)
        print 'WARNING: Shape %s %s contains an invalid polygon' %(facid, src_id)
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

def write_locations(type_df, source_type):
    '''
    Write the locations file
    '''
    out_cols = ['state','region_cd','facid','src_id','type','col','row','utmx','utmy','utm_zone',
      'lon','lat']
    fname = os.path.join(os.environ['WORK_PATH'], 'locations', '%s_locations.csv' %source_type)
    type_df.sort().to_csv(fname, columns=out_cols, index=False)

def calc_state_rh(emis):
    '''
    Calculate the state release height from the emissions

    The release height is calculated for an entire state based on the function:
      8.4*(sum of c1&c2 HAPS)/(sum of all CMV HAPS) + 20*(sum of c3 HAPS)/(sum of all CMV HAPS)

    Although the release height is the same for the entire state, it is returned as a dataframe
      of release heights by full state and county level FIPS
    '''
    emis['ann_value'] = emis['ann_value'] * emis['ure']
    emis['state'] = emis['region_cd'].str[:2]
    c1c2 = emis.ix[emis['scc'] == 'c1c2', ['state','ann_value']].groupby('state').sum()
    c3 = emis.ix[emis['scc'] == 'c3', ['state','ann_value']].groupby('state').sum()
    rh = c1c2.join(c3, how='outer', lsuffix='_c1c2', rsuffix='_c3').reset_index()
    rh.fillna(0, inplace=True)
    rh['sum'] = rh['ann_value_c1c2'] + rh['ann_value_c3']
    rh['rh'] = (8.4 * rh['ann_value_c1c2']/rh['sum']) + (20. * rh['ann_value_c3']/rh['sum'])
    return rh[['state','rh']].copy().drop_duplicates()

def write_parameters(type_df, state_rh, source_type):
    '''
    Write the parameters file
    
    The maximum number of vertices for a polygon needs to be calculated first
    The coordinates for each polygon is written to a single line
    '''
    param_df = pd.merge(type_df, state_rh, on='state', how='left')
    param_df['sz'] = param_df['rh']/2.15 
    fname = os.path.join(os.environ['WORK_PATH'], 'parameters', '%s_area_params.csv' %source_type)
    maxverts = param_df['numvert'].max()
    coord_cols = []
    [coord_cols.extend(['utmx','utmy']) for vert_num in xrange(1, maxverts + 1)]
    ll_cols = []
    [ll_cols.extend(['lon','lat']) for vert_num in xrange(1, maxverts + 1)]
    out_cols = ['state','region_cd','facid','src_id','type','area','rh','numvert','sz'] + coord_cols\
      + ll_cols
    with open(fname,'w') as f:
        f.write('%s\n' %','.join(out_cols))
        for facid in list(param_df['facid'].drop_duplicates()):
            for src_id in list(param_df.ix[param_df['facid'] == facid, 'src_id'].drop_duplicates()):
                src_df = param_df[(param_df['facid']==facid) & (param_df['src_id']==src_id)].copy()
                numverts = src_df['numvert'].values[0]
                out_line = [src_df['state'].values[0], src_df['region_cd'].values[0], facid, 
                  src_id, src_df['type'].values[0], str(src_df['area'].values[0]),
                  str(src_df['rh'].values[0]), str(numverts),str(src_df['sz'].values[0])]
                for ix, row in src_df.iterrows():
                    out_line.extend([str(row['utmx']),str(row['utmy'])])
                if numverts < maxverts:
                    for x in xrange(maxverts - numverts):
                        out_line.extend(['',''])
                for ix, row in src_df.iterrows():
                    out_line.extend([str(row['lon']),str(row['lat'])])
                if numverts < maxverts:
                    for x in xrange(maxverts - numverts):
                        out_line.extend(['',''])
                f.write('%s\n' %','.join(out_line))
        
def write_temporal(emis, type_df, source_type):
    '''
    Write the standard non-hourly non-daily profile based temporal profiles
    The profiles are based on the toxicity weighted emissions.
    '''
    temp_df = pd.merge(type_df[['state','facid','src_id']], emis[['facid','src_id','scc',
      'ann_value','ure']], on=['facid','src_id'], how='left')
    temp_df['ann_value'] = temp_df['ann_value'] * temp_df['ure']
    c1c2 = temp_df.ix[temp_df['scc'] == 'c1c2', ['state','ann_value']].groupby('state').sum()
    c3 = temp_df.ix[temp_df['scc'] == 'c3', ['state','ann_value']].groupby('state').sum()
    scalars = c1c2.join(c3, how='outer', lsuffix='_c1c2', rsuffix='_c3').reset_index()
    scalars.fillna(0, inplace=True)
    # These are provided c3 monthly factors
    c3_facs = (0.077, 0.072, 0.077, 0.08, 0.086, 0.087, 0.09, 0.091, 0.085, 0.084, 0.084, 0.087)
    for scalar_num, c3_factor in enumerate(c3_facs):
        scalars['Scalar%s' %(scalar_num + 1)] = 12. * (((scalars['ann_value_c1c2'] * (1/12.)) + \
          (scalars['ann_value_c3'] * c3_factor)) / (scalars['ann_value_c1c2'] + scalars['ann_value_c3']))
    temp_df = pd.merge(temp_df[['state','facid','src_id']].drop_duplicates(), scalars, on='state', 
      how='left')
    temp_df['qflag'] = 'MONTH'
    fname = os.path.join(os.environ['WORK_PATH'],'temporal','%s_temporal.csv' %source_type)
    cols = ['state','facid','src_id','qflag'] + ['Scalar%s' %x for x in xrange(1,13)]
    temp_df.to_csv(fname, index=False, columns=cols)

def write_emis(emis, source_type):
    '''
    Write the "AFTER-AERMOD" Emissions files

    These are the HAP emissions by facility and source ID
    '''
    emis['state'] = emis['region_cd'].str[:2]
    emis = emis[['state','facid','src_id','smoke_name','ann_value']].copy().groupby(['state',
      'facid','src_id','smoke_name'], as_index=False).sum()
    emis['run_group'] = 'CMV'
    emis['source_group'] = source_type
    fname = os.path.join(os.environ['WORK_PATH'], 'emis', '%s_emis.csv' %source_type)
    emis.to_csv(fname, columns=['state','run_group','facid','src_id','source_group','smoke_name',
      'ann_value'], index=False)

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

def main():
    var_list = ('GRIDDESC','REGION_IOAPI_GRIDNAME','EMISINV_A','HAPS_LIST','BASE_YEAR',
      'SECTOR','WORK_PATH','POLY_FILE','STATE_FIPS')
    check_env_vars(var_list)
    inv_list = get_inv_list()
    state_fips = os.environ['STATE_FIPS']
    emis = read_emissions(inv_list, os.environ['HAPS_LIST'], state_fips)
    shapes = read_shapes(os.environ['POLY_FILE'], state_fips)
    # Drop shapes that have no emissions
    shapes = pd.merge(shapes, emis[['facid','src_id','scc']].drop_duplicates(['facid','src_id']), 
      on=['facid','src_id'], how='left')
    shapes = shapes[shapes['scc'].notnull()].copy()
    shapes.drop('scc', axis=1, inplace=True)
    met_grid = Grid(os.environ['REGION_IOAPI_GRIDNAME'], os.environ['GRIDDESC'])
    cell_poly_df = get_met_cell(shapes, met_grid)
    shapes.drop(['col','row','cell'],axis=1,inplace=True)
    shapes = pd.merge(shapes, cell_poly_df, on=['facid','src_id'], how='left')
    shapes['centroid_lon'] = met_grid.colrow_to_ll(shapes['col'].astype('f') + 0.5, 
      shapes['row'].astype('f') + 0.5)['lon'].astype('f')
    utm = UTM()
    shapes['utm_zone'] = shapes['centroid_lon'].apply(utm.get_zone)
    shapes[['utmx','utmy']] = shapes[['lon','lat','utm_zone']].apply(lambda x: \
      pd.Series(utm.get_coords(x['lon'],x['lat'],x['utm_zone'])), axis=1)
    shapes['type'] = 'AREAPOLY'
    state_rh = calc_state_rh(emis[['region_cd','scc','ann_value','ure']].copy())
    '''
    Ports and underway get separate output files except for the emissions files
    '''
    shapes['state'] = shapes['region_cd'].str[:2]
    for source_type in ('ports','underway'):
        print 'Running for %s' %source_type
        shapes_type = shapes[shapes['facid'].str.startswith(source_type[0].upper())].copy()
        if shapes_type.empty:
            print 'WARNING: No %s sources found' %source_type
        else:
            emis_type = emis[emis['facid'].str.startswith(source_type[0].upper())].copy()
            write_emis(emis_type.copy(), source_type)
            write_locations(shapes_type.copy().drop_duplicates(['region_cd','facid','src_id',
              'type','col','row','utm_zone']), source_type)
            write_parameters(shapes_type, state_rh, source_type) 
            write_temporal(emis_type, shapes_type, source_type)

if __name__ == '__main__':
	main()




