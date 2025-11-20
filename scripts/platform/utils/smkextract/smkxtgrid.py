#!/nas/longleaf/home/tranhuy/software/pkg/miniconda3/envs/AERMOD/bin/python

# Import the function from griddesc2shp.py
from griddesc2shp import extract_grid
import geopandas as gpd
from shapely.geometry import Polygon
import pyproj
import argparse
import json
import pandas as pd
import os
import csv
import warnings

# SMK Extract
def get_intersecting_counties(domain_gdf, county_shp_path):
	"""
	Returns a list of unique county IDs or names that intersect with the domain_gdf.
	"""
	counties = gpd.read_file(county_shp_path)
	# Ensure both are in the same CRS
	if counties.crs != domain_gdf.crs:
		counties = counties.to_crs(domain_gdf.crs)
	# Find counties that intersect with the domain (any grid cell)
	intersecting = gpd.sjoin(counties, domain_gdf, how="inner", predicate="intersects")
	# Try to find a standard county identifier
	county_id_col = None
	for col in ["GEOID", "geoid", "COUNTYFP", "NAME", "name"]:
		if col in intersecting.columns:
			county_id_col = col
			break
	if county_id_col:
		unique_counties = intersecting[[county_id_col]].drop_duplicates()
		return unique_counties[county_id_col].tolist()
	else:
		# If no standard column, return all index values
		return intersecting.index.tolist()

def create_domain_gdf(griddesc_path, grid_id):
	"""
	Create a GeoPandas dataframe for the specified domain using griddesc_simplified.txt.
	Returns: GeoDataFrame with the full grid polygons and attributes.
	"""
	# Use extract_grid to get projection and grid params
	coord_params, grid_params = extract_grid(griddesc_path, grid_id)
	# Projection parameters
	proj_type, p_alpha, p_beta, p_gamma, x_cent, y_cent = coord_params
	# Grid parameters
	_, xorig, yorig, xcell, ycell, ncols, nrows, _ = grid_params
	# Define the Lambert Conformal projection (spherical earth)
	proj_str = (
		f"+proj=lcc +lat_1={p_alpha} +lat_2={p_beta} +lat_0={y_cent} "
		f"+lon_0={x_cent} +a=6370000.0 +b=6370000.0 +x_0=0 +y_0=0 +units=m +no_defs"
	)
	lcc_proj = pyproj.Proj(proj_str)
	wgs84_proj = pyproj.Proj(proj='latlong', datum='WGS84')
	transformer = pyproj.Transformer.from_proj(lcc_proj, wgs84_proj, always_xy=True)
	# Build polygons for each cell
	features = []
	rows_attr = []
	cols_attr = []
	centers_x = []
	centers_y = []
	centers_lon = []
	centers_lat = []
	for r in range(1, int(nrows)+1):
		y0 = yorig + (r-1) * ycell
		y1 = y0 + ycell
		for c in range(1, int(ncols)+1):
			x0 = xorig + (c-1) * xcell
			x1 = x0 + xcell
			pts_proj = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
			pts_ll = [transformer.transform(px, py) for (px, py) in pts_proj]
			poly = Polygon(pts_ll + [pts_ll[0]])
			features.append(poly)
			rows_attr.append(r)
			cols_attr.append(c)
			cx = x0 + 0.5 * xcell
			cy = y0 + 0.5 * ycell
			clon, clat = transformer.transform(cx, cy)
			centers_x.append(cx)
			centers_y.append(cy)
			centers_lon.append(clon)
			centers_lat.append(clat)
	gdf = gpd.GeoDataFrame(
		{
			'name': [grid_id]*len(features),
			'ROWNUM': rows_attr,
			'COLNUM': cols_attr,
			'center_x_m': centers_x,
			'center_y_m': centers_y,
			'center_lon': centers_lon,
			'center_lat': centers_lat
		},
		geometry=features,
		crs='EPSG:4326'
	)
	return gdf

def get_ff10_type(path):
	"""
	Determines the FF10 file type based on the identifier in the first line.
	Returns: string indicating FF10 type.
	"""
	# Default headers if not found in file
	default_headers = {
		'FF10_POINT':['country_cd', 'region_cd', 'tribal_code', 'facility_id', 'unit_id', 'rel_point_id', 'process_id', 'agy_facility_id', 'agy_unit_id', 'agy_rel_point_id', 'agy_process_id', 'scc', 'poll', 'ann_value', 'ann_pct_red', 'facility_name', 'erptype', 'stkhgt', 'stkdiam', 'stktemp', 'stkflow', 'stkvel', 'naics', 'longitude', 'latitude', 'll_datum', 'horiz_coll_mthd', 'design_capacity', 'design_capacity_units', 'reg_codes', 'fac_source_type', 'unit_type_code', 'control_ids', 'control_measures', 'current_cost', 'cumulative_cost', 'projection_factor', 'submitter_id', 'calc_method', 'data_set_id', 'facil_category_code', 'oris_facility_code', 'oris_boiler_id', 'ipm_yn', 'calc_year', 'date_updated', 'fug_height', 'fug_width_xdim', 'fug_length_ydim', 'fug_angle', 'zipcode', 'annual_avg_hours_per_year', 'jan_value', 'feb_value', 'mar_value', 'apr_value', 'may_value', 'jun_value', 'jul_value', 'aug_value', 'sep_value', 'oct_value', 'nov_value', 'dec_value',"jan_pctred","feb_pctred","mar_pctred","apr_pctred","may_pctred","jun_pctred","jul_pctred","aug_pctred","sep_pctred","oct_pctred","nov_pctred","dec_pctred","comment"],
		'FF10_NONPOINT':['country_cd', 'region_cd', 'tribal_code', 'census_tract_cd', 'shape_id', 'scc', 'emis_type', 'poll', 'ann_value', 'ann_pct_red', 'control_ids', 'control_measures', 'current_cost', 'cumulative_cost', 'projection_factor', 'reg_codes', 'calc_method', 'calc_year', 'date_updated', 'data_set_id', 'jan_value', 'feb_value', 'mar_value', 'apr_value', 'may_value', 'jun_value', 'jul_value', 'aug_value', 'sep_value', 'oct_value', 'nov_value', 'dec_value', 'jan_pctred','feb_pctred','mar_pctred','apr_pctred','may_pctred','jun_pctred','jul_pctred','aug_pctred','sep_pctred','oct_pctred','nov_pctred','dec_pctred','comment'],
		'FF10_ACTIVITY':['country_cd', 'region_cd', 'tribal_code', 'census_tract_cd', 'shape_id', 'scc', 'CD', 'MSR', 'activity_type', 'ann_parm_value', 'calc_year', 'date_updated', 'data_set_id', 'jan_value', 'feb_value', 'mar_value', 'apr_value', 'may_value', 'jun_value', 'jul_value', 'aug_value', 'sep_value', 'oct_value', 'nov_value', 'dec_value', 	'comment'],
		'FF10_HOURLY_POINT':['country_cd', 'region_cd', 'tribal_code', 'facility_id', 'unit_id', 'rel_point_id', 'process_id', 'scc', 'poll', 'op_type_cd', 'calc_method', 'date_updated', 'date', 'daytot', 'hrval0', 'hrval1', 'hrval2', 'hrval3', 'hrval4', 'hrval5', 'hrval6', 'hrval7', 'hrval8', 'hrval9', 'hrval10', 'hrval11', 'hrval12', 'hrval13', 	'comment'],
		'FF10_NONROAD':['country_cd', 'region_cd', 'tribal_code', 'census_tract_cd', 'shape_id', 'scc', 'emis_type', 'poll', 'ann_value', 'ann_pct_red', 'control_ids', 'control_measures', 'current_cost', 'cumulative_cost', 'projection_factor', 'reg_codes', 'calc_method', 'calc_year', 'date_updated', 'data_set_id', 'jan_value', 'feb_value', 'mar_value', 'apr_value', 'may_value', 'jun_value', 'jul_value', 'aug_value', 'sep_value', 'oct_value', 'nov_value', 'dec_value','jan_pctred','feb_pctred','mar_pctred','apr_pctred','may_pctred','jun_pctred','jul_pctred','aug_pctred','sep_pctred','oct_pctred','nov_pctred','dec_pctred','comment'],
		'FF10_DAILY_POINT':['country_cd','region_cd','tribal_code','facility_id','unit_id','rel_point_id','process_id','scc','poll','op_type_cd','calc_method','date_updated','monthnum','monthtot','dayval1','dayval2','dayval3','dayval4','dayval5','dayval6','dayval7','dayval8','dayval9','dayval10','dayval11','dayval12','dayval13','dayval14','dayval15','dayval16','dayval17','dayval18','dayval19','dayval20','dayval21','dayval22','dayval23','dayval24','dayval25','dayval26','dayval27','dayval28','dayval29','dayval30','dayval31','comment']
	}

	with open(path, 'r') as f:
		first_line = f.readline()
		if 'FF10_POINT' in first_line:
			ff10_fmt = 'FF10_POINT'
		elif 'FF10_NONPOINT' in first_line:
			ff10_fmt = 'FF10_NONPOINT'
		elif 'FF10_ACTIVITY' in first_line:
			ff10_fmt = 'FF10_ACTIVITY'
		elif 'FF10_HOURLY_POINT' in first_line:
			ff10_fmt = 'FF10_HOURLY_POINT'
		elif 'FF10_NONROAD' in first_line:
			ff10_fmt = 'FF10_NONROAD'
		elif 'FF10_DAILY_POINT' in first_line:
			ff10_fmt = 'FF10_DAILY_POINT'
		else:
			warnings.warn("Unknown FF10 file format identifier in first line.")
			ff10_fmt = None

		if ff10_fmt is not None:
			expected_headers = default_headers[ff10_fmt]
		# Reset file pointer to beginning
		f.seek(0)
	
	# Find index of "region_cd" column in header
	for icol, col_name in enumerate(expected_headers):
		if col_name in ["region_cd","region","fips"]:
			region_cd_idx = icol
			break

	return ff10_fmt, expected_headers, region_cd_idx

def read_ff10(path):
	"""
	Reads a FF10_NONPOINT CSV file, skipping comment lines and using the correct header.
	Returns: pandas DataFrame
	"""
	"""
	Reads a FF10_NONPOINT CSV file, skipping comment lines and using the first non-comment line as header.
	Returns: pandas DataFrame
	"""
	
	ff10_fmt, expected_headers, expected_region_idx = get_ff10_type(path)

	header_line = None
	region_cd_idx = None

	with open(path, 'r') as f:
		for i, line in enumerate(f):
			if not line.lstrip().startswith('#') and line.strip() != '' and 'region' in line:
				header_line = i
				# Check if found header matches expected headers
				found_headers = [h.strip() for h in line.strip().split(',')]
				if found_headers != expected_headers:
					# Raise warning if headers do not match expected, but still use found headers
					warnings.warn("Header mismatch in FF10 file.")
				#print(f"Header line identified: {line.strip()}", flush=True)
				break
	if header_line is not None:
		# Read the CSV file using the identified header line
		df = pd.read_csv(
			path,
			skiprows=header_line,
			header=0,
			comment='#',
			dtype=object,
			low_memory=False,
		)
	if header_line is None:
		warnings.warn(f"No header line found in file. Using default headers for file format: {ff10_fmt}")
		# Used default headers
		df = pd.read_csv(
			path,
			names=expected_headers,
			header=None,
			comment='#',
			dtype=object,
			low_memory=False,
		)

	# Skip metadata rows so pandas uses the detected header row for column names
	# Find index of "region_cd" column in header
	for icol, col_name in enumerate(df.columns):
		if col_name in ["region_cd","region","fips"]:
			region_cd_idx = icol
			break
	
	if region_cd_idx is None:
		warnings.warn("No region_cd column found in FF10 file.")
		region_cd_idx = expected_region_idx
	
	# Attach region_cd_idx as an attribute for later use
	df.region_cd_idx = region_cd_idx
	return df

def filter_ff10_by_counties(input_path, output_path, region_cd_list):
	"""
	Writes a filtered FF10 file with only rows where region_cd is in region_cd_list.
	Retains all comment lines and the original header.
	Adds a comment line with timestamp and description at the top.
	"""
	import datetime
	ff10_fmt, expected_headers, expected_region_idx = get_ff10_type(input_path)
	header_line = None
	today = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	comment_line = f"# Extracted from {input_path} on {today} by filtering region_cd\n"
	with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
		comment_lines = []
		data_lines = []
		header_line = None
		region_cd_idx = None
		for line in fin:
			if line.lstrip().startswith('#') or line.strip() == '':
				comment_lines.append(line)
				continue

			if header_line is None:
				header_line = line
				header = [h.strip() for h in line.strip().split(',')]
				for icol, col_name in enumerate(header):
					if col_name in ["region_cd","region","fips"]:
						region_cd_idx = icol
						expected_region_idx = region_cd_idx
						break
				if region_cd_idx is None:
					header_line = None

			# Only write data lines where region_cd is in the list
			fields = line.strip().split(',')
			if len(fields) > expected_region_idx:
				rc = fields[expected_region_idx].strip('"')
				if rc in region_cd_list:
					data_lines.append(line)

		# Write all comment lines first
		for cl in comment_lines:
			fout.write(cl)
		# Write the new comment line after all other comments
		fout.write(comment_line)
		# Write header and data
		if header_line is not None:
			fout.write(header_line)
		for dl in data_lines:
			fout.write(dl)


def load_config(config_path):
	"""
	Load configuration values from a JSON file.
	Returns: dict with resolved paths for inputs/outputs/county_shp/griddesc and grid/ff10 settings.
	"""
	if not os.path.isfile(config_path):
		raise FileNotFoundError(f"Configuration file not found: {config_path}")
	with open(config_path, 'r') as cfg_file:
		config = json.load(cfg_file)
	required_keys = ['inputs', 'outputs', 'sector', 'county_shp', 'griddesc_path', 'grid_id']
	missing = [key for key in required_keys if key not in config]
	if missing:
		raise KeyError(f"Missing required key(s) in configuration: {', '.join(missing)}")
	inputs_dir = config['inputs']
	outputs_dir = config['outputs']
	sector_mapping = config['sector']
	sector_skip = config.get('sector_skip', [])
	county_shp_path = config['county_shp']
	griddesc_path = config['griddesc_path']
	grid_id = config['grid_id']
	for key, value in [('inputs', inputs_dir), ('outputs', outputs_dir), ('county_shp', county_shp_path), ('griddesc_path', griddesc_path), ('grid_id', grid_id)]:
		if not isinstance(value, str):
			raise TypeError(f"Configuration key '{key}' must be a string.")
	if not isinstance(sector_mapping, dict):
		raise TypeError("Configuration key 'sector' must map sectors to lists of filenames.")
	for sector, files in sector_mapping.items():
		if not isinstance(files, list):
			raise TypeError(f"Configuration entry for sector '{sector}' must be a list of filenames.")
	if not isinstance(sector_skip, list) or not all(isinstance(item, str) for item in sector_skip):
		raise TypeError("Configuration key 'sector_skip' must be a list of strings when provided.")
	config_dir = os.path.dirname(os.path.abspath(config_path))
	if not os.path.isabs(inputs_dir):
		inputs_dir = os.path.normpath(os.path.join(config_dir, inputs_dir))
	if not os.path.isabs(outputs_dir):
		outputs_dir = os.path.normpath(os.path.join(config_dir, outputs_dir))
	if not os.path.isabs(county_shp_path):
		county_shp_path = os.path.normpath(os.path.join(config_dir, county_shp_path))
	if not os.path.isabs(griddesc_path):
		griddesc_path = os.path.normpath(os.path.join(config_dir, griddesc_path))
	return {
		'inputs': inputs_dir,
		'outputs': outputs_dir,
		'sector': sector_mapping,
		'sector_skip': sector_skip,
		'county_shp': county_shp_path,
		'griddesc_path': griddesc_path,
		'grid_id': grid_id,
	}

if __name__ == "__main__":
	default_config_name = "smkxtgrid.json"
	parser = argparse.ArgumentParser(description="Filter FF10 inventories to counties intersecting a target grid.")
	parser.add_argument("--config", default=None, help="Path to JSON configuration with required inputs/outputs/sector/grid settings.")
	args = parser.parse_args()
	if args.config:
		config_path = args.config
	else:
		script_dir = os.path.dirname(os.path.abspath(__file__))
		config_path = os.path.join(script_dir, default_config_name)
	config = load_config(config_path)
	county_shp = config['county_shp']
	griddesc_path = config['griddesc_path']
	grid_id = config['grid_id']
	inputs = config['inputs']
	outputs = config['outputs']
	sector_mapping = config['sector']
	skip_sectors = set(config.get('sector_skip', []))
	print(f"Using configuration {config_path}", flush=True)
	print(f"Inputs directory resolved to {inputs}", flush=True)
	print(f"Outputs directory resolved to {outputs}", flush=True)
	print(f"County shapefile resolved to {county_shp}", flush=True)
	print(f"griddesc path resolved to {griddesc_path}", flush=True)
	print(f"Grid ID set to {grid_id}", flush=True)
	os.makedirs(outputs, exist_ok=True)

	print(f"Creating domain {grid_id}...\n", flush=True)
	domain = create_domain_gdf(griddesc_path, grid_id)
	print(domain.head(), flush=True)

	print("\nFinding intersecting counties...", flush=True)
	intersecting_counties = get_intersecting_counties(domain, county_shp)
	print(f"\nCounties intersecting the domain:")
	print(intersecting_counties, flush=True)

	# Read emission input data to pandas DataFrame
	for sector,fnames in sector_mapping.items():
		if sector in skip_sectors:
			print(f"Skipping sector '{sector}' per configuration.", flush=True)
			continue
		sector_output = os.path.join(outputs, sector)
		os.makedirs(sector_output, exist_ok=True)
		for fname in fnames:
			ff10_path = os.path.join(inputs, sector, fname)
			print(f"\nReading emission input data from {ff10_path}...", flush=True)
			emis_df = read_ff10(ff10_path)
			print(emis_df.head(), flush=True)
			ff10_counties = list(emis_df['region_cd'].unique())
			#print(ff10_counties)

			# Filter ff10 file to only intersecting counties
			ff10_base = os.path.splitext(os.path.basename(ff10_path))[0]
			ff10_ext = os.path.splitext(ff10_path)[1]
			output_ff10 = os.path.join(sector_output, f"{ff10_base}_{grid_id}{ff10_ext}")
			print(f"\nWriting filtered FF10 file to {output_ff10}...", flush=True)
			filter_ff10_by_counties(ff10_path, output_ff10, [str(rc) for rc in intersecting_counties])
			print("Done.", flush=True)
