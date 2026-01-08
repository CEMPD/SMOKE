#!/nas/longleaf/home/tranhuy/software/pkg/miniconda3/envs/AERMOD/bin/python

import argparse
import json
import os
import csv
import warnings

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
		'FF10_DAILY_POINT':['country_cd','region_cd','tribal_code','facility_id','unit_id','rel_point_id','process_id','scc','poll','op_type_cd','calc_method','date_updated','monthnum','monthtot','dayval1','dayval2','dayval3','dayval4','dayval5','dayval6','dayval7','dayval8','dayval9','dayval10','dayval11','dayval12','dayval13','dayval14','dayval15','dayval16','dayval17','dayval18','dayval19','dayval20','dayval21','dayval22','dayval23','dayval24','dayval25','dayval26','dayval27','dayval28','dayval29','dayval30','dayval31','comment'],
		'FF10_DAILY_NONPOINT':['country_cd','region_cd','tribal_code','census_tract','shape_id','tbd','emis_type','scc','poll','op_type_cd','calc_method','date_updated','monthnum','monthtot','dayval1','dayval2','dayval3','dayval4','dayval5','dayval6','dayval7','dayval8','dayval9','dayval10','dayval11','dayval12','dayval13','dayval14','dayval15','dayval16','dayval17','dayval18','dayval19','dayval20','dayval21','dayval22','dayval23','dayval24','dayval25','dayval26','dayval27','dayval28','dayval29','dayval30','dayval31','comment']
	}

	# Default return if unknown
	expected_headers = None

	with open(path, 'r') as f:
		first_line = f.readline()
		# Order is important here
		if 'FF10_DAILY_POINT' in first_line or 'FF10_POINT_DAILY' in first_line:
			ff10_fmt = 'FF10_DAILY_POINT'
		elif 'FF10_DAILY_NONPOINT' in first_line:
			ff10_fmt = 'FF10_DAILY_NONPOINT'			
		elif 'FF10_ACTIVITY' in first_line:
			ff10_fmt = 'FF10_ACTIVITY'
		elif 'FF10_HOURLY_POINT' in first_line:
			ff10_fmt = 'FF10_HOURLY_POINT'
		elif 'FF10_NONROAD' in first_line:
			ff10_fmt = 'FF10_NONROAD'
		elif 'FF10_POINT' in first_line:
			ff10_fmt = 'FF10_POINT'		
		elif 'FF10_NONPOINT' in first_line:
			ff10_fmt = 'FF10_NONPOINT'				
		else:
			warnings.warn("Unknown FF10 file format identifier in first line.")
			ff10_fmt = None

		if ff10_fmt is not None:
			expected_headers = default_headers[ff10_fmt]
		# Reset file pointer to beginning
		f.seek(0)

	return ff10_fmt, expected_headers

def parse_ff10_file(file_path):
	"""
	Parse an FF10 file preserving raw lines while providing parsed values for logic.

	Returns a dict with:
	- 'fmt': FF10 type string or None
	- 'expected_headers': list of expected column names or None
	- 'comments': list[str] of comment lines (with no trailing newlines)
	- 'header': raw header line as a string (no trailing newline) or None
	- 'header_cols': list[str] of column names used to interpret the data rows
	- 'rows': list of tuples (raw_line:str without trailing newline, values:list[str]) for each data row
	"""
	fmt, expected_headers = get_ff10_type(file_path)
	comments = []
	header_raw = None
	header_cols = None
	rows = []

	def split_csv_line(line):
		# Use Python CSV reader to properly handle commas/quotes in a single line
		return next(csv.reader([line], skipinitialspace=False))

	with open(file_path, 'r') as f:
		lines = f.readlines()

	# Iterate and capture comments until first non-comment/non-empty
	i = 0
	while i < len(lines):
		line = lines[i].rstrip('\n')
		if line.strip() == '':
			i += 1
			continue
		if line.lstrip().startswith('#'):
			comments.append(line)
			i += 1
			continue
		# First non-comment, non-empty line
		cand_vals = [c.strip() for c in split_csv_line(line)]
		cand_lower = [c.lower() for c in cand_vals]
		# Heuristic: if looks like header (contains 'region_cd' column name), treat as header
		if 'region_cd' in cand_lower:
			header_raw = line
			header_cols = cand_vals
			i += 1
		else:
			# No explicit header row; use expected headers if available
			header_cols = expected_headers if expected_headers is not None else None
			# Do NOT consume this line as header; treat as data
			# fall through without incrementing i here; we'll parse data rows next
			pass
		break

	# If we broke out without encountering any non-comment lines
	if i >= len(lines):
		return {
			'fmt': fmt,
			'expected_headers': expected_headers,
			'comments': comments,
			'header': header_raw,
			'header_cols': header_cols,
			'rows': rows
		}

	# Process remaining lines as data rows
	while i < len(lines):
		line = lines[i].rstrip('\n')
		if line.strip() == '':
			i += 1
			continue
		if line.lstrip().startswith('#'):
			# Comments in the middle should be preserved if we ever need, but
			# treat as comments not data. Append to comments to ensure they are retained.
			comments.append(line)
			i += 1
			continue
		values = split_csv_line(line)
		rows.append((line, values))
		i += 1

	return {
		'fmt': fmt,
		'expected_headers': expected_headers,
		'comments': comments,
		'header': header_raw,
		'header_cols': header_cols,
		'rows': rows
	}

def filter_file_by_col(file_path, config):
	"""
	Filter rows in file based on column match.
	"""
	parsed = parse_ff10_file(file_path)
	header_cols = parsed['header_cols']
	
	if not header_cols:
		warnings.warn(f"No header cols found for {file_path}, cannot filter by column.")
		return {'comments': parsed['comments'], 'header': parsed['header'], 'lines': []}
	
	col = config.get('filter_col')
	start_val = config.get('start_val')
	end_val = config.get('end_val')
	if isinstance(config.get('filtered_val'), list):
		filtered_val = set(config.get('filtered_val'))
	else:
		filtered_val = None

	name_to_idx = {name.lower(): idx for idx, name in enumerate(header_cols)}
	col_idx = name_to_idx.get(col.lower()) if col else None

	if col_idx is None:
		warnings.warn(f"Column '{col}' not found in {file_path}. Available: {header_cols}")
		return {'comments': parsed['comments'], 'header': parsed['header'], 'lines': []}

	lines_out = []
	for raw_line, values in parsed['rows']:
		if col_idx >= len(values):
			continue
		
		val = values[col_idx]
		keep = False

		# Discrete check
		if filtered_val:
			if val in filtered_val:
				keep = True
		
		# Range check
		if not keep and (start_val is not None and end_val is not None):
			# Try float comparison
			try:
				fval = float(val)
				sval = float(start_val)
				eval_ = float(end_val)
				if sval <= fval <= eval_:
					keep = True
			except ValueError:
				# String comparison
				if str(start_val) <= val <= str(end_val):
					keep = True
		
		if keep:
			lines_out.append(raw_line)

	return {
		'comments': parsed['comments'],
		'header': parsed['header'],
		'lines': lines_out
	}

def load_config(config_path):
	if not os.path.isfile(config_path):
		raise FileNotFoundError(f"Config not found: {config_path}")
	with open(config_path, 'r') as f:
		config = json.load(f)
	
	inputs_dir = config.get('inputs', '')
	outputs_dir = config.get('outputs', '')
	config_dir = os.path.dirname(os.path.abspath(config_path))

	# Handle sector-based file selection
	files = []
	if 'sector' in config and isinstance(config['sector'], dict):
		sector_skip = set(config.get('sector_skip', []))
		for sector, file_list in config['sector'].items():
			if sector in sector_skip:
				continue
			if isinstance(file_list, list):
				# Append sector subdirectory to filenames
				files.extend([os.path.join(sector, f) for f in file_list])
	elif 'files' in config:
		# Fallback to old list format if present
		files = config['files']

	def resolve_path(p):
		if os.path.isabs(p): return p
		base = os.path.join(config_dir, inputs_dir, p) if inputs_dir else os.path.join(config_dir, p)
		return os.path.normpath(base)

	resolved_files = [resolve_path(f) for f in files]
	
	if outputs_dir and not os.path.isabs(outputs_dir):
		outputs_dir = os.path.normpath(os.path.join(config_dir, outputs_dir))
	
	config['files'] = resolved_files
	config['outputs'] = outputs_dir
	return config

if __name__ == "__main__":
	default_config_name = "smkxtcol.json"
	parser = argparse.ArgumentParser(description="Filter FF10 files by column values.")
	parser.add_argument("--config", default=None, help="Path to config file")
	args = parser.parse_args()

	if args.config:
		config_path = args.config
	else:
		script_dir = os.path.dirname(os.path.abspath(__file__))
		config_path = os.path.join(script_dir, default_config_name)

	config = load_config(config_path)
	outdir = config['outputs']
	os.makedirs(outdir, exist_ok=True)

	for fpath in config['files']:
		res = filter_file_by_col(fpath, config)
		
		base = os.path.basename(fpath)
		outfile = os.path.join(outdir, f"filtered_{base}")
		
		with open(outfile, 'w') as f:
			# Write original comments
			for c in res['comments']:
				f.write(c + '\n')
			col = config.get('filter_col')
			if config.get('filtered_val'):
				filter_info = f"values: {config.get('filtered_val')}"
			elif config.get('start_val') is not None and config.get('end_val') is not None:
				filter_info = f"range: {config.get('start_val')} to {config.get('end_val')}"
			else:
				filter_info = "unknown criteria"
			
			f.write(f"# Filtered by column '{col}' with {filter_info}\n")
			
			if res['header']:
				f.write(res['header'] + '\n')
			
			for line in res['lines']:
				f.write(line + '\n')
		
		print(f"Wrote {len(res['lines'])} lines to {outfile}")
