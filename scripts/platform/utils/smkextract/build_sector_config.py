#!/nas/longleaf/home/tranhuy/software/pkg/miniconda3/envs/AERMOD/bin/python

import argparse
import json
import os
import re
import sys

SETENV_PATTERN = re.compile(r'^setenv\s+(\S+)\s+"([^"]*)"')


def parse_run_script(script_path):
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"Run script not found: {script_path}")
    sector = None
    inventory_paths = []
    with open(script_path, 'r') as script_file:
        for raw_line in script_file:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            match = SETENV_PATTERN.match(line)
            if not match:
                continue
            key, value = match.groups()
            key = key.strip()
            value = value.strip()
            if key == 'SECTOR':
                sector = value
            elif key.startswith('EMISINV_') or key.startswith('EMISDAY_') or key.startswith('EMISHOUR_'):
                # Treat both inventory and day-specific emission variables the same
                inventory_paths.append(value)
    if sector is None:
        raise ValueError("SECTOR not defined in run script.")
    if not inventory_paths:
        raise ValueError("No EMISINV_* or EMISDAY_* or 'EMISHOUR_*' entries found in run script.")
    unique_files = []
    for path in inventory_paths:
        filename = os.path.basename(path)
        if filename and filename not in unique_files:
            unique_files.append(filename)
    if not unique_files:
        raise ValueError("EMISINV_* entries did not yield any filenames.")
    return sector, unique_files


def update_config(config_path, sector, filenames):
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(config_path, 'r') as cfg_file:
        config = json.load(cfg_file)
    sector_list = config.get('sector')
    if not isinstance(sector_list, dict):
        raise TypeError("Configuration file must contain a 'sector' object.")
    existing = sector_list.get(sector, [])
    if isinstance(existing, list):
        merged = existing + [fn for fn in filenames if fn not in existing]
    else:
        merged = filenames
    sector_list[sector] = merged
    config['sector'] = sector_list
    with open(config_path, 'w') as cfg_file:
        json.dump(config, cfg_file, indent=2)
        cfg_file.write('\n')


def load_runscripts(json_path):
    if not os.path.isfile(json_path):
        raise FileNotFoundError(f"Runscript manifest not found: {json_path}")
    with open(json_path, 'r') as rs_file:
        data = json.load(rs_file)
    if not isinstance(data, dict):
        raise TypeError("Runscript manifest must be a JSON object mapping subsectors to lists of scripts.")
    manifest = {}
    for subsector, scripts in data.items():
        if not isinstance(scripts, list):
            raise TypeError(f"Runscript list for subsector '{subsector}' must be a list.")
        manifest[subsector] = scripts
    return manifest


def main():
    parser = argparse.ArgumentParser(description='Populate smkxtgrid.json sector entries from run scripts.')
    default_config = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'smkxtgrid.json')
    default_runscripts = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'emp_runscripts.json')
    parser.add_argument('--runscripts', default=default_runscripts, help='Path to emp_runscripts.json manifest.')
    parser.add_argument('--config', default=default_config, help='Path to configuration JSON file for smkxtgrid.py')
    args = parser.parse_args()
    try:
        manifest = load_runscripts(args.runscripts)
        updates = 0
        manifest_dir = os.path.dirname(os.path.abspath(args.runscripts))
        for subsector, scripts in manifest.items():
            for script_path in scripts:
                if not os.path.isabs(script_path):
                    script_path = os.path.normpath(os.path.join(manifest_dir, script_path))
                sector, filenames = parse_run_script(script_path)
                update_config(args.config, sector, filenames)
                print(f"Processed subsector '{subsector}' script {script_path} -> sector '{sector}' ({len(filenames)} file(s)).")
                updates += 1
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"Completed updates for {updates} run script(s). Configuration saved to {args.config}.")


if __name__ == '__main__':
    main()
