# smkscript.py

## Overview
`smkscript.py` is a utility script for processing and grouping environment variable assignments in tcsh/csh run scripts used in EPA Emission Modeling Platform (EMP). It reads a configuration JSON file to determine how to group and annotate variables, then rewrites the EMP input runscript with improved structure and optional path validation.

## Features
- Groups `setenv`/`set` variable assignments by user-defined sections.
- Preserves original order of non-matched lines.
- Annotates grouped variables under section headers.
- Optionally checks file paths and flags missing or duplicate assignments.

## Usage
```tcsh
python smkscript.py [--config PATH_TO_CONFIG] [--check-path] [--dry-run]
```
- `--config`: Path to the JSON config file (default: `smkscript_config.json` in the script directory).
- `--check-path`: Validate file paths referenced in variable assignments; alternative defined by `check_path` in configuration file
- `--dry-run`: Show planned grouping summary without writing output.

## Requirements
- Python 3.7 or newer
- No external dependencies (uses only Python standard library)

## Configuration
The script requires a JSON config file with the following keys:
- `runscript_in`: Path to the input tcsh/csh script
- `runscript_out`: Path to the output script
- `sections`: Dictionary mapping section names to lists of variable prefixes
- `directory_definitions`: Path to a definitions script for environment variables

## Example
See `smkscript_config.json` for a template configuration.

## Special Case
To validate file paths in output runscript after modification; change configuration file (e.g., `smkscript_config.json`) so that `runscript_in` and `runscript_out` both point to same the modified runscript, then rerun smkscript.py with `--check-path` (or `check_path`: true in configuration file)

## Author
Contact: [Huy Tran: tranhuy@email.unc.edu]
