# smkextract

This utility facilitates filtering and processing emission inventory files (FF10 format) from existing EPA's Emission Modeling Platform (EMP) to match a target modeling grid and intersecting counties. The target modeling grid is required to be a subset of modeling grids that the EMP is developed for (e.g., 12US1, 12US2, 36US3, etc.)

This utilities three component python scripts: `smkxtgrid.py` which is the main script for filtering inventory files, and auxilarly `build_sector_config.py` script that facilitate inputs for `smkxtgrid.py`; `griddesc2shp.py` is an dependent script for generating shapefile from grid parameters from GRIDDESC file

# A. build_sector_config.py

## Overview

`build_sector_config.py` is a utility script designed to automate the population of sector entries in a configuration JSON file for `smkxtgrid.py`(by default `smkxtgrid.json`) for emission inventory files filtering. It parses a set of EMP run scripts to extract sector names and emission inventory file references, then updates the configuration file accordingly. This helps streamline the setup and management of emission sector data for large-scale air quality modeling projects.

`build_sector_config.py` requires a manifest of EMP runscripts in form of JSON (e.g., `emp_runscripts.json`) and a target JSON configuration input file of `smkxtgrid.py` that it will make update to (e.g., `smkxtgrid.json`). In a nutshell, `build_sector_config.py` read each EMP runscript to extract the `SECTOR` variable and all emission inventory file paths defined by variables like `EMISINV_*`, `EMISDAY_*`, or `EMISHOUR_*`, and then update `smkxtgrid.json` with list of found emission inventory files for each found sector accordingly.

## Usage

```sh
python build_sector_config.py [--runscripts <path/to/emp_runscripts.json>] [--config <path/to/smkxtgrid.json>]
```

- `--runscripts`: Path to the manifest JSON file listing run scripts (default: `emp_runscripts.json` in the script directory).
- `--config`: Path to the configuration JSON file to update (default: `smkxtgrid.json` in the script directory).

## Requirements

- **Python 3.x**
- Standard Python libraries: `argparse`, `json`, `os`, `re`, `sys`
- Access to the run scripts and configuration files referenced in the arguments or defaults

## Example

To update `smkxtgrid.json` using the default manifest:

```sh
python build_sector_config.py
```

To specify custom paths:

```sh
python build_sector_config.py --runscripts /path/to/emp_runscripts.json --config /path/to/smkxtgrid.json
```

## Notes

- The script expects run scripts to use `setenv` statements for variable assignment.
- Only unique filenames from emission inventory variables are added to each sector.
- The configuration file must contain a top-level `sector` object (dictionary).

# B. smkxtgrid.py

## Overview

`smkxtgrid.py` is the main script for filtering and processing emission inventory files (FF10 format) to match a target modeling grid and intersecting counties.

## Features

- **Grid and County Intersection:** Reads a grid description file and a county shapefile to determine which counties intersect the modeling domain.
- **FF10 Filtering:** Filters FF10 emission inventory files to include only records for counties intersecting the grid.
- **Flexible Configuration:** Uses a JSON configuration file (`smkxtgrid.json`) to specify input/output directories, grid, shapefile, and sector-to-file mappings.
- **Batch Processing:** Processes all sectors and files listed in the configuration, with support for skipping specified sectors.
- **Output:** Writes filtered FF10 files to the specified output directory, appending the grid ID to filenames.

## Usage

```sh
python smkxtgrid.py --config <path/to/smkxtgrid.json>
```

- `--config`: Path to the configuration JSON file (default: `smkxtgrid.json` in the script directory).

## Requirements

- **Python 3.x**
- Python packages: `geopandas`, `shapely`, `pyproj`, `pandas`
- Standard Python libraries: `argparse`, `json`, `os`, `csv`, `warnings`
- Access to the required grid description file, county shapefile, and emission inventory files as specified in the configuration

## Example

To filter all inventories for a grid and county set defined in `smkxtgrid.json`:

```sh
python smkxtgrid.py --config smkxtgrid.json
```

## Notes

- The script expects FF10 files to be in CSV format and to contain a region/county code column.
- The configuration file must specify all required paths and sector mappings.
- Output files are written to the output directory, with the grid ID appended to filenames.

# C. griddesc2shp.py

## Overview

`griddesc2shp.py` is a utility script for converting a SMOKE/CMAQ GRIDDESC file into a GeoPackage (.gpkg) polygon layer, either as a full grid (one polygon per cell) or as a single extent polygon. It is useful for visualizing or spatially analyzing the modeling domain and for generating shapefiles compatible with GIS tools.

## Features

- **GRIDDESC Parsing:** Reads SMOKE/CMAQ GRIDDESC files, supporting both coordinate and grid definition sections.
- **Flexible Output:** Generates either a full grid (with row/column attributes) or a single extent polygon.
- **Projection Handling:** Converts cell coordinates from native Lambert Conformal Conic (LCC) projection to WGS84 lon/lat.
- **Customizable Output:** Allows user to specify output file, layer name, and whether to include cell center attributes.
- **Command-Line Interface:** Supports a variety of CLI options for input, output, and processing mode.

## Usage

```
python griddesc2shp.py -g <path/to/GRIDDESC.txt> -i <GRID_ID> [options]
```

**Key options:**
- `-g/--griddesc`: Path to the GRIDDESC file (required)
- `-i/--grid-id`: Grid/domain name/ID to extract (required)
- `-o/--output`: Output GeoPackage file path (optional)
- `--layer-name`: Layer name inside GeoPackage (optional)
- `--extent-only`: Output only a single extent polygon (flag)
- `--no-centers`: Skip center coordinate attributes (flag)

## Requirements

- **Python 3.8+**
- Python packages: `geopandas`, `shapely`, `pyproj`
- Standard Python libraries: `re`, `os`, `argparse`
- Access to a valid SMOKE/CMAQ GRIDDESC file

## Example

To create a full grid GeoPackage for grid `12LISTOS`:

```sh
python griddesc2shp.py -g griddesc.txt -i 12LISTOS
```

To create a single extent polygon:

```sh
python griddesc2shp.py -g griddesc.txt -i 12LISTOS --extent-only
```

## Notes

- The script expects the GRIDDESC file to contain a coordinate section terminated by `' '  !  end coords.`
- By default, uses a spherical Earth for LCC projection (can be toggled in the script).
- Output is a GeoPackage (.gpkg) file, suitable for use in GIS software.



# D. Author
Huy Tran: tranhuy@email.unc.edu