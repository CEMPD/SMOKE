# smkextract

This utility facilitates filtering and processing emission inventory files (FF10 format) from existing EPA's Emission Modeling Platform (EMP) to match a target modeling grid, intersecting counties, or specific column values.

The toolkit comprises five component python scripts:
1. `build_sector_config.py`: Automates configuration setup by parsing EMP run scripts.
2. `smkxtgrid.py`: Filters inventory files based on a modeling grid (GRIDDESC) and county intersection.
3. `smkxtgis.py`: Filters inventory files based on intersection with a custom shapefile/domain.
4. `smkxtcol.py`: Filters inventory files based on specific column values or ranges.
5. `griddesc2shp.py`: Converts GRIDDESC files to shapefiles/GeoPackages for visualization.

# Setup Instructions

These scripts require Python 3 and several geospatial and data processing libraries. It is recommended to use a Conda environment.

## 1. Create and Activate Conda Environment

```sh
conda create -n smkextract python=3.9
conda activate smkextract
```

## 2. Install Dependencies

Install the required packages (`geopandas`, `shapely`, `pyproj`, `pandas`) using `conda-forge`:

```sh
conda install -c conda-forge geopandas shapely pyproj pandas
```

*Note: Standard Python libraries used include `argparse`, `json`, `os`, `re`, `sys`, `csv`, `warnings`.*

---

# A. build_sector_config.py

## Overview

`build_sector_config.py` is a utility script designed to automate the population of sector entries in a configuration JSON file for `smkxtgrid.py` (by default `smkxtgrid.json`) or other extraction scripts. It parses a set of EMP run scripts to extract sector names and emission inventory file references, then updates the configuration file accordingly. This helps streamline the setup and management of emission sector data for large-scale air quality modeling projects.

`build_sector_config.py` requires a manifest of EMP runscripts in form of JSON (e.g., `emp_runscripts.json`) and a target JSON configuration input file that it will make update to (e.g., `smkxtgrid.json`). In a nutshell, `build_sector_config.py` read each EMP runscript to extract the `SECTOR` variable and all emission inventory file paths defined by variables like `EMISINV_*`, `EMISDAY_*`, or `EMISHOUR_*`, and then update the target configuration with list of found emission inventory files for each found sector accordingly.

## Usage

```sh
python build_sector_config.py [--runscripts <path/to/emp_runscripts.json>] [--config <path/to/smkxtgrid.json>]
```

- `--runscripts`: Path to the manifest JSON file listing run scripts (default: `emp_runscripts.json` in the script directory).
- `--config`: Path to the configuration JSON file to update (default: `smkxtgrid.json` in the script directory).

## Example

To update `smkxtgrid.json` using the default manifest:

```sh
python build_sector_config.py
```

To specify custom paths:

```sh
python build_sector_config.py --runscripts /path/to/emp_runscripts.json --config /path/to/smkxtgrid.json
```

# B. smkxtgrid.py

## Overview

`smkxtgrid.py` is the main script for filtering and processing emission inventory files (FF10 format) to match a target modeling grid and intersecting counties using a GRIDDESC file.

## Features

- **Grid and County Intersection:** Reads a grid description file (GRIDDESC) and a county shapefile to determine which counties intersect the modeling domain.
- **FF10 Filtering:** Filters FF10 emission inventory files to include only records for counties intersecting the grid.
- **Flexible Configuration:** Uses a JSON configuration file (`smkxtgrid.json`) to specify input/output directories, grid, shapefile, and sector-to-file mappings.
- **Batch Processing:** Processes all sectors and files listed in the configuration, with support for skipping specified sectors.
- **Output:** Writes filtered FF10 files to the specified output directory, appending the grid ID to filenames.

## Usage

```sh
python smkxtgrid.py --config <path/to/smkxtgrid.json>
```

- `--config`: Path to the configuration JSON file (default: `smkxtgrid.json` in the script directory).

## Example

To filter all inventories for a grid and county set defined in `smkxtgrid.json`:

```sh
python smkxtgrid.py --config smkxtgrid.json
```

# C. smkxtgis.py

## Overview

`smkxtgis.py` filters FF10 emission inventory files based on their geographical intersection with a user-provided input shapefile (domain), rather than a GRIDDESC file. This is useful for analyzing emissions within a specific custom region (e.g., a basin, a state subset, or a custom polygon).

## Features

- **Spatial Intersection:** Determines intersecting counties based on a custom input shapefile (e.g., `PermianBasin_Extent.shp`) and a county boundary shapefile.
- **FF10 Filtering:** Filters inventory records to retain only those valid for the intersecting counties.
- **Configurable Inputs:** Uses a JSON configuration file to manage multiple sectors and files.

## Usage

```sh
python smkxtgis.py --config <path/to/smkxtgis.json>
```

- `--config`: Path to the configuration JSON file (default: `smkxtgis.json`).

## Configuration (`smkxtgis.json`)

The configuration file requires the following keys:
- `inputs`: Base directory for input files.
- `outputs`: Directory where filtered files will be saved.
- `county_shp`: Path to the US county shapefile (e.g., `.gpkg` or `.shp`).
- `input_shp`: Path to the domain shapefile defining the area of interest.
- `sector`: Dictionary mapping sector names to lists of filenames.
- `sector_skip`: (Optional) List of sectors to skip.

## Example

```sh
python smkxtgis.py --config smkxtgis.json
```

# D. smkxtcol.py

## Overview

`smkxtcol.py` is a general-purpose utility for filtering FF10 emission files based on specific column values or value ranges. This is useful for extracting data for specific facility IDs, SCCs, or regions without needing a spatial file.

## Features

- **Column-Based Filtering:** Filter rows where a specified column matches a list of values (`filtered_val`) or falls within a range (`start_val` to `end_val`).
- **Structure Preservation:** Retains original file headers and comments.
- **Flexible Data Handling:** Supports various FF10 formats by dynamically detecting headers and processing rows.

## Usage

```sh
python smkxtcol.py --config <path/to/smkxtcol.json>
```

- `--config`: Path to the configuration JSON file (default: `smkxtcol.json`).

## Configuration (`smkxtcol.json`)

Required keys:
- `inputs`: Directory for input files.
- `outputs`: Directory for output files.
- `sector` or `files`: List of files to process.
- `filter_col`: The name of the column to apply the filter on (e.g., `region_cd`).
- `filtered_val`: (Optional) A list of exact values to keep.
- `start_val` / `end_val`: (Optional) A range of values to keep.

## Example

To filter files for specific NY counties (FIPS starting with 36...):

```sh
python smkxtcol.py --config smkxtcol.json
```

# E. griddesc2shp.py

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

## Example

To create a full grid GeoPackage for grid `12LISTOS`:

```sh
python griddesc2shp.py -g griddesc.txt -i 12LISTOS
```

To create a single extent polygon:

```sh
python griddesc2shp.py -g griddesc.txt -i 12LISTOS --extent-only
```

# F. Author
Huy Tran: tranhuy@email.unc.edu
