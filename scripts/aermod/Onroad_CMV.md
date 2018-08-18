# Onroad and CMV Sector Processing

## Requirements

Python 3.5+ with the packages:

* Pandas 0.20.1+
* Numpy 1.12.1+
* Pyproj 1.9.5.1+
* GDAL 2.2.0+

Functions may not be backwards compatible with python 2.7 or earlier versions of these python packages.

Compatibility can be tested by building the package under scripts/aermod/smk2ae:

`python setup.py build`

## Onroad

### Setup and Description

Onroad processing supports five general AERMOD run groups: HDON, LDON, HDOFF, LDOFF, and HOTEL. These run groups are defined by source category (SCC) in the SOURCE_GROUPS ancillary file. The SCCs associated with the specified run groups are selected from the onroad FF10.

Table 1.1. Format of the SOURCE_GROUPS comma delimited ancillary file.

Position|Name|Type|Description
-|-|-|-
1|source_group|char|AERMOD source group name
2|run_group|char|AERMOD run group name
3|nata_cat|char|NATA category
4|nei_cat|char|NEI category
5|scc|char(20)|Source Category Code

The grid resolutions are specified in the processing script described below. For 12 km CONUS domains, the run groups HDON, LDON, and HOTEL generate 4 km AERMOD helper files; the RUN_GROUPS environment variable is specified as HDON4, LDON4, and HOTEL4 respectively. The off-network groups, HDOFF and LDOFF, generate 12 km AERMOD helper files for 12 km CONUS domains. These run groups are specified as HDOFF12 and LDOFF12.

For non-CONUS state specific domains, the resolution is the same across all run groups. These run groups are specified by the name of the run group followed by the resolution and state-specific suffix: eg. HDON9AK, LDOFF3PR, HDOFF3HI, etc. All five run groups may not exist within smaller domains depending on which SCCs are available. When no data is available within a selected run group no helper files will be returned for that run group.

A monthly nonpoint flat file (FF10) ([format](https://www.cmascenter.org/smoke/documentation/4.5/html/ch08s02s04.html#d0e37919)) is used as the source of emissions data by county FIPS code, SCC, and pollutant. This processor was designed to support inventories derived from the output of an annual SMOKE-MOVES run. Multiple inventories may be specified using iterative letters as a suffix to the EMISINV variable.

Hourly county temporal profiles are developed for use with this processor based on hourly SMOKE-MOVES emissions reports by county and SCC. These reports are aggregated up to each run group for PM2.5 and benzene. Hourly scalars by county and run group are then calculated using PM2.5 for the HDON, HDOFF, and HOTEL run groups while benzene is used for the LDON and LDOFF run groups. The resulting scalars are then converted directly to SMOKE-AERMOD temporal helper files for county and run group using the supplemental onroad_temporal tool. SMOKE-MOVES emissions reports are not available for non-CONUS domains. National averages are used for the non-CONUS temporal profiles.

A list of environment variables is given in the tables below. All environment variables should be properly defined for successful completion.

The processing script is packaged with SMOKE4.5 and later under the scripts/aermod/smk2ae distribution path. The script library requirements are given at the top of this document. No additional modification of the scripts or libraries are needed.

### Inputs

Table 2.1. SMOKE-AERMOD onroad input file locations

Environment Variable|File Description
-|-
GROUP_PARAMS|Area release parameters (release height and sigma z) by run group
SOURCE_GROUPS|SCC to source and run group cross-reference
EMISINV_A|Area FF10 containing monthly emissions. May have multiple FF10 inputs by iterating the letter in the suffix (eg. EMISINV_B, EMISINV_C)
INVTABLE|SMOKE compatible inventory table specifying which pollutants to run. The NAME is the name of the pollutant that will be in the emis helper files and is also referred to as “SMOKE name”.
SRGPRO|Domain specific SMOKE-compatible gridding surrogates (4 km for HDON4, LDON4, and HOTEL4)
SRGDESC|Domain specific SMOKE-compatible surrogate description file
COSTCY|SMOKE-compatible country, state, county definition file
GRIDDESC|SMOKE and IOAPI compatible grid description file that contains all grids and projections to be run.
AGREF|SMOKE-compatible gridding surrogate cross-reference
ATREF|Temporal cross-reference. Not used for onroad runs, but environment variable definition is required to run.
ATPRO_HOURLY,<br>ATPRO_MONTHLY,<br>ATPRO_WEEKLY|Temporal profile files. Not used for onroad runs, but environment variable definition is required to run.

### Parameters

Table 2.2 SMOKE-AERMOD onroad parameters

Environment Variable|Parameter Description
-|-
RUN_GROUPS|Comma delimited list of run groups to process. If set to empty, defaults to process all run groups found in EMISINV. (ex. “HDON4,LDOFF12”)
STATE_FIPS|Comma delimited list of state FIPS to process. If set to empty, defaults to process all states found in EMISINV. (ex. “01,12,48”)
RUN_GROUP_SUFFIX|Output suffix for run group. This is either the cell size for CONUS domains (12 or 4) or the cell size and state for non-CONUS (9AK, 3HI, 3PR).
BASE_YEAR|Current run year as a 4 character string (eg. 2014)
WORK_PATH|Root tree for AERMOD helper output files.
REGION_IOAPI_GRIDNAME|Output grid for AERMOD helper files. Must subset in selected surrogates.

### Running

A driver script should be used to define the environment variables described in tables 2.1 and 2.2. Each of the specified input and parameter environment variables is required for successful processing. Once the environment variables are defined the script `aermod.py` is called from the driver script to initiate processing. This script is available in distributions of SMOKE 4.5 and later under: scripts/aermod/smk2ae/scripts

Due to the fine scale the 4 km CONUS run groups require more RAM than the 12 km or non-CONUS run groups. On most systems HDON4 and LDON4 in particular should be run individually rather than as a group run. Alternatively, a subset of states may be run to reduce memory usage using the STATE_FIPS environment variable.

### Output

Emissions by run group and resolution are output to the “emis” subdirectory under the WORK_PATH using a naming scheme: [GRID RESOLUTION]_[RUN_GROUP]_emis.csv (eg. 12_NPLO12_emis.csv)

Locations and parameters files are found under the respective “locations” and “parameters” subdirectories with a similar name scheme that omits the leading grid resolution.

Temporal helper files are not generated by this script and use the supplemental conversion to generate county hourly temporal files by run group. The supplemental conversion utility generates county hourly temporal profiles by run group and state.

A special county to grid cross-reference file is output under the “xwalk” subdirectory. This file is used to link grid cells to the county of greatest spatial overlap. This file is used in conjunction with the temporal helper files to associate a temporal profile with a grid cell.

## CMV

### Setup and Description

Commercial marine vessel processing produces a single run group, CMV, with multiple source groups as defined in the source groups file (table 1.1). Unlike other AERMOD run groups the CMV helper files are based on polygons rather than grid cells. A polygon facility ID based on the polygon shape ID and the county FIPS is used to cross-reference the FF10 inventory data to the polygon vertices.

The polygon file (table 3.1) is a comma delimited file that contains the vertices of all the CMV polygon vertices, both port and underway. These are simplified polygons based on NEI port and underway shapes. For ports, the simplified polygons are the same as the NEI port shapes.  For underway, simplified polygons that had been previously developed were used.  There may be multiple polygons per NEI underway shape. The emissions for each simplified polygon are assumed proportional to the simplified polygon area comprising underway shape.

Table 3.1. Format of the POLY_FILE comma delimited ancillary file.

Position|Name|Type|Description
-|-|-|-
1|facid|char|Facility ID including polygon number and county FIPS
2|lon|float|Longitude of polygon vertex
3|lat|float|Latitude of polygon vertex
4|numvert|int|Number of vertices in the polygon
5|area|float|Area of polygon

The polygon file facility ID is a concatenation of the five-character county FIPS and the polygon shape ID associated with the port or underway shape. The facility ID is generated from the region_cd and shape_id fields in the FF10 to cross reference county source group emissions to the polygon. For port sources the facility ID is the letter “P” concatenated with the shape_id, followed by the letter “F”, and ending with the region_cd. For underway sources the facility ID is the letter “U” concatenated with the zero padded four-character shape_id, followed by the letter “F”, the region_cd, the letter “P”, and ending with a polygon number.

An annual nonpoint flat file (FF10) is used to obtain emissions by county FIPS, shape ID, SCC, and pollutant. The shape_id field is required in order to form valid facility IDs to cross reference to the polygons as described above. Multiple inventories may be specified using iterative letters as a suffix to the EMISINV variable.

The processing script is packaged with SMOKE4.5 and later under the scripts/aermod/smk2ae distribution path. The script library requirements are given at the top of this document. No additional modification of the scripts or libraries are needed.

### Inputs

Table 4.1. SMOKE-AERMOD CMV input file locations

Environment Variable|File Description
-|-
POLY_FILE|Comma delimited file of lat/lon vertices and area of simplified port and underway shape
SOURCE_GROUPS|SCC to source and run group cross-reference
EMISINV_A|Area FF10 containing monthly emissions. May have multiple FF10 inputs by iterating the letter in the suffix (eg. EMISINV_B, EMISINV_C)
INVTABLE|SMOKE compatible inventory table with SMOKE names
COSTCY|SMOKE-compatible country, state, county definition file
GRIDDESC|SMOKE and IOAPI compatible grid description file that contains all grids and projections to be run.
ATREF|Temporal cross-reference by county and source code (SCC).
ATPRO_HOURLY,<br>ATPRO_MONTHLY,<br>ATPRO_WEEKLY|Temporal profile files for hourly (diurnal), monthly, and weekly profiles. Profile codes must match those assigned to the county and source codes specified in the cross-reference. For CMV only the monthly file is used, but all variables must be assigned.

### Parameters

Table 4.2 SMOKE-AERMOD CMV parameters

Environment Variable|Parameter Description
-|-
STATE_FIPS|Comma delimited list of state FIPS to process. If set to empty, defaults to process all states found in EMISINV. (ex. “01,12,48”)
WORK_PATH|Root tree for AERMOD helper output files.
REGION_IOAPI_GRIDNAME|Output grid for AERMOD helper files. Must subset in selected surrogates.

### Running

A driver script should be used to define the environment variables described in tables 4.1 and 4.2 above. Each of the specified input and parameter environment variables is required for successful processing. Once the environment variables are defined the script `cmv_smk2ae.py` is called in the driver script to initiate processing. This script is available in distributions of SMOKE 4.5 and later under: scripts/aermod/smk2ae/scripts

### Output

Emissions for all supplied facilities in the run group are found in the “emis” subdirectory under the WORK_PATH. The emissions file is typically named “CMV_emis.csv”

Locations, parameters, and temporal files are found under the respective “locations”, “parameters”, and temporal subdirectories. Each helper file follows a naming scheme that is typically CMV [helper type]. For example, CMV_parameters.csv
