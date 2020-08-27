# Introduction

These scripts are designed to create "helper" files that can then be transformed into AERMOD inputs. For most sectors, SMOKE is used to process the emissions inventory and create report files. These reports are then post-processed to produce the helper files.

The scripts are a mix of Perl (version 5.8+), Python (version 3.5+), and shell scripts. The Python scripts require the following packages to be installed:

* Pandas 0.20.1+
* Numpy 1.12.1+
* Pyproj 1.9.5.1+
* GDAL 2.2.0+

The Perl dependencies are included in the Git repository:

* Date::Simple
* Geo::Coordinates::UTM
* Text::CSV

# Design Documents

* [CMV_PORT_UNDERWAY_AERMOD_INPUT_FILES.pdf](https://github.com/CEMPD/SMOKE/blob/master/scripts/aermod/docs/CMV_PORT_UNDERWAY_AERMOD_INPUT_FILES.pdf)
* [NONPOINT_ONROAD_AERMOD_INPUT_FILE.pdf](https://github.com/CEMPD/SMOKE/blob/master/scripts/aermod/docs/NONPOINT_ONROAD_AERMOD_INPUT_FILE.pdf)
* [POINT_AIRPORT_AERMOD_INPUT_FILE.pdf](https://github.com/CEMPD/SMOKE/blob/master/scripts/aermod/docs/POINT_AIRPORT_AERMOD_INPUT_FILE.pdf)

# Input Data Files

SMOKE flatfiles input for SMOKE4AERMOD for the 2014 NATA application and other inputs for the 2014 NATA modeling platform can be downloaded from EPA's modeling platform website: https://www.epa.gov/air-emissions-modeling/2014-version-71-platform  

# Sectors

## Point sectors
* [ptnonipm](#ptnonipm-sector)
* [ptegu](#ptegu-sector)

## Airports
* [ptairport](#ptairport-sector)

## Area sectors
* [nonpoint](#nonpoint-sector)
* [nonroad](#nonroad-sector)
* [np_oilgas](#np_oilgas-sector)
* [rwc and ag](#rwcag-sectors)
* [onroad](#onroad-sector)

## Commercial marine vessel sources
* [cmv](#cmv-sector)

## ptnonipm sector

The SMOKE programs Smkinven, Grdmat, and one day of Temporal (to assign temporal profiles) need to be run for the ptnonipm inventory. When running Smkinven, the INVTABLE file should be used to restrict the inventory to sources with HAP or diesel PM emissions.

### Smkreport configuration

To group the inventory sources into AERMOD sources, run Smkreport with the following REPCONFIG file:

    SMK_SOURCE P

    /NEWFILE/ REPORT1

    /CREATE REPORT/
        AERMOD POINT PTNONIPM
    /END/
    
    /NEWFILE/ REPORT2
    
    /CREATE REPORT/
        BY SOURCE
    /END/

This creates a report where the inventory sources are grouped by:

* facility
* facility source type
* temporal profile (month-of-year, day-of-week, and hour-of-day)
* emission release point type
* stack parameters
* fugitive parameters
* latitude and longitude

The report also includes the grid cell, coordinates in Lambert and UTM projections, and the facility name. Smkreport also creates a crosswalk file, listing inventory sources that contribute to each record in the report output.

A second report gives emissions for each inventory source.

### Post-processing

The Perl script ptnonipm.pl reads the Smkreport files and temporal profile files to create the AERMOD helper outputs:

    locations/ptnonipm_location.csv
    parameters/ptnonipm_fug_srcparam.csv
    parameters/ptnonipm_point_srcparam.csv
    temporal/ptnonipm_temporal.csv
    xwalk/ptnonipm_process_releasept_emis.csv
    xwalk/ptnonipm_srcid_xwalk.csv

The script uses environment variables to locate the input files ($REPORT, $REP\_XWALK, $REP\_SRC, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY) and the directory where the outputs will be written ($OUTPUT\_DIR). The shell script run_ptnonipm.sh sets up the environment variables and runs ptnonipm.pl.

After processing both the ptnonipm and ptegu sectors, the outputs need to be combined use the Perl script combine_point.pl. This script uses the environment variable $OUTPUT_DIR to locate the individual sector output files, and write the combined files:

    locations/point_location.csv
    parameters/point_point_srcparam.csv
    parameters/point_fug_srcparam.csv
    xwalk/point_combined_process_releasept_emis.csv
    xwalk/point_srcid_xwalk.csv

The shell script run_combine_point.sh sets up the $OUTPUT_DIR environment variable and runs combine_point.pl.

## ptegu sector

When running Smkinven, the environment variable FLOW\_RATE\_FACTOR should not be set. The environment variable OUTPUT\_LOCAL\_TIME should be set to Y.

The utility program convert\_phour is used to read the NetCDF PHOUR file created by Smkinven, and output a text file containing the hourly factors for the CEM sources. To compile convert\_phour, run `make`.

The PHOUR file output from Smkinven needs to contain all hours of the year, at least 365*24 timesteps. convert\_phour determines the first time step to read from PHOUR using the YEAR environment variable. The first time step will be January 1st of the specified year, hour 0.

### Smkreport configuration

    SMK_SOURCE P

    /NEWFILE/ REPORT1

    /CREATE REPORT/
        AERMOD POINT PTEGU
    /END/
    
    /NEWFILE/ REPORT2
    
    /CREATE REPORT/
        BY SOURCE
    /END/

This uses the same grouping as the ptnonipm sector. Sources that use hourly CEM data will have their temporal profile ID reported as 'HR' + source ID.

### Post-processing

The Perl script ptegu.pl reads the Smkreport files, temporal profile files, and text version of the PHOUR file to create the AERMOD helper outputs:

    locations/ptegu_location.csv
    parameters/ptegu_fug_srcparam.csv
    parameters/ptegu_point_srcparam.csv
    temporal/<facility_id>_<state_code>_hourly.csv
    xwalk/ptegu_process_releasept_emis.csv
    xwalk/ptegu_srcid_xwalk.csv

The shell script run_ptegu.sh runs both convert\_phour and ptegu.pl. convert\_phour uses the environment variables $PHOUR, $PHOUR\_OUT, and $YEAR. ptegu.pl uses $REPORT, $REP\_XWALK, $REP\_SRC, $PHOUR\_OUT, $YEAR, $PTPRO\_MONTHLY, $PTPRO\_DAILY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY\_WINTER, $PTPRO\_HOURLY\_SUMMER, and $OUTPUT\_DIR.

After running the ptegu sector, the outputs need to be combined with the outputs from the ptnonipm sector. See the [ptnonipm](#ptnonipm-sector) section for more information.

## Point sectors quality assurance

The Python script point_qa.py generates quality assurance outputs based on the post-processed AERMOD helper files. The QA outputs contain comparisons of AERMOD sources between the helper files and comparisons against the inventory files:

    point_counts.csv
    point_emis_qa.csv
    point_fug_param_qa.csv
    point_point_param_qa.csv
    point_poll_qa.csv
    point_srcid_qa.csv
    temporal_point_check.csv

The Python script uses environment variables $CASE, $PROJECT_ROOT, $WORK_PATH, $INVTABLE, $SECTOR, and $EMISINV_[A-Z] to locate the AERMOD helper and inventory files. The script may be run for each sector, but it is recommended to run against the combined point file with all point inventories included.

## ptairport sector

Airports should be separated from the ptnonipm inventory and processed as a separate sector through SMOKE. Like ptnonipm, Smkinven, Grdmat, and one day of Temporal need to be run.

### Smkreport configuration

    SMK_SOURCE P

    /NEWFILE/ REPORT1

    /CREATE REPORT/
        AERMOD POINT
    /END/

This report groups inventory sources by:

* facility
* facility source type
* temporal profile (month-of-year, day-of-week, and hour-of-day)

The report also includes the grid cell, coordinates in lat-lon, Lambert, and UTM projections, and the facility name.

### Post-processing

The Perl script ptairport.pl reads the Smkreport file, temporal profile files, and the runway data file to create the AERMOD helper outputs:

    locations/airport_line_locations.csv
    locations/airport_nonrunway_locations.csv
    parameters/airport_line_params.csv
    parameters/airport_nonrunway_params.csv
    temporal/airport_line_temporal.csv
    temporal/airport_nonrunway_temporal.csv
    xwalk/airport_srcid_emis.csv

The shell script run_ptairport.sh sets up the environment variables needed by ptairport.pl ($REPORT, $RUNWAYS, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY, and $OUTPUT\_DIR).

### Runway data file

The runway data file lists start and end coordinates and widths for runways corresponding to each airport facility. A default file [runway_source_airports_endpoints_2017_02dec2019_v0.csv](https://raw.githubusercontent.com/CEMPD/SMOKE/master/scripts/aermod/runway_source_airports_endpoints_2017_02dec2019_v0.csv) is provided.

### Quality assurance

The Python script airport_qa.py generates quality assurance outputs based on the post-processed AERMOD helper files. The QA outputs contain comparisons of AERMOD sources between the helper files and comparisons against the inventory files:

    airport_counts.csv
    airport_emis_qa.csv
    airport_locs_check.csv
    temporal_airport_check.csv
    airport_poll_qa.csv
    airport_srcid_qa.csv

The Python script uses environment variables $CASE, $PROJECT_ROOT, $WORK_PATH, $INVTABLE, $SECTOR, and $EMISINV_[A-Z] to locate the AERMOD helper and inventory files. A second Python script nonrunway_loc.py generates an airport_locs_check.csv for nonrunway airports. 

## nonpoint sector

For the nonpoint sector, Smkinven, Grdmat, and one day of Temporal (to assign temporal profile) need to be run before creating the custom AERMOD report from Smkreport.

### Smkreport configuration

    SMK_SOURCE A

    /NEWFILE/ REPORT1

    /CREATE REPORT/
        AERMOD NONPOINT
    /END/

This creates a report where the inventory sources are grouped by:

* grid cell
* county
* SCC
* temporal profile (month-of-year, day-of-week, and hour-of-day)

The report also includes the latitude and longitude of each grid cell corner.

### Post-processing

The Perl script nonpt.pl reads the Smkreport file, temporal profile files, and run group configuration files to create the AERMOD helper outputs:

    nonpt_locations.csv
    nonpt_area_params.csv
    nonpt_temporal.csv
    nonpt_emis.csv

The shell script run_nonpt.sh sets up the environment variables needed by nonpt.pl ($REPORT, $SOURCE_GROUPS, $GROUP_PARAMS, $ATPRO_MONTHLY, $ATPRO_WEEKLY, $ATPRO_HOURLY, and $OUTPUT_DIR).

### Run group configuration files

The source groups file assigns each SCC to a source group and a run group. A default file [hem_aermod_groups_2017_11jun2020_nf_v2.csv](https://github.com/CEMPD/SMOKE/blob/master/scripts/aermod/hem_aermod_groups_2017_11jun2020_nf_v2.csv) is provided.

The group parameters file provides the release height and initial vertical dispersion (sigma z) values for each run group. Both values are in meters. A default file [nonpoint_rungroup_stack_parameters_11jun2020_v4.csv](https://github.com/CEMPD/SMOKE/blob/master/scripts/aermod/nonpoint_rungroup_stack_parameters_11jun2020_v4.csv) is provided.

## nonroad sector

After running Smkinven, Grdmat, and Temporal, Smkreport should be run for each month of the year. This will produce a set of 12 custom AERMOD reports with emissions for each month.

### Smkreport configuration

Same as the [nonpoint sector](#nonpoint-sector) Smkreport configuration.

### Post-processing

The Perl script nonroad.pl reads the monthly Smkreport files, temporal profile files, and run group configuration files to create the AERMOD helper outputs:

    nonroad_locations.csv
    nonroad_area_params.csv
    nonroad_temporal.csv
    nonroad_emis.csv

The shell script run_nonroad.sh sets up the environment variables needed by nonroad.pl ($REPORT_JAN, $REPORT_FEB, ..., $REPORT_DEC, $SOURCE_GROUPS, $GROUP_PARAMS, $ATPRO_MONTHLY, $ATPRO_WEEKLY, $ATPRO_HOURLY, and $OUTPUT_DIR).

### Run group configuration files

Same as the [nonpoint sector](#nonpoint-sector) configuration files.

## np_oilgas sector

The np_oilgas sector post-processing scripts assume that a 4 km grid was used in the SMOKE processing, and that the 4 km cells should be reported within a 12 km grid.

### Smkreport configuration

Same as the [nonpoint sector](#nonpoint-sector) Smkreport configuration.

### Post-processing

The Perl script np_oilgas.pl reads the Smkreport file, temporal profile files, and run group configuration files to create the AERMOD helper outputs:

    np_oilgas_locations.csv
    np_oilgas_area_params.csv
    np_oilgas_temporal.csv
    np_oilgas_emis.csv

The shell script run_np_oilgas.sh sets up the environment variables needed by np_oilgas.pl ($REPORT, $SOURCE_GROUPS, $GROUP_PARAMS, $ATPRO_MONTHLY, $ATPRO_WEEKLY, $ATPRO_HOURLY, and $OUTPUT_DIR).

### Run group configuration files

Same as the [nonpoint sector](#nonpoint-sector) configuration files.

## rwc/ag sectors

The rwc and ag sectors read day-of-month temporal profiles in addition to month-of-year, day-of-week, and hour-of-day profiles. Currently, the rwc/ag post-processing requires that each individual source use the same hour-of-day profile for every day of the week. The rwc/ag post-processing is set up for the year 2014.

### Smkreport configuration

Same as the [nonpoint sector](#nonpoint-sector) Smkreport configuration.

### Post-processing

The Perl script rwc.pl reads the Smkreport file, temporal profile files, and run group configuration files to create the AERMOD helper outputs:

    rwc_locations.csv
    rwc_area_params.csv
    rwc_temporal.csv
    rwc_county_to_gridcell.csv
    rwc_emis.csv

The shell script run_rwc.sh sets up the environment variables needed by rwc.pl ($REPORT, $SOURCE_GROUPS, $GROUP_PARAMS, $ATPRO_MONTHLY, $ATPRO_DAILY, $ATPRO_WEEKLY, $ATPRO_HOURLY, and $OUTPUT_DIR).

### Run group configuration files

Same as the [nonpoint sector](#nonpoint-sector) configuration files.

## onroad sector

The onroad sector uses annual emissions reports generated by SMOKE-MOVES. The Python script aermod.py reads the SMOKE-MOVES reports, grid data, and run group configuration files to create the AERMOD helper inputs. More details can be found in the [onroad and CMV sector documentation](Onroad_CMV.md).

## Area sectors quality assurance

The Python script area_qa.py generates quality assurance outputs based on the post-processed AERMOD helper files. The QA outputs contain comparisons of AERMOD sources between the helper files and comparisons against the inventory files. A set of QA files are generated for each run group:

    rungroup_counts.csv
    rungroup_emis_qa.csv
    rungroup_grid_comparison.csv
    rungroup_poll_qa.csv
    rungroup_scc_qa.csv
    rungroup_srcid_qa.csv
    rungroup_temporal_check.csv

The Python script uses environment variables $REGION_IOAPI_GRIDNAME, $WORK_PATH, $INVTABLE, $GROUP_PARAMS, $SOURCE_PARAMS, and $EMISINV_[A-Z] to locate the AERMOD helper and inventory files. The script is run for a set of inventories that cover all run groups in the platform sector.

## cmv sector

The CMV sector is not processed through SMOKE since the helper files are based on polygons rather than grid cells. The Python script cmv_smk2ae.py reads a polygon file, an annual nonpoint flat file (FF10), temporal profile files, and run group configuration files to create the AERMOD helper outputs. More details can be found in the [onroad and CMV sector documentation](Onroad_CMV.md).

### Quality assurance

The Python script cmv_qa.py generates quality assurance outputs based on the post-processed AERMOD helper files. The QA outputs contain comparisons of AERMOD sources between the helper files and comparisons against the inventory files:

    CMV_counts.csv
    CMV_emis_qa.csv
    CMV_poll_qa.csv
    CMV_scc_qa.csv
    CMV_srcid_qa.csv
    CMV_temporal_check.csv

The Python script uses environment variables $WORK_PATH, $INVTABLE, and $EMISINV_[A-Z] to locate the AERMOD helper and inventory files. The script is run for a set for all CMV inventories.
