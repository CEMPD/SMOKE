# Sectors

* [ptnonipm](#ptnonipm-sector)
* [ptairport](#ptairport-sector)
* [ptegu](#ptegu-sector)
* [nonpoint](#nonpoint-sector)
* [nonroad](#nonroad-sector)
* [np_oilgas](#np_oilgas-sector)
* [rwc](#rwc-sector)

## ptnonipm sector

The SMOKE programs Smkinven, Grdmat, and one day of Temporal (to assign temporal profiles) need to be run for the ptnonipm inventory. When running Smkinven, the INVTABLE file should be used to restrict the inventory to sources with HAP or diesel PM emissions.

### Smkreport configuration

To group the inventory sources into AERMOD sources, run Smkreport with the following REPCONFIG file:

    SMK_SOURCE P

    /NEWFILE/ REPORT1

    /CREATE REPORT/
        AERMOD POINT PTNONIPM
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

### Post-processing

The Perl script ptnonipm.pl reads the Smkreport files and temporal profile files to create the AERMOD helper outputs:

    locations/ptnonipm_location.csv
    parameters/ptnonipm_fug_srcparam.csv
    parameters/ptnonipm_point_srcparam.csv
    temporal/ptnonipm_temporal.csv
    xwalk/ptnonipm_srcid_emis.csv
    xwalk/ptnonipm_srcid_xwalk.csv

The script uses environment variables to locate the input files ($REPORT, $REP\_XWALK, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY) and the directory where the outputs will be written ($OUTPUT\_DIR). The shell script run_ptnonipm.sh sets up the environment variables and runs ptnonipm.pl.

After processing both the ptnonipm and ptegu sectors, the outputs need to be combined use the Perl script combine_point.pl. This script uses the environment variable $OUTPUT_DIR to locate the individual sector output files, and write the combined files:

    locations/point_location.csv
    parameters/point_point_srcparam.csv
    parameters/point_fug_srcparam.csv
    xwalk/point_srcid_emis.csv
    xwalk/point_srcid_xwalk.csv

The shell script run_combine_point.sh sets up the $OUTPUT_DIR environment variable and runs combine_point.pl.

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

The runway data file lists start and end coordinates and widths for runways corresponding to each airport facility. A default file [Final_2014_Runway_Source_Airports_Endpoints_4.csv](https://raw.githubusercontent.com/CEMPD/SMOKE/master/scripts/aermod/Final_2014_Runway_Source_Airports_Endpoints_4.csv) is provided.

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

This uses the same grouping as the ptnonipm sector. Sources that use hourly CEM data will have their temporal profile ID reported as 'HR' + source ID.

### Post-processing

The Perl script ptegu.pl reads the Smkreport file, temporal profile files, and text version of the PHOUR file to create the AERMOD helper outputs:

    locations/ptegu_location.csv
    parameters/ptegu_fug_srcparam.csv
    parameters/ptegu_point_srcparam.csv
    temporal/<facility_id>_<state_code>_hourly.csv
    xwalk/ptegu_srcid_emis.csv
    xwalk/ptegu_srcid_xwalk.csv

The shell script run_ptegu.sh runs both convert\_phour and ptegu.pl. convert\_phour uses the environment variables $PHOUR, $PHOUR\_OUT, and $YEAR. ptegu.pl uses $REPORT, $REP\_XWALK, $PHOUR\_OUT, $YEAR, $PTPRO\_MONTHLY, $PTPRO\_DAILY, $PTPRO\_HOURLY\_WINTER, $PTPRO\_HOURLY\_SUMMER, and $OUTPUT\_DIR.

After running the ptegu sector, the outputs need to be combined with the outputs from the ptnonipm sector. See the [ptnonipm](#ptnonipm-sector) section for more information.

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

## rwc sector

The rwc sector reads day-of-month temporal profiles in addition to month-of-year, day-of-week, and hour-of-day profiles. Currently, the RWC post-processing requires that each individual source use the same hour-of-day profile for every day of the week. The RWC post-processing is set up for the year 2014.

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
