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
* temporal profile (month-of-year, day-of-week, and hour-of-day)
* emission release point type
* stack parameters
* fugitive parameters
* latitude and longitude

The report also includes the grid cell, coordinates in Lambert and UTM projections, and the facility name.

### Post-processing

The Perl script ptnonipm.pl reads the Smkreport file and temporal profile files to create the AERMOD helper outputs:

    point_location.csv
    point_point_srcparam.csv
    point_fug_srcparam.csv
    point_temporal.csv
    point_srcid_emis.csv

The script uses environment variables to locate the input files ($REPORT, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY) and the directory where the outputs will be written ($OUTPUT\_DIR). The shell script run_ptnonipm.sh sets up the environment variables and runs ptnonipm.pl.

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
* temporal profile (month-of-year, day-of-week, and hour-of-day)

The report also includes the grid cell, coordinates in lat-lon, Lambert, and UTM projections, and the facility name.

### Post-processing

The Perl script ptairport.pl reads the Smkreport file, temporal profile files, and the runway data file to create the AERMOD helper outputs:

    airport_line_locations.csv
    airport_line_params.csv
    airport_line_temporal.csv
    airport_nonrunway_locations.csv
    airport_nonrunway_params.csv
    airport_nonrunway_temporal.csv
    airport_srcid_emis.csv

The shell script run_ptairport.sh sets up the environment variables needed by ptairport.pl ($REPORT, $RUNWAYS, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY, and $OUTPUT\_DIR).

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

    point_location.csv
    point_point_srcparam.csv
    point_fug_srcparam.csv
    point_temporal.csv
    <facility_id>_hourly.csv
    point_srcid_emis.csv

The shell script run_ptegu.sh runs both convert\_phour and ptegu.pl. convert\_phour uses the environment variables $PHOUR, $PHOUR\_OUT, and $YEAR. ptegu.pl uses $REPORT, $PHOUR\_OUT, $PTPRO\_MONTHLY, $PTPRO\_WEEKLY, $PTPRO\_HOURLY, and $OUTPUT\_DIR.
