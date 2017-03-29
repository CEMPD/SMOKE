#!/bin/bash

# location of PHOUR file created by Smkinven
export PHOUR=$INVOPD/phour.ptegu.2011ek.ncf

# location of AERMOD-source report created by Smkreport
export REPORT=$REPSCEN/pointgt.aermod.ptegu.12US2.2011ek.2011001.rpt

# locations of temporal profile files
export PTPRO_MONTHLY=$GE_DAT/noncem_2014_tpro_monthly_csv_19jul2016_v0
export PTPRO_DAILY=$GE_DAT/noncem_2014_tpro_daily_07_21_2016_19jul2016_v0
# used for Jan, Feb, Mar, Apr, Oct, Nov, Dec
export PTPRO_HOURLY_WINTER=$GE_DAT/ptegu_2014_tpro_hourly_winter_07_21_2016_csv_21jul2016_v0
# used for May through Sep
export PTPRO_HOURLY_SUMMER=$GE_DAT/ptegu_2014_tpro_hourly_summer_07_21_2016_csv_21jul2016_v0

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

# location of intermediate PHOUR_OUT text file
export PHOUR_OUT=$OUTPUT_DIR/phour.ptegu.2011ek.txt

export YEAR=2011

export PROMPTFLAG=N
$SCRIPTS/aermod/convert_phour

$SCRIPTS/aermod/ptegu.pl
