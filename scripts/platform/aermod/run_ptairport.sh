#!/bin/bash

# location of AERMOD-source report created by Smkreport
export REPORT=$REPSCEN/pointgt.aermod.ptairport.12US2.2011ek.2011001.rpt

# location of runways data file
export RUNWAYS=$SCRIPTS/aermod/runway_source_airports_endpoints_2017_02dec2019_v0.csv

# locations of temporal profile files
export PTPRO_MONTHLY=$GE_DAT/amptpro_general_2011platform_tpro_monthly_20nov2015_v2
export PTPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
export PTPRO_HOURLY=$GE_DAT/amptpro_general_2011platform_tpro_hourly_13nov2014_v1

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

$SCRIPTS/aermod/ptairport.pl
