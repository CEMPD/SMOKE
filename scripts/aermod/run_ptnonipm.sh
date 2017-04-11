#!/bin/bash

# location of AERMOD-source report created by Smkreport
export REPORT=$REPSCEN/pointgt.aermod.ptnonipm.12US2.2011ek.2011001.rpt
export REP_XWALK=$REPSCEN/pointgt.aermod.ptnonipm.12US2.2011ek.2011001.rpt_src_crosswalk.txt

# locations of temporal profile files
export PTPRO_MONTHLY=$GE_DAT/amptpro_general_2011platform_tpro_monthly_20nov2015_v2
export PTPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
export PTPRO_HOURLY=$GE_DAT/amptpro_general_2011platform_tpro_hourly_13nov2014_v1

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

$SCRIPTS/aermod/ptnonipm.pl
