#!/bin/bash

# location of PHOUR file created by Smkinven
export PHOUR=$INVOPD/phour.ptegu.2011ek.ncf

# location of AERMOD-source report created by Smkreport
export REPORT=$REPSCEN/pointgt.aermod.ptegu.12US2.2011ek.2011001.rpt

# locations of temporal profile files
export PTPRO_MONTHLY=$GE_DAT/amptpro_general_2011platform_tpro_monthly_20nov2015_v2
export PTPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
export PTPRO_HOURLY=$GE_DAT/amptpro_general_2011platform_tpro_hourly_13nov2014_v1

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod/ptegu

# location of intermediate PHOUR_OUT text file
export PHOUR_OUT=$OUTPUT_DIR/phour.ptegu.2011ek.txt

export YEAR=2011

export PROMPTFLAG=N
$SCRIPTS/aermod/convert_phour

$SCRIPTS/aermod/ptegu.pl
