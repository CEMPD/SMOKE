#!/bin/bash

# location of PHOUR file created by Smkinven
export PHOUR01=$INVOPD/phour.ptegu.jan.2011ek.ncf
export PHOUR02=$INVOPD/phour.ptegu.feb.2011ek.ncf
export PHOUR03=$INVOPD/phour.ptegu.mar.2011ek.ncf
export PHOUR04=$INVOPD/phour.ptegu.apr.2011ek.ncf
export PHOUR05=$INVOPD/phour.ptegu.may.2011ek.ncf
export PHOUR06=$INVOPD/phour.ptegu.jun.2011ek.ncf
export PHOUR07=$INVOPD/phour.ptegu.jul.2011ek.ncf
export PHOUR08=$INVOPD/phour.ptegu.aug.2011ek.ncf
export PHOUR09=$INVOPD/phour.ptegu.sep.2011ek.ncf
export PHOUR10=$INVOPD/phour.ptegu.oct.2011ek.ncf
export PHOUR11=$INVOPD/phour.ptegu.nov.2011ek.ncf
export PHOUR12=$INVOPD/phour.ptegu.dec.2011ek.ncf

# location of AERMOD-source reports created by Smkreport
export REPORT=$REPSCEN/pointgt.aermod.ptegu.12US2.2011ek.2011001.rpt
export REP_XWALK=$REPSCEN/pointgt.aermod.ptegu.12US2.2011ek.2011001.rpt_src_crosswalk.txt
export REP_SRC=$REPSCEN/point.aermod.ptegu.2011ek.rpt

# locations of temporal profile files
export PTPRO_MONTHLY=$GE_DAT/noncem_2014_tpro_monthly_csv_19jul2016_v0
export PTPRO_DAILY=$GE_DAT/noncem_2014_tpro_daily_07_21_2016_19jul2016_v0
export PTPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
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
