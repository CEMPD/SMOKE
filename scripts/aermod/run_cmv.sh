#!/bin/bash

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

# location of source report created by Smkreport
export REPORT=$REPSCEN/rep_cmv_12_2016fh_16j_aermod.txt

# locations of source groups data files
export GROUP_PARAMS=$SCRIPTS/aermod/nonpoint_rungroup_stack_parameters_11jun2020_v4.csv
export SOURCE_GROUPS=$SCRIPTS/aermod/hem_aermod_groups_2017_11jun2020_nf_v2.csv
export PORT_POLYGONS=$SCRIPTS/aermod/cmv_ports_v2_040718_fornata_12may20_02jun2020_v0.csv

export YEAR=2017
export HOURLY_POLL=VOC

$SCRIPTS/aermod/cmv.pl

# location of cmv PHOUR files created by Smkinven
export PHOUR01=$INVOPD/phour_cmv_12_jan_2016fh_16j.ncf
export PHOUR02=$INVOPD/phour_cmv_12_feb_2016fh_16j.ncf
export PHOUR03=$INVOPD/phour_cmv_12_mar_2016fh_16j.ncf
export PHOUR04=$INVOPD/phour_cmv_12_apr_2016fh_16j.ncf
export PHOUR05=$INVOPD/phour_cmv_12_may_2016fh_16j.ncf
export PHOUR06=$INVOPD/phour_cmv_12_jun_2016fh_16j.ncf
export PHOUR07=$INVOPD/phour_cmv_12_jul_2016fh_16j.ncf
export PHOUR08=$INVOPD/phour_cmv_12_aug_2016fh_16j.ncf
export PHOUR09=$INVOPD/phour_cmv_12_sep_2016fh_16j.ncf
export PHOUR10=$INVOPD/phour_cmv_12_oct_2016fh_16j.ncf
export PHOUR11=$INVOPD/phour_cmv_12_nov_2016fh_16j.ncf
export PHOUR12=$INVOPD/phour_cmv_12_dec_2016fh_16j.ncf

# location of intermediate source list created by cmv.pl
export CMV_SRC_LIST=$OUTPUT_DIR/temporal/CMV_12_source_list.csv

export PROMPTFLAG=N

# remove any existing files
rm -f $OUTPUT_DIR/temporal/CMVU_12_*_hourly.csv
rm -f $OUTPUT_DIR/temporal/CMVP_12_*_hourly.csv

$SCRIPTS/aermod/cmv_output_hourly
