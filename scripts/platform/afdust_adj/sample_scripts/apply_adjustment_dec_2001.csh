#!/bin/csh


set ndays = 10
set ct = 1
set dateout = 20011222
set dates = ( 20021207 20021208 20021224 20021225 20021226 20021203 20021203 \
              20021207 20021208 20021202 )

while ($ct <= $ndays)
set yy = `echo $dateout |cut -c3-4`
set mm = `echo $dateout |cut -c5-6`
set dd = `echo $dateout |cut -c7-8`
set datein = $dates[$ct]
setenv WORK /amber/work/pou/afdust
setenv INFILE $WORK/emis_mole_afdust_meth2_specPMFINE_${datein}_36US1_cmaq_cb05_soa_B71_2002af.ncf
setenv METCRO2D $WORK/METCRO2D_${yy}${mm}${dd}
setenv OUTFILE $WORK/emis_mole_afdust_meth2_adjusted_specPMFINE_${dateout}_36US1_cmaq_cb05_soa_B71_2002af.ncf
setenv LOGFILE apply_precip_adj.$dateout.log
if ( -e $LOGFILE) then
  rm -f $LOGFILE
endif
if ( -e $OUTFILE) then
  rm -f $OUTFILE
endif
/amber/home/pou/codes/afdust/src/apply_precip_adj.x
@ ct = ($ct + 1)
@ dateout = ($dateout + 1)

end
