#!/bin/csh


set ndays = 31
set ct = 1
set dateout = 20020701
set dates = ( 20020708 20020709 20020709 20020704 20020705 20020713 20020714 \
              20020708 20020709 20020709 20020709 20020709 20020713 20020714 \
	      20020708 20020709 20020709 20020709 20020709 20020713 20020714 \
	      20020708 20020709 20020709 20020709 20020709 20020713 20020714 \
	      20020708 20020709 20020709 )

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
