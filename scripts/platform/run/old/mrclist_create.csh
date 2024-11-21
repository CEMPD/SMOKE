#!/bin/csh

set indir = /garnet/oaqps/tier3_prop/Tier3Base2005_20110520/output
set filestamp = 20may2011_tier3base2005

set origdir = $cwd

cd $indir

foreach moves_type (rateperdistance ratepervehicle rateperprofile)

  set outfile = $indir/mrclist.$moves_type.$filestamp.summed.lst
  if (-e $outfile) rm $outfile
  
  echo "# FIPS, fuelmonth, EFtable" >> $outfile
  echo "FIPS, fuelmonth, EFtable" >> $outfile

  set filelist = (`ls ${moves_type}_smoke_postprocess_*.summed.csv`)
  
  foreach f ($filelist)
  
    set fips = `echo $f | cut -f4 -d'_'`
    
    set fuelmo = `echo $f | cut -f5 -d'_' | cut -f1 -d'.'`
  
    echo "$fips,$fuelmo,$f" >> $outfile
  
  end # foreach f
  
end # foreach moves_type

cd $origdir
