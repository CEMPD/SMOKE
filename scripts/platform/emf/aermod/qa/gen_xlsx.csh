#!/bin/csh -f

set aerpath = /work/EMIS/em_v8/2017_HAPCAP/2017_hapcap_aermod/aermod
set case = 2017_hapcap_aermod

#set groups = (NPLO9AK NPHI9AK NONRD9AK)
set groups = (point_combined)

foreach rg ($groups)
    $aerpath/qa/scripts/gen_xlsx.py $rg $case $aerpath
end
