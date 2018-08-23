#!/bin/csh -f
# Combine the weekly hourly local rungroup reports into single rungroup files for processing

set replist = (`ls /work/EMIS/users/cvy/sas/WO150.8/onroad/nata_hourly/2014fd_nata_onroad_hourly_*_localtime_rungroup.csv`)

# PM2_5 based profile groupings
foreach rungroup (HDON4 HDOFF12 HOTEL4 HDOFF4)
    echo $rungroup
    set outfile = ${rungroup}_localtime.csv
    echo "fips,rungroup,date,hour,emis" > $outfile
    foreach rep ($replist)
        grep "${rungroup}" $rep | cut -d',' -f1-5 >> $outfile
    end
end

# Benzene based profile groupings
foreach rungroup (LDON4 LDOFF12)
    echo $rungroup
    set outfile = ${rungroup}_localtime.csv
    echo "fips,rungroup,date,hour,emis" > $outfile
    foreach rep ($replist)
        grep "${rungroup}" $rep | cut -d',' -f1-4,6 >> $outfile
    end
end


