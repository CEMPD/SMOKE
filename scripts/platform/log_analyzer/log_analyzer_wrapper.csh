#!/bin/csh
#SBATCH --export=NONE

foreach sector (onroad)
foreach rate (RPH)

set scriptdir = /work/EMIS/smoke/smoke4.9/scripts/log_analyzer
set logdir = /work/MOD3DEV/callen05/P105/P105_2020ha_CRACMM/intermed/$sector/$rate/logs

$scriptdir/log_analyzer.py -k $scriptdir/known_messages.txt --list_unknowns -e 1 -f $scriptdir/log_analyzer_wrapper_output_${sector}_${rate}.txt -d $logdir

end
end
