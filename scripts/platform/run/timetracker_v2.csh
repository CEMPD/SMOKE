#!/bin/csh -f

# Script to put date/time information about SMOKE to a file
# If file already exists, it will look for entries in the file
# that are already present, and overwrite those lines with
# the new data.

# Arguments:
#   Intialize file?
#   Output file name
#   Start date/time
#   Program
#   Processing date

# Output format: CSV file
# Output columns:
#   Run type
#   Job Name (if applicable)
#   Case abbrev
#   Sector
#   Program
#   Processing date
#   Start date
#   Start time
#   End date
#   End time

# Get end time (want to do this immediately).
set enddt =  `date +%m/%d/%Y,%T`

# parse out script arguments
switch ( $#argv )
   case 2:
      set initialize = $argv[1]
      set outfile    = $argv[2]
   breaksw

   case 5: 
      set initialize = $argv[1]
      set outfile    = $argv[2]
      set startdt    = $argv[3]
      set program    = $argv[4]
      set procdate   = $argv[5]
   breaksw

   default 
      echo "SCRIPT ERROR: Incorrect number of arguments for timetracker.csh"
      echo "              The correct sytax for calling that script is:"
      echo " "
      echo "    timetracker.csh <initialize?> <outfile> <start date/time> <program> <processing date>"
      echo " "
      echo "              If <initialize?> = Y, the the first two arguments are expected.
      echo "              Otherwise, all five arguments are expected.
      exit ( 1 )
   breaksw

endsw

# If 

# If setting for intializing timelog file is set, then...
if ( $initialize == Y ) then

   # Whether the time tracking log is overwritten depends on whether 
   #    the AUTO_DELETE_LOG option is set in the main script.
   # If the timelog file exists, then...
   if ( -e $TIMELOG ) then

      # If autodelete feature is turned on, try to delete log file  
      if ( $AUTO_DELETE_LOG == Y ) then
         echo "SCRIPT NOTE: Deleting existing time log file:"
         echo "             $TIMELOG"
         echo "             Since file exists and AUTO_DELETE_LOG = Y"

         /bin/rm -rf $TIMELOG

         if ( $status == 1 ) then 
            echo "SCRIPT ERROR: Could not delete time log prior to writing it"
            exit ( 1 )
         endif

      # If autodelete feature is not turned on...
      else
         # Try to delete log with using the prompt.
         echo "Delete existing time log file?  If no, script will end."
         /bin/rm -i $TIMELOG

         # If the file still exists, then exit with an error
         if ( -e $TIMELOG ) then
            echo "SCRIPT ERROR: Time log needs to be removed prior to rerunning"
            echo "              this script!"
            exit ( 1 )
         endif
      endif
   endif

   # If haven't yet exited, then initialize timelog file with the header row
   echo RUN_TYPE, JOB_NAME, CASE, SECTOR, PROGRAM, PROCDATE, BEGDATE, BEGTIME, ENDDATE, ENDTIME > $outfile

# If not initializing file, then writing out timing info...
else

   if ( $#argv != 5 ) then
      echo "SCRIPT ERROR: timetracker.csh run with Initialize? = N, but without 5 arguments"
      exit ( 1 )
   endif

   # Separate date and times into separate variables for output
   set startdate = `echo $startdt | cut -d, -f1`
   set starttime = `echo $startdt | cut -d, -f2`
   set enddate   = `echo $enddt | cut -d, -f1` 
   set endtime   = `echo $enddt | cut -d, -f2`

   # Reformat processing date into mm/dd/yyyy
   if ( $procdate != inv ) then
      set dd   = `echo $procdate | cut -c7-8`
      set mm   = `echo $procdate | cut -c5-6`
      set yyyy = `echo $procdate | cut -c1-4`

      set procdate = $mm/$dd/$yyyy
   endif

   # Check to see if there is already a line present with the primary keys
   #    for this dataset, and set internal flag accordingly.
   # Added a cut to the beginning to only get the first six fields as the
   #   last four are not relevant to this check and can produce false positives 
   #   when running timetracker for the same procdate as the current real date.
   #   Changed on 5/17/2010 - bte 
   #   Changed line 133 to use # of words in row rather than status for verification.  Reduces false positives.  3/21/2013 - bte
   set nrow = `cut -d',' -f1-6 "$TIMELOG" | /bin/grep -n "$SECTOR" | grep "$EMF_JOBNAME" | grep ", $program," | grep "$procdate" | cut -d":" -f1` 
   echo $nrow

   # If grep found a matching line, then delete it
   if ( $#nrow != 0 ) then
      # If just one entry in nrow, then delete the line
      if ( $#nrow == 1 ) then
         sed "${nrow}d" $TIMELOG > $TIMELOG.tmp
         mv -f $TIMELOG.tmp $TIMELOG
         echo "SCRIPT NOTE: timetracker is replacing line $nrow of the TIMELOG file"

      # Otherwise, given an error and abort with failure
      else
         echo "SCRIPT ERROR: timetracker script found multiple existing entries"
         echo "              for the following primary keys in the file already:"
         echo "      Sector   = $SECTOR"
         echo "      Job Name = $EMF_JOBNAME"
         echo "      Program  = $program"
         echo "      Run Date = $procdate"
         echo "      This should not happen, so script does not know how to"
         echo "      replace the entries for these keys"
         exit ( 1 )
      endif

   endif       

   # Create string for this line of data and store it
   # If running this from the EMF, then include that information on each line.
   # Otherwise, indicated that it's interactive
   if ( $?EMF_JOBID ) then
       echo EMF, \"$EMF_JOBNAME\", \"$CASE\", \"$SECTOR\", $program, $procdate, \
            $startdate, $starttime, $enddate, $endtime >> $outfile

   else

       echo "interactive, n/a," $CASE, $SECTOR, $program, $procdate, \
            $startdate, $starttime, $enddate, $endtime >> $outfile

   endif

endif

exit ( 0 )

