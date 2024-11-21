#!/bin/tcsh -f
#
# This script automatically copies model-ready emissions to /asm for backup
# purposes. Both sector and merged emissions can be backed up using this
# script.
#
# Script created by : C. Allen, CSC, December 2008
#
#*********************************************************************
# 07 Nov 2013 update: remove --template options from aput commands,
# no longer needed as of today
#
# 05 Dec 2013: added support for a new optional parameter:
# ASM_TEMPLATE = use this ASM template, e.g. "romo:bigcost"
#   As of today, romo:smallcost creates a second copy that lives for 6 months;
#   romo:bigcost does the same, but for 12 months. For 2018ed and 2011ed_ussa,
#   we were asked to use the romo:bigcost template for all archives.
# Also added support for alternate ASM archive directories. For example,
#   for tr_o3, we're required to archive model-ready files under
#   /asm/ROMO/transportrules/[project]/smoke_out/ instead of the usual
#   /asm/ROMO/[project]. This will be determined on the fly by /asm/ROMO/README.
#   This does not affect intermediate files; those will go under 
#   /asm/ROMO/em_[platform]/[project]/ regardless.
#
# 01 Feb 2016: Converted to sol, which required the following changes:
# 1) Sol batch nodes cannot aput, so instead ssh to interactive node atmost for
#    all aput commands
# 2) Sol backups work best when the aput commands are queued (aput -q), which
#    means we can no longer generate the report confirming that all files have
#    been archived
# In addition, we are now archiving the reports/inv and reports/smkmerge directories.
# Reports under reports/inv are zipped on disk and then archived one-by-one.
# Reports under reports/smkmerge are tarred, one .tar per sector, then zipped and archived.

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

# Initialize exit status
set exitstat = 0

switch ( $#argv )
   case 0:
      echo "SCRIPT ERROR: Script requires an argument for the grid name."
      echo " "
      echo "  This script expects to be called with one of the following argument lists:"
      echo "     <grid abbrv>"
      echo "     <grid abbrv> <label>"
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list"
      echo " "
      echo "  Example:"
      echo "     <script name> 36US1"
      echo "              This example runs the script for the 36US1 grid."
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: asm backup script did not receive any arguments" -t "e"
      exit( 1 )
endsw

# Get the first option for the grid abbreviation
setenv GRID "$argv[1]"

# Set TLABEL if available
if ( $#argv >= 2 ) then
   setenv TLABEL "$argv[2]"
endif

# Set ASM_TEMPLATE default, if not defined
# See description of this parameters in the 05 Dec 2013 revision note 
#   at the top of this script
if (! $?ASM_TEMPLATE) setenv ASM_TEMPLATE ""
echo "ASM template parameter = $ASM_TEMPLATE"

## Set up scripting environment variables prior to calling the Assigns file
setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

source $ASSIGNS_FILE

## 7 May 2014: Added support for sensitivity tool runs
if ( ! $?RUN_SOURCESENS ) then
    setenv RUN_SOURCESENS N
endif

## If source sensitivity, need source sector override file
if ( $RUN_SOURCESENS == Y ) then
    echo "Running archives for source sensitivity"
    set sector_override = $IMD_ROOT/mrggrid/source_sector_override_${CASE}_${GRID}.txt
    if ( ! -e $sector_override ) then
	echo "SCRIPT ERROR: source sector override file doesn't exist "
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: " -t "e"
	set exitstat = 1
	goto end_of_script
    endif

    ## Unique list of sectors to override (ie. set to source sens case over parent case)
    set source_sense_sectors = `sort -u $sector_override`
    echo "SCRIPT NOTE: Sectors to archive from source sensitivity case: $source_sense_sectors"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Sectors to archive from source sensitivity case: $source_sense_sectors"
endif

### List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
#set sectorlist_parser  = $SCRIPTS/run/sectorlist_parser.py

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

## Set naming label
set namelabel = ${CASE}_${GRID}
if ( $?TLABEL ) then
  set namelabel = ${namelabel}_$TLABEL
endif

# Months list; used later for file counting purposes
set allmos = (01 02 03 04 05 06 07 08 09 10 11 12)

# Set ARCHIVE_ALL_SECTORS to Y if not defined already
# If set to N, sectors from older cases will not be zipped/archived,
# under the assumption that they have already been archived. If this assumption
# is false, then the file count check will fail and the job will fail to alert the user
if ( ! $?ARCHIVE_ALL_SECTORS ) then
   setenv ARCHIVE_ALL_SECTORS Y
endif

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for $SECTOR" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
#echo $sectorlist_parser >> $LOGS/helper_scripts_list$suffix

## Define archival report file, which will list each sector that was archived along with file counts
setenv ASMLOG $REPOUT/asm_backup/archive_report_$namelabel.txt
if (! -e $REPOUT/asm_backup) then
   mkdir -p $REPOUT/asm_backup
endif
if (-e $ASMLOG) then
   rm -f $ASMLOG
endif

# Header for ASMLOG, which will be in .csv format
echo "case, sector, grid, speciation, files_on_disk, files_archived" >> $ASMLOG

## Set Time Log filename and initialize file
setenv TIMELOG $LOGS/timelog_$namelabel.txt

## If TIMELOG_YN = N, don't do timetracker (default is Y)
if (! $?TIMELOG_YN) setenv TIMELOG_YN Y

# Only initialize TIMELOG if it doesn't already exist, since the timeracker 
#   can now delete/add entries to prevent duplicates
if ( ! -e $TIMELOG && $TIMELOG_YN != N ) then 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Initializing Time Log" -x $timetracker  ## log w/ EMF server
   $timetracker Y $TIMELOG
   if ( $status != 0 ) then
	echo "ERROR: running timetracker"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running timetracker to initialize time log" -t "e" -x $timetracker
	exit ( 1 )
   endif
endif

# If timelog turned off, unset TIMELOG variable, which will turn timetracker off in the rest of the scripts
if ($TIMELOG_YN == N) unsetenv TIMELOG

## Check that environment variables for required inputs are set and exist
if ( ! $?SECTORLIST ) then
   echo "SCRIPT ERROR: Environment variable SECTORLIST is not defined"
   echo "              but is required by asm_backup_emf.csh"
   exit ( 1 )
else
   if ( ! -e $SECTORLIST ) then
      echo "SCRIPT ERROR: SECTORLIST file is not found:"
      echo "              "$SECTORLIST
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: SECTORLIST file is not found" -t "e"
      exit( 1 )
   endif
endif

## Set the EMF_PERIOD to the year
setenv EMF_PERIOD $YEAR

## Grab group name and (if necessary) romo subproject from EMF_QUEUE_OPTIONS
#  EMF_QUEUE_OPTIONS will probably never exceed 16 fields
#  First step: determine where in EMF_QUEUE_OPTIONS '-q' and '-A' appear
setenv GROUP_NAME "null"
setenv ROMO_SUBPROJ "null"
set found_q = 0
set found_A = 0

# Search EMF_QUEUE_OPTIONS field-by-field to find each
# The idea is, if '-q' or '-A' is found, then set a flag
# Then next iteration, set the following field to the appropriate variable
## C. Allen: There is probably an easier/better way to do this, but this works just fine
foreach x (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

   set tmp = `echo $EMF_QUEUE_OPTIONS | cut -f$x -d' '`

   if ( $found_q ) then
      setenv GROUP_NAME $tmp
      set found_q = 0
      continue
   endif

   if ( $found_A ) then
      setenv ROMO_SUBPROJ $tmp
      set found_A = 0
      continue
   endif
   
   if ( "$tmp" == "-q" ) then
      set found_q = 1
   endif
   if ( "$tmp" == "-A" ) then
      set found_A = 1
   endif
   
end # foreach

# If couldn't find group name, exit
if ( $GROUP_NAME == "null" ) then
    echo "SCRIPT ERROR: Unable to extract group name from EMF_QUEUE_OPTIONS"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Unable to extract group name from EMF_QUEUE_OPTIONS" -t "e"
    exit ( 1 )
endif

## UPDATE 3/16/2009: If running on garnet, the queue name is workq, but we obviously can't archive to
#  /asm/WORKQ/ or use a "workq" ASM template. The assumption is that any garnet job is a romo job,
#  since garnet is supposed to be an OAQPS-only machine.
#  To determine whether or not we're on garnet, look up the hostname. Grab the first six characters
#  since there are multiple garnet nodes (garnet01, garnet02, etc).
set hostname_16 = `hostname | cut -c1-6`
if ( $hostname_16 == garnet ) then
   setenv GROUP_NAME romo
endif

# If group = romo and couldn't find romo subproject, exit
if ( $GROUP_NAME == "romo" && $ROMO_SUBPROJ == "null" ) then
    echo "SCRIPT ERROR: Unable to extract romo subproject from EMF_QUEUE_OPTIONS"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Unable to extract romo subproject from EMF_QUEUE_OPTIONS" -t "e"
    exit ( 1 )
endif

## Define output directory root for sector backups
#  18 May 2016: Now that we're running everything on /sol, we can use PROJECT_ROOT for this.
#  Group name is either 3rd or 4th field, delimited by /, of PROJECT_ROOT (depending on whether it's defined as /sol/work or just /work)

#set group_name_caps = `echo $GROUP_NAME | tr '[:lower:]' '[:upper:]'`

set group_name_caps = `echo $PROJECT_ROOT | cut -f3 -d'/'`
if ($group_name_caps == "work") set group_name_caps = `echo $PROJECT_ROOT | cut -f4 -d'/'`

set asm_root = /asm/$group_name_caps

set group_name_lower = `echo $group_name_caps | tr '[:upper:]' '[:lower:]'`

# Model-ready files go to /asm/ROMO if group is EMIS
if ($group_name_caps == "EMIS") then
  set asm_root2 = /asm/ROMO
else
  set asm_root2 = $asm_root
endif

## If asm_root doesn't exist, exit
## Can't do this on sol because /asm isn't visible from batch nodes
#if ( ! -e $asm_root ) then
#    echo "SCRIPT ERROR: Directory $asm_root does not exist or is unreachable"
#    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Directory $asm_root does not exist or is unreachable" -t "e"
#    exit ( 1 )
#endif

# Get number of lines in sector file
set nlines = `cat $SECTORLIST | wc -l`

# Can't read /asm/ROMO/README from batch nodes, so copy it locally first
ssh atmost cp /asm/ROMO/README ~/asm_readme_temp

## Loop through sectors to back up using SECTORLIST
set nl = 0
while ( $nl < $nlines )

   @ nl = $nl + 1
   if ( $nl > 1 ) then   # skip header line

      # Grab sector names and cases; no need for other fields unless CASE = "none"
      # Since not all runs are annual, we don't want to hardwire file counts using the merge type
      # Update 8/25/10 (C. Allen) reverted to old method of parsing SECTORLIST lines,
      #   without using sectorlist parser python script, in order to allow for situation
      #   where we want to merge sectors with like names from different cases.
      set line = ( `head -$nl $SECTORLIST | tail -1 | sed 's/""/"null"/g' | sed 's/,/ /g' | sed 's/"//g'` )
      set sector_name  = $line[1]
      set sector_case  = $line[2]

     ## 7 May 2014: Added support for sensitivity tool runs
     ## If this is a source sensitivity, override 
     ## the sectorlist's case (typically PARENT_CASE) with
     ## source sensitivity case if this sector is found in the 
     ## source sector override file
     if ( $RUN_SOURCESENS == Y ) then
     	grep -xqi $sector_name $sector_override
     	if ( $status == 0 ) then
	    ## found this exact sector name in override file
	    echo "Overriding case for sector $sector_name, setting to $CASE"
	    set sector_case = $CASE
     	else
	    echo "Not overriding case sector $sector_name, set to $sector_case"
     	endif
     endif

      ## C. Allen update 05/15/2009: Support for varying speciations modified due to new "mergesector" column
      #  speciation; if it's not there, then set to default $SPC
      #  Updated again on 8/25/2010 to get away from sectorlist_parser.py
      set sector_spc = $SPC
      if ( $#line >= 7 ) then
	 if ( $line[7] != "null" ) set sector_spc = $line[7]
      endif

      ## project root; if it's not there, then set to default $PROJECT_ROOT
      set sector_project    = $PROJECT_ROOT
      if ( $#line >= 9 ) then
         if ( $line[9] != "null" && $line[9] != "") set sector_project = $line[9]
      endif
      
     set firstdir = `echo $sector_project | cut -f2 -d'/'`
     if ($firstdir == "sol") then
       set projectx = `echo $sector_project | cut -f6 -d'/'`
       set tempx = `echo $sector_project | cut -f1,3- -d'/'`
       set sector_project = $tempx # remove /sol from directory names; otherwise aput won't work
     else
       set projectx = `echo $sector_project | cut -f5 -d'/'`
     endif
     echo "projectx = $projectx"
     
     # Remove "sol" from OUT_ROOT, if necessary; otherwise aput won't work
     set tempx = `echo $OUT_ROOT | cut -f2 -d'/'`
     if ($tempx == "sol") then
       set out_rootx = `echo $OUT_ROOT | cut -f1,3- -d'/'`
     else
       set out_rootx = $OUT_ROOT
     endif
     echo "out_rootx = $out_rootx"

     # Remove "sol" from IMD_ROOT, if necessary; otherwise aput won't work
     set tempx = `echo $IMD_ROOT | cut -f2 -d'/'`
     if ($tempx == "sol") then
       set imd_rootx = `echo $IMD_ROOT | cut -f1,3- -d'/'`
     else
       set imd_rootx = $IMD_ROOT
     endif
     echo "imd_rootx = $imd_rootx"

      # For "one-day" emissions files such as ocean chlorine and volcanic mercury, skip
      if ( $sector_case == none ) then
         
         echo "SCRIPT NOTE: Skipping $sector_name emissions file; this is a single file which has likely already been archived"
      
      else

         # If case is set using an environment variable, evaluate it
         # All new SECTORLISTs should be explicit, but older ones may use $CASE, for example
         echo $sector_case | /bin/grep -q \\$   # Search for $ sign in setting to ID environment var
         if ( $status == 0 ) then
            set tmpcase = `echo $sector_case | sed 's/\$//'`
            set sector_case = `env | /bin/grep "^$tmpcase=" | cut -d\= -f2`
         endif

         set asm_dir1 = $asm_root/em_${PLATFORM}/$projectx/$sector_case/premerged
	 set asm_dir_rep = $asm_root/em_${PLATFORM}/$projectx/$sector_case/reports
	 
	 # 5 Dec 2013: Search /asm/ROMO/README to see if [projectx] should be archived in 
	 #   /asm/ROMO/[top-level project]/[projectx] instead of /asm/ROMO/[projectx]
	 #  
	 grep ^$asm_root2/$projectx ~/asm_readme_temp > /dev/null
	 if ($status == 0) then # alternate ASM directory found
	    set new_asm_root = `grep ^$asm_root2/$projectx ~/asm_readme_temp | cut -f2 -d'>' | sed 's/ //g'`
	    set asm_dir2 = $new_asm_root/smoke_out/$sector_case/$GRID/$sector_spc
	 else # normal ASM directory
	    set asm_dir2 = $asm_root2/$projectx/smoke_out/$sector_case/$GRID/$sector_spc
	 endif

         ## C. Allen update 02/03/2010: The removal of Smkinven-created $sector_project/$sector_case/intermed/$CASE/$SECTOR/import_tmp* temporary files 
	 #    doesn't work properly in smk_run_v7.csh, so instead we will remove these temporary files here.
	 rm -fv $sector_project/$sector_case/intermed/$sector_name/tmp/import_tmp*

         ## C. Allen update 04/14/2009: Split apart by month to help avoid "argument list too long" errors

         if ( $ARCHIVE_ALL_SECTORS == N && $sector_case != $CASE ) then

            echo "SCRIPT NOTE: Skipping $GRID emissions in $sector_case/$sector_name due to ARCHIVE_ALL_SECTORS = N"

            # Are files already on asm zipped or unzipped? In some cases, both zipped and unzipped copies may exist on ASM,
	    # and double-counting them would result in a Fail. (This should never happen with anything archived using this script,
	    # but there could be some instances of this floating around in other ASM directories).
	    # C. Allen 9/15/11: I don't think this is needed anymore, and I'm not sure it worked properly in the first place. Commenting out.
	    
#	    # This foreach loop, repeated several times throughout the script, splits up file lists by month to avoid
#	    # "argument list too long" errors with ls, gzip, and aput.
#	    set asmfiles_zipped = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $asm_dir1/${sector_name}/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf.gz | wc -l`
#  	       @ asmfiles_zipped = $asmfiles_zipped + $i1
#	       set i1 = `ls -1 $asm_dir2/${sector_name}/inln_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf.gz | wc -l`
#  	       @ asmfiles_zipped = $asmfiles_zipped + $i1
#	    end # foreach
#	    if ( $asmfiles_zipped > 0 ) then
#	       set ext = ncf.gz
#	    else
#	       set ext = ncf
#	    endif

#	    # 27 Sep 2013: Changed "inln" to "*", so that SGINLN and STACK_GROUPS_OUT files for source apportionment
#	    # are counted
#            set diskfiles = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $sector_project/$sector_case/intermed/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf* | wc -l`
#  	       @ diskfiles = $diskfiles + $i1
#	       set i1 = `ls -1 $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ diskfiles = $diskfiles + $i1
#	    end # foreach
#
#	    set asmfiles = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $asm_dir1/${sector_name}/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf* | wc -l`
#  	       @ asmfiles = $asmfiles + $i1
#	       set i1 = `ls -1 $asm_dir2/${sector_name}/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ asmfiles = $asmfiles + $i1
#	    end # foreach
#
#            if ( $diskfiles > $asmfiles ) then
#               echo "ERROR: Only $asmfiles of $diskfiles already archived on ASM for sector $sector_name"
#               $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Only $asmfiles of $diskfiles already archived on ASM for sector $sector_name" -t "e"
#               set exitstat = 1
#            else
#               echo "SCRIPT NOTE: $asmfiles of $diskfiles already archived to ASM for sector $sector_name"
#   	    endif

         else

            echo "SCRIPT NOTE: Archiving $GRID emissions in $sector_case/$sector_name..."

            # Check RUNSET to see if this sector is to be skipped
            # This can be done to avoid lengthy gzips of files that have previously been archived
	    ## C. Allen: I originally considered adding this feature, but decided not to

            # gzip emis all files before copying to asm, but NOT INLN files (9/15/11 update, otherwise we may interfere with an ongoing model run)
	    # 01 Feb 2016: also zip Smkmerge reports and add to .tar file that will be archived
	    foreach mon ($allmos)
               gzip -fv $sector_project/$sector_case/premerged/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf
#               gzip -fv $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/inln_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf
               pushd $sector_project/$sector_case/reports/smkmerge
	       gzip -fv $sector_name/rep_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*txt
	       tar -rvf ${sector_name}_${GRID}.tar $sector_name/rep_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*txt.gz
	       popd
            end # foreach
#	    echo "$OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/inln_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf"

            # Sometimes, emissions files are softlinked to other directories, necessary for a variety of reasons.
	    # In these instances, the files will not be able to be zipped, but they can still be archived.
	    # Thus, we must allow for the case where the files are still unzipped on disk.
	    # C. Allen 9/15/11: I don't think this is needed anymore. Commenting out.
	    
#	    set remaining = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $asm_dir1/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf | wc -l`
#  	       @ remaining = $remaining + $i1
#	       set i1 = `ls -1 $asm_dir2/$sector_name/inln_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf | wc -l`
#  	       @ remaining = $remaining + $i1
#	    end # foreach
#
#	    if ( $remaining > 0 ) then
#	       set ext = ncf
#	    else
#	       set ext = ncf.gz
#	    endif

            set startdt = `date +%m/%d/%Y,%T` # for timetracker
      
            # Script supports use of romo subtemplates using the subproject name, but not for subtemplates for
	    # non-romo jobs (if such subtemplates exist)
	    # 12/5/13: Templates not used anymore, except when explicitly set via ASM_TEMPLATE parameter (e.g. romo:bigcost).
	    #   Commenting out, except for a section to support non-romo templates.
	    
#            if ( $GROUP_NAME == romo ) then
#	       # At one time, there was a discrepancy between the names of the queue subprojects and the corresponding ASM templates.
#	       # As of 12/18/2012, the names are now synced up.
#	       if ( $ROMO_SUBPROJ == plateval ) then
#   	          set asm_template = romo:plat_eval
#	       else if ( $ROMO_SUBPROJ == naaqs ) then
#   	          set asm_template = romo:naaqs_di
#	       else if ( $ROMO_SUBPROJ == mobile ) then
#   	          set asm_template = romo:naaqs_ss
#	       else
#   	          set asm_template = romo:${ROMO_SUBPROJ}
#	       endif
#	    else
#	       set asm_template = ${GROUP_NAME}
#	    endif
            if ($ASM_TEMPLATE == "" && $group_name_lower != "romo" && $group_name_lower != "emis") then
               setenv ASM_TEMPLATE $group_name_lower
            endif
	    
            # Set command line --template parameter, unless no template
	    if ($ASM_TEMPLATE != "") then
	       set asm_temp = "--template=$ASM_TEMPLATE"
	    else
	       set asm_temp = ""
	    endif
	    
            # If these files have already been backed up to ASM, aput will skip over these files
            foreach mon ($allmos)
	    
               # If model-ready files are on /terra, then aput won't work. Must copy to /garnet and archive from there.
	       # 9/15/11 update C. Allen: if on terra, might as well gzip the garnet copy before archiving
	       #                          if not on terra, archive directly if already zipped. if not already zipped: copy, zip, then archive
	       # 9/27/13 update C. Allen: changed "inln" to "*" so that SGINLN and STACK_GROUPS_OUT files for source
	       #   apportionment are archived
               if ($OUT_ROOT == "/terra/work/ROMO" || $OUT_ROOT == "/sol/work/ROMO") then           
	          echo "$asm_dir1/$sector_name"
  	          ssh atmost aput -q -a $asm_dir1/$sector_name $asm_temp $sector_project/$sector_case/premerged/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf*
	          cp -pv $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* $sector_project/$sector_case/premerged
	          cp -pv $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/stack_groups_${sector_name}_*ncf* $sector_project/$sector_case/premerged
		  gzip -v $sector_project/$sector_case/premerged/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf
		  ssh atmost aput -q -d -a $asm_dir2/$sector_name $asm_temp $sector_project/$sector_case/premerged/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf*
		  ssh atmost aput -q -d -a $asm_dir2/$sector_name $asm_temp $sector_project/$sector_case/premerged/stack_groups_${sector_name}_*ncf*

		  # Archive onroad RPD/RPP/RPV source apportionment files
		  if (-e $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/RPD) then
		    foreach rpx (RPD RPP RPV RPH)
		      cp -pv $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/$rpx/*_${rpx}_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* $sector_project/$sector_case/premerged
		      gzip -v $sector_project/$sector_case/premerged/*_${rpx}_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf
  		      ssh atmost aput -q -d -a $asm_dir2/$sector_name/$rpx $asm_temp $sector_project/$sector_case/premerged/*_${rpx}_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf.gz
		    end
		  endif

	       else
  	          ssh atmost aput -q -a $asm_dir1/$sector_name $asm_temp $sector_project/$sector_case/premerged/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf*
		  set t1 = `ls -1 $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf | wc -l`
		  if ($t1 > 0) then
		    cp -pv $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* $sector_project/$sector_case/premerged
		    gzip -v $sector_project/$sector_case/premerged/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf
  		    ssh atmost aput -q -d -a $asm_dir2/$sector_name $asm_temp $sector_project/$sector_case/premerged/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf.gz
		  else
                    ssh atmost aput -q -a $asm_dir2/$sector_name $asm_temp $out_rootx/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf.gz
		  endif
		  
               endif	       
	       
            end # foreach

            # Aput .tar file containing Smkmerge reports
	    echo "archiving $sector_project/$sector_case/reports/smkmerge/${sector_name}_${GRID}.tar"
            ssh atmost aput -q -d -a $asm_dir_rep/smkmerge $asm_temp $sector_project/$sector_case/reports/smkmerge/${sector_name}_${GRID}.tar 
	    
            # update time log
	    if ($TIMELOG_YN != N) then
               $timetracker N $TIMELOG $startdt $sector_name $ESDATE
               if ( $status != 0 ) then
                   echo "ERROR: Problem calling timetracker from asm_backup script"
                   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from asm_backup script" -x $timetracker  -t "e" -p $EMF_PERIOD ## log w/ EMF server
                   exit ( 1 )
               endif 
	    endif

#            # Make sure all files were copied; aput status return isn't always reliable
#            # If there is a mismatch in the number of files, flag an error,
#            # but DO NOT exit the script so that the other sectors are copied
#            set diskfiles = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $sector_project/$sector_case/intermed/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf* | wc -l`
#  	       @ diskfiles = $diskfiles + $i1
#	       set i1 = `ls -1 $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ diskfiles = $diskfiles + $i1
#	       set i1 = `ls -1 $OUT_ROOT/$projectx/smoke_out/$sector_case/$GRID/$sector_spc/$sector_name/RP?/*_RP?_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ diskfiles = $diskfiles + $i1
#	    end # foreach
#
#	    set asmfiles = 0
#	    foreach mon ($allmos)
#	       set i1 = `ls -1 $asm_dir1/$sector_name/emis_mole_${sector_name}_????${mon}??_${GRID}_${sector_spc}_${sector_case}*ncf* | wc -l`
#  	       @ asmfiles = $asmfiles + $i1
#	       set i1 = `ls -1 $asm_dir2/$sector_name/*_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ asmfiles = $asmfiles + $i1
#	       set i1 = `ls -1 $asm_dir2/$sector_name/RP?/*_RP?_${sector_name}_????${mon}??_${GRID}_*${sector_case}*ncf* | wc -l`
#  	       @ asmfiles = $asmfiles + $i1
#	    end # foreach
#
#            if ( $diskfiles > $asmfiles ) then
#               echo "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for sector $sector_name"
#               $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for sector $sector_name" -t "e"
#               set exitstat = 1
#            else
#               echo "SCRIPT NOTE: $asmfiles of $diskfiles successfully copied to ASM for sector $sector_name"
#   	    endif
	    
         endif # ARCHIVE_ALL_SECTORS
            
#         # Write information to ASMLOG
#         echo "$sector_case, $sector_name, $GRID, $sector_spc, $diskfiles, $asmfiles" >> $ASMLOG
      
      endif # case = none
   
   endif # nl > 1
   
end # while nl < nlines
      
# 5 Dec 2013: Search /asm/ROMO/README to see if [projectx] should be archived in 
#   /asm/ROMO/[top-level project]/[projectx] instead of /asm/ROMO/[projectx]
#  
grep ^$asm_root2/$PROJECT ~/asm_readme_temp > /dev/null
if ($status == 0) then # alternate ASM directory found
   set new_asm_root = `grep ^$asm_root2/$PROJECT ~/asm_readme_temp | cut -f2 -d'>' | sed 's/ //g'`
   set asm_dir2 = $new_asm_root/smoke_out/$CASE/$GRID/$SPC
else # normal ASM directory
   set asm_dir2 = $asm_root2/$PROJECT/smoke_out/$CASE/$GRID/$SPC
endif

# Merged emissions must be zipped before running this script
# This policy prevents zipping emissions if a model run is in progress, interrupting the model run
# Ideally, emissions will stay far enough ahead of modeling so that archiving happens before modeling begins
#gzip -fv $OUTPUT/emis_mole_all_*_${GRID}_${SPC}_${CASE}.ncf

set startdt = `date +%m/%d/%Y,%T` # for timetracker
     
      
# Script supports use of romo subtemplates using the subproject name, but not for subtemplates for
# non-romo jobs (if such subtemplates exist)
# 12/5/13: Templates not used anymore, except when explicitly set via ASM_TEMPLATE parameter (e.g. romo:bigcost).
#   Commenting out, and adding a separate section to support non-romo templates

#if ( $GROUP_NAME == romo ) then
#   # At one time, there was a discrepancy between the names of the queue subprojects and the corresponding ASM templates.
#   # As of 12/18/2012, the names are now synced up.
#  if ( $ROMO_SUBPROJ == plateval ) then
#     set asm_template = romo:plat_eval
#  else if ( $ROMO_SUBPROJ == naaqs ) then
#     set asm_template = romo:naaqs_di
#  else if ( $ROMO_SUBPROJ == mobile ) then
#     set asm_template = romo:naaqs_ss
#  else
#      set asm_template = romo:${ROMO_SUBPROJ}
#   endif
#endif
if ($ASM_TEMPLATE == "" && $group_name_lower != "romo" && $group_name_lower != "emis") then
   setenv ASM_TEMPLATE $group_name_lower
endif

# Set command line --template parameter, unless no template
if ($ASM_TEMPLATE != "") then
   set asm_temp = "--template=$ASM_TEMPLATE"
else
   set asm_temp = ""
endif
	    
# If these files have already been backed up to ASM, aput will skip over these files

# If model-ready files are on /terra, then aput won't work. Must copy to /garnet and archive from there.
# 9/15/11 update: might as well gzip the copies, right?
set out_disk = `echo $OUT_ROOT | cut -f2 -d'/'`
if ($out_disk == "terra" || $out_disk == "sol") then      
   cp -pv $OUTPUT/emis_mole_all_*_${GRID}_*_${CASE}*ncf* $OUTPUT/stack_groups_* $IMD_ROOT
   gzip -v $IMD_ROOT/emis_mole_all_*_${GRID}_*_${CASE}*ncf
   ssh atmost aput -q -d -a $asm_dir2 $asm_temp $imd_rootx/emis_mole_all_*_${GRID}_*_${CASE}*ncf.gz
   ssh atmost aput -q -d -a $asm_dir2 $asm_temp $imd_rootx/stack_groups_*
else # if any unzipped files: copy all, gzip copies, aput -d. if no unzipped files, aput directly
   set t1 = `ls -1 $OUTPUT/emis_mole_all_*_${GRID}_*_${CASE}*ncf | wc -l`
   if ($t1 > 0) then
     cp -pv $OUTPUT/emis_mole_all_*_${GRID}_*_${CASE}*ncf* $OUTPUT/stack_groups_* $IMD_ROOT
     gzip -v $IMD_ROOT/emis_mole_all_*_${GRID}_*_${CASE}*ncf
     ssh atmost aput -q -d -a $asm_dir2 $asm_temp $imd_rootx/emis_mole_all_*_${GRID}_*_${CASE}*ncf.gz
     ssh atmost aput -q -d -a $asm_dir2 $asm_temp $imd_rootx/stack_groups_*
   else
     ssh atmost aput -q -a $asm_dir2 $asm_temp $out_rootx/emis_mole_all_*_${GRID}_*_${CASE}*ncf.gz 
     ssh atmost aput -q -a $asm_dir2 $asm_temp $out_rootx/stack_groups_*
   endif
endif	       



# update time log
if ($TIMELOG_YN != N) then
   $timetracker N $TIMELOG $startdt $sector_name $ESDATE
   if ( $status != 0 ) then
      echo "ERROR: Problem calling timetracker from asm_backup script"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from asm_backup script" -x $timetracker  -t "e" -p $EMF_PERIOD ## log w/ EMF server
      exit ( 1 )
   endif  
endif

## Make sure all files were copied; aput status return isn't always reliable
## If there is a mismatch in the number of files, flag an error,
## but DO NOT exit the script so that the other sectors are copied
#set diskfiles = `ls -1 $OUTPUT/emis_mole_all_*_${GRID}_${SPC}_${CASE}.ncf* | wc -l`
#set asmfiles = `ls -1 $asm_dir2/emis_mole_all_*_${GRID}_${SPC}_${CASE}.ncf.gz | wc -l`
#if ( $diskfiles > $asmfiles ) then
#   echo "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for merge"
#   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for merge" -t "e"
#   set exitstat = 1
#else
#   echo "SCRIPT NOTE: $asmfiles of $diskfiles successfully copied to ASM for merge"
#endif
#
## Write information to ASMLOG
#echo "$CASE, merge, $GRID, $SPC, $diskfiles, $asmfiles" >> $ASMLOG
#
## C. Allen updated 04/01/2009 to archive STACK_GROUPS files as well, since these are needed for inline modeling
## This should work for date-specific STACK_GROUPS files, as in ptfire
## The files are small, so I see no need to zip
#set diskfiles = `ls -1 $OUTPUT/stack_groups_*_${GRID}_*.ncf* | wc -l`
#if ( $diskfiles > 0 ) then
#   if ($out_disk == "terra" || $out_disk == "sol") then           
#      cp -pv $OUTPUT/stack_groups_*_${GRID}_*.ncf* $IMD_ROOT
#      ssh atmost aput -q -d -a $asm_dir2 $asm_temp $IMD_ROOT/stack_groups_*_${GRID}_*.ncf*
#   else
#      ssh atmost aput -q -a $asm_dir2 $asm_temp $OUTPUT/stack_groups_*_${GRID}_*.ncf*
#   endif	       
#endif
#set asmfiles = `ls -1 $asm_dir2/stack_groups_*_${GRID}_*.ncf* | wc -l`
#if ( $diskfiles > $asmfiles ) then
#   echo "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for merge"
#   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Only $asmfiles of $diskfiles successfully copied to ASM for STACK_GROUPS" -t "e"
#   set exitstat = 1
#else
#   echo "SCRIPT NOTE: $asmfiles of $diskfiles successfully copied to ASM for STACK_GROUPS"
#endif
#
## Write information to ASMLOG
#echo "$CASE, stack_groups, $GRID, $SPC, $diskfiles, $asmfiles" >> $ASMLOG

# Archive reports/inv

reports_inv:

# Remove "sol" from REPOUT, if necessary; otherwise aput won't work
set tempx = `echo $REPOUT | cut -f2 -d'/'`
if ($tempx == "sol") then
  set repoutx = `echo $REPOUT | cut -f1,3- -d'/'`
else
  set repoutx = $REPOUT
endif
echo "repoutx = $repoutx"

pushd $REPOUT/inv
gzip -v rep*t
tar -rvf inv_state.tar rep*state.txt.gz
tar -rvf inv_state_scc.tar rep*state_scc.txt.gz
tar -rvf inv_county.tar rep*county.txt.gz
tar -rvf inv_county_scc.tar rep*county_scc.txt.gz
tar -rvf inv_${GRID}.tar rep*${GRID}.txt.gz
tar -rvf inv_prof.tar rep*prof.txt.gz
tar -rvf inv_naics_oris.tar rep*naics.txt.gz rep*oris.txt.gz
ssh atmost aput -q -d -a $asm_root/em_${PLATFORM}/$PROJECT/$CASE/reports/inv $repoutx/inv/*.tar
popd

# Label for the end of the script, used during script abort
## C. Allen: This label is never actually used in this script
end_of_script:

## Register time log
#echo "SCRIPT NOTE: Registering time log"
#$EMF_CLIENT -k $EMF_JOBKEY -F $TIMELOG -T "SMOKE time log (External)" -N "SMOKE timelog $namelabel" -O "Timelog $namelabel (External)"

## Register ASM log
#echo "SCRIPT NOTE: Registering archive log"
#$EMF_CLIENT -k $EMF_JOBKEY -F $ASMLOG -T "ASM Archive Report (CSV)" -N "ASM archive log $namelabel" -O "ASM archive log $namelabel"

rm ~/asm_readme_temp

## Ending of script
#
exit( $exitstat )
