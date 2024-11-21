#!/bin/tcsh -fx
#
# $Version$
# $Path$
# $Date$
#
# This script sets the dates needed for a given script run and
# stores these days into.
#
# This script should be "sourced" (not run) so that it will set the 
#  "SMK_RUN_DATES" environment variable
#
# Script created by : M. Houyoux, EPA
# Last edited : December, 2006
# Updated by M. Houyoux on 12/6/2006 to correctly set holiday flag for
#     spin-up dates of the previous year.
#
#*********************************************************************

# ----------------------------------------------------------------------------
#
# NOTES:
#
#  Cases for which calling script may need dates:
#	1) Average weekday representing the entire year
#	2) Ave. Weekday, Sat, Sun, Monday for entire year
#	3) Ave. Monday through Sunday for entire year
#	4) Average weekday per season
#	5) Average Weekday, Sat, Sun, Monday per season.
#	6) Average Monday through Sunday per season
#	7) Average weekday per for select months
#	8) Average Weekday, Sat, Sun, Monday per month.
#	9) Average Monday through Sunday per month
#	10) Average weekday per month
#	11) Average Weekday, Sat, Sun, Monday per month.
#	12) Average Monday through Sunday per month
#	13) Run all days individually

# For 1-9, can be with and without holidays specified separately.
# 10 does not need to ID holidays because they would be modeled separately no
#   matter what.

# ----------------------------------------------------------------------------

# Check that controlling variables have been defined, and if not, set
#    default values.

set eflag = 0
if ( $?YEAR ) then
   set year = $YEAR       
   
else
   echo "SCRIPT ERROR: environment variable YEAR is a required input setting"
   echo "       to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1

endif

if ( $?MONTH_ARRAY ) then
   set monset = ( $MONTH_ARRAY )
   
else
   echo "SCRIPT ERROR: environment variable MONTH_ARRAY is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1

endif

if ( $?SPINUP_ARRAY ) then
   set spinup_arr = ( $SPINUP_ARRAY )
   
else
   echo "SCRIPT ERROR: environment variable SPINUP_ARRAY is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1

endif

if ( $?RUN_HOLIDAYS ) then
   if ( $RUN_HOLIDAYS == y ) then
      setenv RUN_HOLIDAYS Y
   endif   
else
   echo "SCRIPT ERROR: environment variable RUN_HOLIDAYS is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1

endif

if ( $?RESTART_JDATE ) then
   set restart = Y
   if ( $RESTART_JDATE == " " ) then
      set restart = N
   endif
else
   set restart = N
   setenv RESTART_JDATE " "
endif

## Note - SPINUP_SET is only used for T_TYPE = "all"

if ( $?T_TYPE ) then
   switch ( $T_TYPE ) 
      case ALL:
      case All:
      case all:
         set t_type = all
      breaksw
      
      case AVEDAY:
      case Aveday:
      case aveday:
         set t_type = aveday
      breaksw
      
      case MWDSS:
      case Mwdss:
      case mwdss:
         set t_type = mwdss
      breaksw
      
      case WEEK:
      case Week:
      case week:
         set t_type = week
      breaksw
            
      default:
      
      echo "SCRIPT ERROR: Setting T_TYPE = "$T_TYPE" is not recognized."
      echo "       Please check setting. Correct values are all, aveday,"
      echo "       mwdss, and week. Please check upper/lower case as well."
      echo " "
      set eflag = 1
   endsw
      
else
   echo "SCRIPT ERROR: environment variable T_TYPE is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1

endif

if ( $?SPINUP_DURATION ) then
else
   setenv SPINUP_DURATION 10    # Default to 10 days spinup duration
endif

# Make sure all executables are available
if ( -e $IOAPIDIR/juldate ) then
else
   echo "SCRIPT ERROR: Program juldate is needed but cannot be found in IOAPIDIR:"
   echo "      "$IOAPIDIR
   echo " "
   set eflag = 1   
endif

if ( -e $IOAPIDIR/gregdate ) then
else
   echo "SCRIPT ERROR: Program gregdate is needed but cannot be found in IOAPIDIR:"
   echo "      "$IOAPIDIR
   echo " "
   set eflag = 1   
endif

if ( -e $IOAPIDIR/datshift ) then
else
   echo "SCRIPT ERROR: Program datshift is needed but cannot be found in IOAPIDIR:"
   echo "      "$IOAPIDIR
   echo " "
   set eflag = 1   
endif

# Set name of output file & remove it if it's already there
if ( $?SMK_RUN_DATES ) then
   if ( -e $SMK_RUN_DATES ) then
      set fileuser = `/bin/ls -l $SMK_RUN_DATES | cut -d" " -f6`
      /bin/rm -rf $SMK_RUN_DATES
      if ( $status == 1 ) then
         echo "SCRIPT ERROR: Do not have permission to overwrite file:"
         echo "              $SMK_RUN_DATES"
         echo "              Contact user "$fileuser" to request permissions."
      endif
   endif

else
   echo "SCRIPT ERROR: Environment variable SMK_RUN_DATES is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1      
endif

# Set name of output file & remove it if it's already there
if ( $?PROCDATES ) then
   if ( -e $PROCDATES ) then
      set fileuser = `/bin/ls -l $PROCDATES | cut -d" " -f6`
      /bin/rm -rf $PROCDATES
      if ( $status == 1 ) then
         echo "SCRIPT ERROR: Do not have permission to overwrite file:"
         echo "              $PROCDATES"
         echo "              Contact user "$fileuser" to request permissions."
      endif
   endif
else
   echo "SCRIPT ERROR: Environment variable PROCDATES is a required input"
   echo "       setting to script set_days.csh, but it is not defined."
   echo " "
   set eflag = 1      
endif

if ( $eflag == 1 ) then
   exit ( 1 )
endif


# Set names of temporary files using the csh process number

set holiday_tmpfile = .holiday_tmpfile_$$
set holiday_intermd = .holiday_intermd_$$
set monset_tmpfile  = .monset_tmpfile_$$
set uniqmon_tmpfile = .uniqmon_tmpfile_$$
set outdates_tmpfile = .outdates_tmpfile_$$

# Define parameters
set days_per_month = (31 28 31 30 31 30 31 31 30 31 30 31)

# Define parameters for each temporal output resolution
switch ( $t_type )
   case all:
      set annual = y
   breaksw

   case aveday:
      set annual = n
      set nday_needed = 1
      set day_names = (Tuesday 0)
   breaksw

   case mwdss:
      set annual = n
      set nday_needed = 4
      set day_names = (Monday Tuesday Saturday Sunday)
   breaksw
   
   case week:
      set annual = n
      set nday_needed = 7
      set day_names = (Monday Tuesday Wednesday Thursday Friday Saturday Sunday)
   breaksw
   
endsw


# Determine if the current year is a leap year, and reset the days per
#   month accordingly
set ndays_in_yr = `$IOAPIDIR/juldate 12 31 $year | grep "," | cut -d, -f2 | cut -c6-8`
if ( $ndays_in_yr == 366 ) then
  set days_per_month[2] = 29
endif

# Create intermediate file of all Julian dates to be modeled as holidays
#   (includes dates for $year in HOLIDAYS file and day after these)
#   Do this regardless of whether holidays are being modeled separately, so that the
#      same dates will be obtained regardless for non-holiday dates.

/bin/rm -rf $holiday_tmpfile
/bin/rm -rf $holiday_intermd

grep $year $HOLIDAYS > $holiday_tmpfile
set nline = `cat $holiday_tmpfile | wc -l`
set n = 0
set first_rec = y
while ( $n < $nline )   # Loop through lines of holidays temporary file
   @ n = $n + 1
   set line = ( `head -$n $holiday_tmpfile | tail -1` )
   set h_mon = $line[2]
   set h_day = $line[3]
   set h_year = $line[4]

   # Check that year matches
   if ( $year != $h_year ) then
      echo "SCRIPT ERROR: year grepped from HOLIDAY file was $year"
      echo "       but year found in holiday temporary file was $h_year."
      echo "       Something strange is going on in set_days.csh."
      exit ( 1 )
   endif

   # Get julian date of this date and write to intermediate file
   set h_jdate = `$IOAPIDIR/juldate $h_mon $h_day $h_year | grep , | cut -d, -f2`
   set h_desc  = ( `$IOAPIDIR/gregdate $h_jdate | grep ,` )
   set h_greg  = `$IOAPIDIR/datshift $h_jdate 0`
   if ( $restart == Y && $h_jdate < $RESTART_JDATE ) then
   else
      if ( $first_rec == y ) then
         echo $h_jdate $h_greg H $h_desc > $holiday_intermd
         set first_rec = n
      else
         echo $h_jdate $h_greg H $h_desc >> $holiday_intermd
      endif
   endif

   # Unless this is the last day of the year, also output the next day
   if ( $h_mon == 12 && $h_day == 31 ) then
   else
      set nextday = `$IOAPIDIR/datshift $h_year$h_mon$h_day 1`
      set h_mon = `echo $nextday | cut -c5-6`
      set h_day = `echo $nextday | cut -c7-8`
      set h_jdate = `$IOAPIDIR/juldate $h_mon $h_day $h_year | grep , | cut -d, -f2`
      set h_desc  = ( `$IOAPIDIR/gregdate $h_jdate | grep ,` )
      set h_greg  = `$IOAPIDIR/datshift $h_jdate 0`

      # Don't output if date is not >= restart date
      if ( $restart == Y && $h_jdate < $RESTART_JDATE ) then
      else
         if ( $first_rec == y ) then
            echo $h_jdate $h_greg H $h_desc > $holiday_intermd
            set first_rec = n
         else
            echo $h_jdate $h_greg H $h_desc >> $holiday_intermd
         endif
      endif

   endif

end  # End loop through holidays temporary file

# Sort holidays intermediate file to ensure only unique entries
mv $holiday_intermd $holiday_tmpfile
sort -u $holiday_tmpfile > $holiday_intermd

# Count number of months (may be duplicates) listed in month array
set nm = 0
foreach m ( $monset )
   @ nm = $nm + 1
end

# Create list of months that need to be modeled & their spinup settings
echo $monset[1]\m $spinup_arr[1] > $monset_tmpfile
set m = 1
while ( $m < $nm )
   @ m = $m + 1
   echo $monset[$m]\m $spinup_arr[$m] >> $monset_tmpfile      
end

# Uniquely sort the months that need to be run & their spinup settings
sort -u $monset_tmpfile > $uniqmon_tmpfile

# If not modeling all days in the year
if ( $annual != y ) then

   # Determine the number of months that need to be run
   set nmon = `cat $uniqmon_tmpfile | wc -l`

   # Loop through unique months that need to be modeled
   set n = 0
   while ( $n < $nmon )
      @ n = $n + 1
      
      set done_month = n
      set mon = `head -$n $uniqmon_tmpfile | tail -1 | cut -d"m" -f1`
      set nset = 0
      set set_stat = 0

      # Initialize arrays for current month
      switch ( $t_type )
	 case aveday
	    set day_jdate = (0 0)
	    set day_count = (0 0)
	 breaksw

	 case mwdss
	    set day_jdate = (0 0 0 0)
	    set day_count = (0 0 0 0)
	 breaksw

	 case week
	    set day_jdate = (0 0 0 0 0 0 0)
	    set day_count = (0 0 0 0 0 0 0)
	 breaksw

      endsw

      
      # Iterate through attempts at setting all dates to non-holidays
      while ( $done_month == n )
      
         # Initialize holiday check	    
	 set all_check = okay
	    
         # Create array of days needed for this month...
         # Iterate through the number of days that need setting
	 set d = 0
	 while ( $d < $nday_needed )
	 
	    @ d = $d + 1
	    
	    # For initial setting
	    if ( $day_jdate[$d] == 0 ) then

               # Initialize day count for the first day needed in the list (after the
	       #    first day, the days will be obtained sequentially after that)
	       # Set to use the second day of the month, because we never want to select 
	       #    the first day, because of the time-zones / month transition.
	       if ( $d == 1 ) then
	          set day = 1
               endif
	       
               # Iterate through days in the month starting with the last day identified
	       while ( $day_jdate[$d] == 0 )
	          @ day = $day + 1

		  set f_name  = `$IOAPIDIR/juldate $mon $day $year | grep , | cut -d, -f1`
		  set f_jdate = `$IOAPIDIR/juldate $mon $day $year | grep , | cut -d, -f2`

		  if ( $f_name == $day_names[$d] ) then
		     set day_jdate[$d] = $f_jdate
		     set day_count[$d] = $day
		  endif

	       end
		
            # If value has been initialized, then reset to the next week, since this date must have
	    #   failed the holidays test
	    else
	    
	       # Reset only if the set status indicates the value is not final
	       if ( $set_stat == 0 ) then
	          set day = $day_count[$d]
	          @ day = $day + 7            # Check the next week
		  set day_jdate[$d] = `$IOAPIDIR/juldate $mon $day $year | grep , | cut -d, -f2`
		  set day_count[$d] = $day
	       endif
	    
	    endif
	 
            # Ensure that the date does not fall on a holiday or have a holiday 
	    #    preceding or following.  Ensure that none of the dates are on the
	    #    first day of the month.
	    # Initialize checking array
	    set check_jdate = (0 0 0)
	    
	    # Set first date
	    if ( $day_jdate[$d] > ${year}001 ) then
	       set jd = $day_jdate[$d]
	       @ jd = $jd - 1
	       set check_jdate[1] = $jd
	    endif
	    
	    # Set second date
	    set check_jdate[2] = $day_jdate[$d]
	    
	    # Set the third date
	    if ( $day_jdate[$d] < $year$ndays_in_yr ) then
	       set jd = $day_jdate[$d]
	       @ jd = $jd + 1
	       set check_jdate[1] = $jd
	    endif

            # Initialize check for this iteration
	    set check_this_iteration = okay

	    # Loop through three days to check
	    set i = 0
	    while ( $i < 3 && $check_this_iteration == okay )
	       @ i = $i + 1
	       if ( $check_jdate[$i] > 0 ) then
   	          grep -q "$check_jdate[$i] " $holiday_intermd
	       
		  # If match found
		  if ( $status == 0 ) then
	             set check_this_iteration = failed
        	  endif
	       endif
	    end	   

	    if ( $check_this_iteration == failed ) then
	       set all_check = failed
	    endif

         end   # end loop through number of dates needed
	 
         # If date is acceptable, then set flag to complete current month.     
	 if ( $set_stat == 0 && $all_check == okay ) then
	    set set_stat = 1
	    set done_month = y
	 endif
         
      end  # Loop through attempts in a month
      
      # Output dates needed to output file
      set d = 0
      while ( $d < $nday_needed )
	 @ d = $d + 1
	 
         set d_desc  = ( `$IOAPIDIR/gregdate $day_jdate[$d] | grep ,` )
         set d_greg  = `$IOAPIDIR/datshift $day_jdate[$d] 0`

         # Don't output if date is not >= restart date
         if ( $restart == Y && $day_jdate[$d] < $RESTART_JDATE ) then
         else
	    if ( $d == 1 && $n == 1 ) then
	       echo $day_jdate[$d] $d_greg"   "$d_desc > $outdates_tmpfile
	    else
	       echo $day_jdate[$d] $d_greg"   "$d_desc >> $outdates_tmpfile
	    endif
         endif
      end
	 
   end # Loop through months
   
   # Append all holiday dates that apply to the months being modeled to the 
   # end of the output and resort it to the final output file.
   if ( $RUN_HOLIDAYS == Y ) then
      cp $outdates_tmpfile $SMK_RUN_DATES
      set n = 0
      while ( $n < $nmon )
         @ n = $n + 1
         set mon = `head -$n $uniqmon_tmpfile | tail -1 | cut -d"m" -f1`
	 if ( $mon < 10 ) then
	    set mn2 = "0"$mon
	 else
	    set mn2 = $mon
	 endif
	 
	 grep " $year$mn2.. H" $holiday_intermd >> $SMK_RUN_DATES

      end
      
      mv $SMK_RUN_DATES $outdates_tmpfile
      sort $outdates_tmpfile > $SMK_RUN_DATES
      
   else
      sort $outdates_tmpfile > $SMK_RUN_DATES
   endif

# If modeling all days in the year (adjusted via spinup or not)
else

   # Loop through all days per year
   set d = 0
   set firstout = 0
   while ( $d < $ndays_in_yr )
      @ d = $d + 1
      
      if ( $d < 10 ) then
         set jdate = ${year}00$d
      else

	  if ( $d < 100 ) then
	     set jdate = ${year}0$d
	  else
             set jdate = $year$d
          endif
      endif
      
      set d_desc = ( `$IOAPIDIR/gregdate $jdate | grep ,` )
      set d_greg = `$IOAPIDIR/datshift $jdate 0`
      
      # Decide whether the month is in the month list
      set mon = `echo $d_greg | cut -c5-6`    
      @ mon = `expr $mon + 0`     ## needed to remove any leading 0 when mon < 10      
      grep -q ^${mon}m $uniqmon_tmpfile
      set mon_stat = $status
      set spinup_setting = ( `grep ^${mon}m $uniqmon_tmpfile | cut -d"m" -f2` )

      # Decide whether day of the month is correct or not
      set day_stat = 0
      set sdate = $jdate
      switch ( $spinup_setting )
	 case sp:
	    set day = `echo $d_greg | cut -c7-8`
	    set cutoff = 0
	    @ cutoff = $days_per_month[$mon] - $SPINUP_DURATION
	    if ( $day <= $cutoff ) then
	       set day_stat = 1 
	    endif
            # Change year of date
            if ( $day_stat == 0 ) then
               set lastyear = $year
               @ lastyear = $lastyear - 1
               set sdate = $jdate     # Reset date for searching holidays file

               # to get correct julian day number (and handle the leap year
               # case, need to convert to gregorian and then back to julian
               # for the previous year and get the new jdate for that.  
               set gmon  = `$IOAPIDIR/datshift $jdate 0 | cut -c5-6`
               set gday  = `$IOAPIDIR/datshift $jdate 0 | cut -c7-8`
               set jdate = `$IOAPIDIR/juldate $gmon $gday $lastyear | grep , | cut -d"," -f2`

               set d_desc = ( `$IOAPIDIR/gregdate $jdate | grep ,` )
               set d_greg = `$IOAPIDIR/datshift $jdate 0`
            endif
	 breaksw

         case s:
	    set day_tmp = `echo $d_greg | cut -c7-8`
            set day = `expr $day_tmp + 0`  # removes leading 0 for values < 10		
	    set cutoff = 0
	    @ cutoff = $days_per_month[$mon] - $SPINUP_DURATION
	    if ( $day <= $cutoff ) then
	       set day_stat = 1 
	    endif
	 breaksw

	 case p:
	    set day_tmp = `echo $d_greg | cut -c7-8`
            set day = `expr $day_tmp + 0`  # removes leading 0 for values < 10		
	    set cutoff = 0
	    @ cutoff = $days_per_month[$mon] - $SPINUP_DURATION
	    if ( $day > $cutoff ) then
	       set day_stat = 1 
	    endif	 
	 breaksw   
      endsw
      
      # If month found in list and day is of interest, then output the dates.
      if ( $mon_stat == 0 && $day_stat == 0 ) then
	 
	 # Learn whether the date is a holiday or not
	 grep -q "$sdate " $holiday_intermd

	 # If match found
	 if ( $status == 0 ) then
	    set d_holiday = H
	 else
            set d_holiday = " "
	 endif

         # Only output if restart is needed and date is after restart date
         if ( $restart == Y && $jdate < $RESTART_JDATE ) then
         else
	    if ( $firstout == 0 ) then
	       set firstout = 1
               echo $jdate $d_greg" $d_holiday "$d_desc > $SMK_RUN_DATES
	    else
               echo $jdate $d_greg" $d_holiday "$d_desc >> $SMK_RUN_DATES 
	    endif
         endif
      endif

   end

endif

# Create PROCDATES file from existing SMK_RUN_DATES file.
set nlines = `cat $SMK_RUN_DATES | wc -l`
set n = 0
while ( $n < $nlines )
   @ n = $n + 1
   set line = ( `head -$n $SMK_RUN_DATES | tail -1` )
   echo $line[2] $G_STTIME $G_RUNLEN >> $PROCDATES
end

chmod g+w $SMK_RUN_DATES $PROCDATES

# Remove all temporary files
/bin/rm -rf $holiday_tmpfile
/bin/rm -rf $holiday_intermd
/bin/rm -rf $monset_tmpfile
/bin/rm -rf $uniqmon_tmpfile
/bin/rm -rf $outdates_tmpfile
