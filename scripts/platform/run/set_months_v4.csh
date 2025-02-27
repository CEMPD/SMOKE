#!/bin/tcsh -f

# This script sets the months needed for a given script run and
# determines which months have spinup
#
# This script should be "sourced" (not run) so that it will set the 
#  "MONTHS_LIST" and "SPINUP_LIST".
#
# spinup codes:
#     a - run all days of this month
#     s - spinup (run last N day of this month)
#     sp - spinup previous year (run last N day of this month)
#     p  - partial month  (run first M - N days of this month)
#      where N is the number of spinup days
#            M is the number of days in this month
#
# SPINUP_MONTH_END determines if script will set partial months
#                Y - no partial months
#                N - partial months (default behaviour)
#
#  Updated:
#    March 2009  A. Zubrow, IE at UNC, added option to run without 
#                                      partial months in spinup mode
#    29 Aug 2012, C. Allen, CSC: Added support for SPINUP_DURATIONs of up to 28 days.
#      Old limit was 20 days, which seemed awfully arbitrary.
#    9 Sep 2024: Let's try 31 days of spinup so that we can seamlessly do a full month.
#########################################################################

set option  = $argv[1]
set add_spin_up = N
if ( $# == 3 ) then
    switch ( $argv[3] )
       case 0: 
          set add_spin_up = N
          breaksw
       case 1:
       case 2:
       case 3:
       case 4:
       case 5:
       case 6:
       case 7:
       case 8:
       case 9:
       case 10:
       case 11:
       case 12:
       case 13:
       case 14:
       case 15:
       case 16:
       case 17:
       case 18:
       case 19:
       case 20:
       case 21:
       case 22:
       case 23:
       case 24:
       case 25:
       case 26:
       case 27:
       case 28:
       case 29:
       case 30:
       case 31:
          set add_spin_up = Y
          breaksw
       case N:
       case n:
       case Y:
       case y:
          echo "SCRIPT ERROR: set_months.csh no longer uses SPINUP value of Y or N"
          echo "              SPINUP value needs to be 0 or a number between 1 and 31 days"
          exit ( 1 )
          breaksw
       default:
          echo "SCRIPT ERROR: Invalid SPINUP value "$argv[3]" when calling set_months.csh"
          echo "              SPINUP value needs to be 0 or a number between 1 and 31 days"
          exit ( 1 )
    endsw
else
    echo "SCRIPT NOTE: set_months.csh defaulting to not using spinup days"
endif

switch ($option)
   case -m:
      set months = ( $argv[2] )
      set qtr_yn = 'N'
      breaksw
   case -q:
      set quarters = ( $argv[2] )
      set qtr_yn = 'Y'
      breaksw
   default:
      echo "SCRIPT ERROR: set_months.csh script requires -m or -q option, but $argv[1] given instead."
      echo "     This script expects to be called using one of the following options:"
      echo "        -m <monthlist>"
      echo "        -q <quarters>"
      echo "     Examples:"
      echo "         <script name> -m '1 2 3' : runs script for Jan, Feb, & Mar"
      echo "         <script name> -q 2       : runs script for the 2nd quarter,"
      echo "                                    including spin-up days specified"
      echo "                                    SPINUP_DURATION environment variable."
      exit( 1 )
endsw

## check if running full months, no partial months parameter is set
if ( ! $?SPINUP_MONTH_END) then
    setenv SPINUP_MONTH_END N
endif

## Running with no partial months for spinup
if ( $SPINUP_MONTH_END == "Y" ) then
    echo "SCRIPT NOTE: Setting months with no partial months"

    ## using quarters
    if ( $qtr_yn == Y ) then
	## using quarters
	set months = " "
	foreach q ( $quarters )
	    switch ( $q )
	     case 0: # for running annual inventory through SMOKE for reporting purposes only (C.Allen added 14 Feb 2019)
	         set months = ( $months 0 ) 
		 breaksw
	     case 1:
                 set months = ( $months 1 2 3 )
                 breaksw
             case 2:
                 set months = ( $months 4 5 6 )
                 breaksw
             case 3:
                 set months = ( $months 7 8 9 )
                 breaksw
             case 4:
                 set months = ( $months 10 11 12 )
                 breaksw
	    endsw
	end
    endif

    ## if spinup, prepend previous month
    if ( $add_spin_up == Y ) then
	## make months an array and calculate previous month
	set months = ( $months )
	@ pm = $months[1] - 1

	## if previous month = 0, reset to 12
	if ( $pm == 0 ) set pm = 12
	set months = ( $pm $months ) 
    endif

    setenv MONTHS_LIST "$months"

    ## initialize spinup array to all
    set spinup = " "
    foreach m  ($months)
	set spinup = ( $spinup a )
    end
    
    ## if spinup, change 1st element to spinup
    if ( $add_spin_up == Y ) then
	set spinup[1] = s
	## if spinup month is previous year, set code
	if ( $months[1] == 12 ) set spinup[1] = sp
    endif

    setenv SPINUP_LIST "$spinup"

    
    exit ( 0 )
endif  ## end of run without partial months

## Running with partial months for spinup (default approach)
## Set up months to be run and spin-up status
if ( $qtr_yn == "N" ) then

   set cnt = 0   
   set monthsnew = " "
   set spinup = " "
   foreach m ( $months )
      if ( $add_spin_up == Y ) then
	 switch ( $m )
	    case 1:
	       set monthsnew = ( $monthsnew 12 $m )
	       set spinup    = ( $spinup    sp  a )
               breaksw
            case 3:
            case 6:
            case 9:
               set monthsnew = ( $monthsnew $m )
               set spinup = ( $spinup p )
               breaksw
            case 4: 
            case 7: 
            case 10: 
               set pm = $m
               @ pm = $pm - 1
               set monthsnew = ( $monthsnew $pm $m )
               set spinup    = ( $spinup    s   a  )
               breaksw
            default:
               set monthsnew = ( $monthsnew $m )
               set spinup    = ( $spinup a )
         endsw

      else 
         set monthsnew = ( $monthsnew $m )
         set spinup    = ( $spinup a )

      endif
   end        #  end loop through months
   set months = ( $monthsnew )

## Set up months and spin-up status based on quarters
else
   set months = " "
   set spinup = " "
   foreach q ( $quarters )
   
      if ( $add_spin_up == Y ) then
         switch ( $q )
             case 1:
                 set months = ( $months 12 1 2 3 )
                 set spinup = ( $spinup sp a a p )
                 breaksw
             case 2:
                 set months = ( $months 3 4 5 6 )
                 set spinup = ( $spinup s a a p )
                 breaksw
             case 3:
                 set months = ( $months 6 7 8 9 )
                 set spinup = ( $spinup s a a p )
                 breaksw
             case 4:
                 set months = ( $months 9 10 11 12 )
                 set spinup = ( $spinup s a  a  a  )
                 breaksw
         endsw
      else
         switch ( $q )
	     case 0: # for running annual inventory through SMOKE for reporting purposes only (C.Allen added 14 Feb 2019)
	         set months = ( $months 0 ) 
                 set spinup = ( $spinup a )
		 breaksw
             case 1:
                 set months = ( $months 1 2 3 )
                 set spinup = ( $spinup a a p )
                 breaksw
             case 2:
                 set months = ( $months 4 5 6 )
                 set spinup = ( $spinup a a p )
                 breaksw
             case 3:
                 set months = ( $months 7 8 9 )
                 set spinup = ( $spinup a a p )
                 breaksw
             case 4:
                 set months = ( $months 10 11 12 )
                 set spinup = ( $spinup a  a  a  )
                 breaksw
         endsw
      endif
   end
endif


setenv MONTHS_LIST "$months"
setenv SPINUP_LIST "$spinup"

exit ( 0 )
