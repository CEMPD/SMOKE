#!/bin/csh -f

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
          set add_spin_up = Y
          breaksw
       case N:
       case n:
       case Y:
       case y:
          echo "SCRIPT ERROR: set_months.csh no longer uses SPINUP value of Y or N"
          echo "              SPINUP value needs to be 0 or a number between 1 and 20 days"
          exit ( 1 )
          breaksw
       default:
          echo "SCRIPT ERROR: Invalid SPINUP value "$argv[3]" when calling set_monhts.csh"
          echo "              SPINUP value needs to be 0 or a number between 1 and 20 days"
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
