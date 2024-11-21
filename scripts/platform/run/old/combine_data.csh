#!/bin/csh -f

# This script creates .lst files for inventory or concatenates 
#     datasets into a single dataset.

# Creation date: 20 jun 2007
# Original author: M. Houyoux

# Purpose: to provide a scriptable way of combining datasets from the EMF,
#     which will be stored in the EMF as separate "inputs" to a case.

# Calling syntax:
#    combine_data.csh [prefix] [outfile] [type] [Y/N monthly]
#
#    Where,
#       prefix = e.v. prefix to search for
#       outfile = e.v. for file to create
#       type = list | concat. Default is "list".
#       Y/N month = Y if e.vs. may use month numbers to denote month-specific
#                   files and N otherwise. N is default if not specified.
#
# Updated October 19, 2007 by M. Houyoux to remove blank lines from inputs
#       when in type=concat

# Check the number of arguments. It must have at least two arguments.
switch ( $#argv )
    case 1:
    case 0:
      echo "SCRIPT ERROR: Script requires at least two arguments."
      echo "  The syntax for calling this script is:"
      echo " "
      echo "     combine_data.csh [prefix] [outfile] [type] [month]"
      echo " "
      echo "         where,"
      echo "            [prefix]  = environment variable prefix to search for"
      echo "            [outfile] = environment variable for file to create"
      echo "            ["type"]    "= \"list\" or \"concat\". Default is \"list\".
      echo "                        This argument is optional."
      echo "            [month]   = the month number of interest for the file."
      echo "                        Use values 01 through 12. Default is not"
      echo "                        month-specific. This argument is optional."
      exit( 1 )

endsw

set scriptname = combine_data.csh

# Create formatted date
set date = `date +%m/%d/%Y`

# Parse out arguments
set valid_month = N
if ( $#argv > 3 ) then
    set month = $argv[4]
    foreach m ( 01 02 03 04 05 06 07 08 09 10 11 12)
       if ( $m == $month ) then
          set valid_month = Y
       endif
    end
    if ( $valid_month == N ) then
       echo SCRIPT ERROR: month argument is invalid value of \"$month\"
       echo "              "Enter valid value of 01 through 12.
       exit ( 1 )
    endif
else
    set month = N
endif

if ( $#argv > 2 ) then
   switch ( $argv[3] )
      case list:
      case concat:
         set type = $argv[3]
         breaksw
      default:
         echo SCRIPT ERROR: Script setting for \"type\" is set to
         echo "              "invalid value of \"$argv[3]\". Valid values
         echo "              "are \"list\" and \"concat\", and these are
         echo "              "case sensitive.
   endsw

endif

set prefix = $argv[1]
set outfile = $argv[2]

# Message      
echo Creating dataset $outfile 
echo "     using script $scriptname"

# Inializations
set eflag = N   # Initialize 
set monname = (jan feb mar apr may jun jul aug sep oct nov dec)
      
# Get list of input environment variables
set envlist = ( `env | grep $prefix | cut -d\= -f1` )
set filelist = ( `env | grep $prefix | cut -d\= -f2` )

# If e.v. MULTIFILE_NAMEBREAK is set > 0, this indicates that EMF external multifiles are in use,
# and therefore the prefix in combination with the _MULTI_ field should be treated specially
# The value of this variable is set to the number of underscores to use to set 
#        prefix of file name for the files in "list" mode.
set multifile_break = 0
if ( $?MULTIFILE_NAMEBREAK ) then
   if ( $MULTIFILE_NAMEBREAK > 0 ) then
      set multifile_break = $MULTIFILE_NAMEBREAK
   endif
endif

# Check if list is empty before proceeding
if ( "$envlist" == "" ) then
   echo "SCRIPT ERROR: No environment variables with prefix $prefix "
   echo "              have been defined."
   exit ( 1 )
endif

# Remove output file
if ( -e $outfile ) then
   /bin/rm -f $outfile
endif

# Temporary list printout
echo Processing environment variables $envlist

# Initializations prior to loop
set n = 0
set firstime = Y
set nmax = 90
set use_ev_array = ( N N N N N N N N N N N N N N N N N N N N \
                     N N N N N N N N N N N N N N N N N N N N \
                     N N N N N N N N N N N N N N N N N N N N \
                     N N N N N N N N N N N N N N N N N N N N \
                     N N N N N N N N N N N N N N N N N N N N \
                     N N N N N N N N N N )

# Loop through environment variables and files to create output file
#    For "list" mode, this loop creates the actual .lst file.
#    For "concat" mode, this loop creates only the header of the file.
foreach ev ( $envlist )

   @ n = $n + 1

   set use_ev = Y        # Use this e.v. in creating the output file (initialize)
   set ev_month = `echo $ev | cut -d_ -f2 | cut -c1-2`   # If present, the value of the month in the e.v.
   set ev_month_yn = N   #  No month in environment variable

   # Determine if contains _MULTI_ in name
   set multi = N
   echo $ev | grep -q _MULTI_ 
   if ( $status == 0 ) then
      set multi = Y
   endif

   #  See whether environment variable contains a month or not
   foreach m ( 01 02 03 04 05 06 07 08 09 10 11 12 )
      if ( $m == $ev_month ) then
         set ev_month_yn = Y
      endif
   end

   # If not using monthly data, then skip any environment variables
   #     with the monthly denotation, unless it has _MULTI_ when namebreak > 0
   if ( $month == N && $ev_month_yn == Y ) then
      if ( $multifile_break > 0 && $multi == Y ) then
      else 
         set use_ev = N
      endif
   endif

   # If using monthly data, then skip any environment variables specified
   #     for a different month.
   if ( $valid_month == Y && $ev_month_yn == Y ) then
      if ( $month != $ev_month ) then
         set use_ev = N
      endif
   endif

   set eflag = N   # Initialize before next section
 
   # If environment variable in list is to be used, the process it, otherwise skip
   if ( $use_ev == Y ) then
     
      # For multifile case, use special approach to set filename, otherwise, just use filelist directly
      if ( $multifile_break > 0 && $multi == Y ) then
         set nameprefix = ( `echo $filelist[$n] | cut -d_ -f1-$multifile_break` )
         echo nameprefix = $nameprefix
         echo month = $month
         echo monname = $monname[$month]
         #note: asterisk in the middle of the file name is to handle CEM file 
         #      naming convention, which includes the date range of the files,
         #      since some files are only for the last day of the month.
         set filename = ( `ls ${nameprefix}_*${monname[$month]}_*` )
      else
         set filename = $filelist[$n]
      endif

      # If data file exists then process it
      if ( -e $filename ) then

         # Store use_ev in array for later use (in another loop)     
         set use_ev_array[$n] = $use_ev

         switch( $type )

            # Create a .lst file (for input to Smkinven)
            case list:

               # insert #LIST header into the output file.
               if ( $firstime == Y ) then
                  set firstime = N
                  echo "#LIST" > $outfile
                  #echo "#DESC This file created by script $scriptname on $date" >> $outfile
               endif

               # insert each path and filename into the output file
               echo $filename >> $outfile

               breaksw

            # In this loop, process only the header lines and output
            #    these to the new file's header.  Will loop again to
            #    actually output the non-header data to the output file.
            case concat:

               # Check that nmax has not been violated
               if ( $n > $nmax ) then
                  echo "SCRIPT ERROR: $scriptname script cannot handle more than"
                  echo "              $nmax files when type=concat is used."
                  set eflag = Y

               # Otherwise, process header
               else
                  echo "#### Header from data file:" >> $outfile
                  echo "####     "$filename >> $outfile
                  grep "^#" $filename | grep -v ^\$ >> $outfile
               endif

               breaksw
            default:
               echo SCRIPT ERROR: Script setting for \"type\" is set to
         endsw

      # If data file doesn't exist, then 
      else
         echo "SCRIPT ERROR: File $ev"
         echo "              does not exist. Will not add to output file."
         set eflag = Y
      endif

   endif    # End if e.v. is used or not

end

# Abort if an error was found in previous loop
if ( $eflag == Y ) then
   /bin/rm -f $outfile
   exit ( 1 )
endif

# For concatenate mode, loop through variables/files again to concatenate contents
#    of data files into outfile (instead of just the headers done in the previous loop).

if ( $type == concat ) then

   set n = 0
   foreach ev ( $envlist )

      @ n = $n + 1
      if ( $use_ev_array[$n] == Y ) then
         echo Including file: $filelist[$n]
         grep -v "^#" $filelist[$n] | grep -v ^\$ >> $outfile
      endif  

   end

endif  #  concat mode 

# Successful script end
exit ( 0 )
