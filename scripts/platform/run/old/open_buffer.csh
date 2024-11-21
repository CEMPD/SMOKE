#!/bin/csh -f
#
# open_buffer.csh - Created January 24, 2008 by Christopher T Allen
# This script finds a display buffer for the calling script to use when
# generating PAVE plots.
# The buffer to be used is returned to the calling script as the exit status.
# Script takes one command line argument: the process ID of the calling script

if ( $#argv != 1 ) then

   echo "***********************************************************"
   echo "SCRIPT ERROR: Script requires one command line argument"
   echo "Example: > open_buffer.csh [PROCESS ID]"
   echo "***********************************************************"
   exit (1)
   
endif

set calling_pid = $argv[1]

## Valid buffers: 3 through 9, and 11
set buffers = ( 3 4 5 6 7 8 9 11 )

## out_buff will be sent back to the calling script; initialize to zero in case
#  no available buffers are found
set out_buff = 0

## This script uses a file called "used_buffers.txt". When a buffer is "claimed",
#  the process ID and buffer number are added to this file.
#  When the script finishes using the buffer, the entry is removed from used_buffers.txt,
#  and the buffer is closed. The buffer will NOT be closed if the buffer was
#  already opened by another user (but is still available).
setenv BUFFER_FILE $SCRIPTS/run/used_buffers.txt

## Make sure BUFFER_FILE exists
if (! -e $BUFFER_FILE) then

  echo "*************************************************"
  echo "ERROR: $BUFFER_FILE does not exist!"
  echo "*************************************************"
  exit (1) # Even though we're using the exit status to return the buffer number,
           # 1 isn't a valid buffer, so we can still use it to signal an error
	   
endif
   
## Loop over each valid buffer
foreach buff ($buffers)

   ## Search BUFFER_FILE for buffer number
   if (`grep "buff:$buff" $BUFFER_FILE | wc -l` > 0) then # Buffer is in BUFFER_FILE
   
      ## Use "ps" to see if the process ID indicated in BUFFER_FILE is actually running
      ## If process is running, then try the next buffer
      set claiming_pid = `grep "buff:$buff" $BUFFER_FILE | cut -d"," -f1`
     
      if (`ps -ef | grep $claiming_pid | wc -l` > 1) then
     
         echo "Buffer $buff is in use; trying next buffer"
         continue
       
      else # Buffer can be used
     
         echo "The process claiming buffer $buff is not active"
       
         ## Remove "dormant" process from BUFFER_FILE
         sed "/$claiming_pid,buff:$buff/d" $BUFFER_FILE > ./x
         mv ./x $BUFFER_FILE
	 
      endif
   endif
    
   ## At this point, the buffer is free to use
       
   ## Add current process to BUFFER_FILE
   echo "$calling_pid,buff:$buff" >> $BUFFER_FILE
   set out_buff = $buff
	 
   ## Use "ps" to see if the buffer is already open. If it is not open,
   ## start it. Either way, use it.
   if (`ps -ef | grep "Xvfb :$buff" | wc -l` == 1) then
      echo "Starting and using buffer $buff"
      Xvfb :$buff -screen 0 1280x1024x24 -ac -fp /usr/X11R6/lib/X11/fonts/misc &
   else
      echo "Using buffer $buff (buffer is already open)"
   endif
	
   break # Leave foreach loop

end #foreach 

exit ($out_buff)
