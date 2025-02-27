#!/bin/tcsh -f
# Script for automatically removing or moving log file

if ( -e $TMPLOG ) then

   #
   ## Ensure automatic deletion variable is initialized
   #
   if ( ! $?AUTO_DELETE_LOG ) then
      set auto_del = 'N'

   else

      if ( $AUTO_DELETE_LOG == 'y' || $AUTO_DELETE_LOG == 'Y' ) then
         set auto_del = 'Y'
      else
         set auto_del = 'N'
      endif

   endif

   #
   ## Check to see if automatic deleting has been set, and if so, delete log
   ## file
   #
   if ( $auto_del == 'Y' ) then
      echo ' '
      echo 'SCRIPT NOTE: Automatically deleting log file.'
      echo '      '$TMPLOG
      echo ' '
      /bin/rm -rf $TMPLOG
   endif

   #
   ## If not removed, move log file to a different name
   #
   if ( -e $TMPLOG ) then

      set stat = 0
      set d = 0
      while ( $stat == 0 )

	 @ d = $d + 1
	 if ( ! -e ${TMPLOG}_$d ) then
             echo ' '
             echo '   Moving log file:'
             echo '      '$TMPLOG
             echo '   To file:'
             echo '      '${TMPLOG}_$d
             mv $TMPLOG ${TMPLOG}_$d
             set stat = 1
	 endif

      end

   endif

endif
  
exit (0)
