
C Version "@(#)$Id$ $Source$ $Date$ 

        LOGICAL FUNCTION GETPKT( NAME, CDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:  
C       Rewind CDEV, then read thru until line with control packet
C       NAME is encountered.
C
C  RETURNS:
C       TRUE if packet found; end-of-file; FALSE otherwise
C
C  PRECONDITIONS REQUIRED:
C       Control file correctly opened on CDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       M3EXIT
C
C  REVISION  HISTORY:
C       Prototype 8/95 by CJC for program CTLAMAT.
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*)   NAME    !  packet name for packet to be skipped
        INTEGER         CDEV    !  unit number for control file
        

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        CHARACTER*240   LINE    !  input line buffer
        INTEGER         IREC    !  current line number
        CHARACTER*80    MESG    !  error message buffer
        INTEGER         IOS     !  READ-status return


C***********************************************************************
C   begin body of subroutine  GETPKT

        REWIND( CDEV )

11      CONTINUE        !  head of loop reading lines in this packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading CONTROL FILE at line ',
     &              IREC, IOS
                CALL M3EXIT( 'GETPKT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .EQ. INDEX( LINE, NAME ) ) GO TO 11  ! to head of read-loop

        GETPKT = .TRUE.         !  if you get to here:  NAME found.
        RETURN
                

900     CONTINUE                !  exit from loop 11
        GETPKT = .FALSE.        !  NAME never found in file CDEV

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94910   FORMAT( A, I7,  '; IOSTAT=', I7 )

        END

