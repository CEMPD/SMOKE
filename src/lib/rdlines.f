
        SUBROUTINE RDLINES( FDEV, DESCRIPT, NLINES, LINES )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the lines of an ASCII file to an array of strings
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         TRIMLEN

        EXTERNAL CRLF, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        INTEGER       FDEV            !  file unit number
        CHARACTER*(*) DESCRIPT        !  file description
        INTEGER       NLINES          !  number of lines in file
        CHARACTER*(*) LINES( NLINES ) !  ASCII lines in file

C...........   Other local variables
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         L, LSAV, L2   !  length indices
        CHARACTER*300   LINE    !  line buffer
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDLINES' ! program name

C***********************************************************************
C   begin body of subroutine RDLINES

        L = LEN( LINES( 1 ) )     

        IREC = 0
        LSAV = 0
11      CONTINUE

            READ( FDEV, 93000, END=22, IOSTAT=IOS ) LINE
            IREC = IREC + 1
 
            L2 = TRIMLEN ( LINE )

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &                 'Error', IOS,  'reading ' // 
     &                 DESCRIPT( 1:TRIMLEN( DESCRIPT ) ) //
     &                 ' file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSEIF( L2 .GT. L ) THEN
                IF( L2 .GT. LSAV ) LSAV = L2

            ELSE

                IF( IREC .LE. NLINES ) THEN
                    LINES( IREC ) = LINE
                END IF

            ENDIF      !  if fip zero, or nn000, or not.

            GO TO  11

22      CONTINUE        !  exit from loop reading FDEV

        IF( IREC .GT. NLINES ) THEN
            WRITE( MESG,94010 ) 'WARNING: ' // 
     &             DESCRIPT( 1:TRIMLEN( DESCRIPT ) ) //
     &             CRLF() // BLANK10 // 'file only read for first ',
     &             NLINES, ' lines of ', IREC, ' total lines.'
            CALL M3MSG2( MESG ) 
        ENDIF

        IF( LSAV .GT. 0 ) THEN
            WRITE( MESG,94010 ) 'ERROR: ' // 
     &             DESCRIPT( 1:TRIMLEN( DESCRIPT ) ) //
     &             CRLF() // BLANK10 // 'file line width is ',
     &             LSAV, ' but allocated string length is ', L
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
        ENDIF

C.........  Rewind file

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
