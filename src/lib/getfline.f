
        INTEGER FUNCTION GETFLINE( IDEV, DESCRIPT )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION: 
C     Counts the number of lines in an ASCII file
C
C  PRECONDITIONS REQUIRED:
C     File opened and unit number provided
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       prototype 10/98 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

        IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        INTEGER       IDEV         ! Unit number for ASCII file
        CHARACTER*(*) DESCRIPT     ! Description of file

C...........   EXTERNAL FUNCTIONS:
        INTEGER    TRIMLEN

        EXTERNAL   TRIMLEN

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       L1           ! length of DESCRIPT
        INTEGER       ICNT         ! line counter
        INTEGER       IOS          ! i/o status

        CHARACTER*1   BUFFER       ! ASCII LINE from X-ref file
        CHARACTER*256 MESG         ! Message buffer

        CHARACTER*16 :: PROGNAME = 'GETFLINE' ! program name

C***********************************************************************
C   begin body of subroutine GETFLINE

        L1 = TRIMLEN( DESCRIPT )

        ICNT = 0

C.........  Loop through lines of file, counting the lines
11      CONTINUE

            READ( IDEV, 93000, IOSTAT=IOS, END=22 ) BUFFER
 
            ICNT = ICNT + 1
 
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG,94010 )
     &                'I/O error', IOS,
     &                'scanning ' // DESCRIPT( 1:L1 ) // 
     &                ' file at line', ICNT
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        GO TO 11

22      CONTINUE

        GETFLINE = ICNT  

        REWIND( IDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END
