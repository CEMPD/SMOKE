
        SUBROUTINE OPENREPOUT( FILNAM, FDEV, MDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C    The OPENREPOUT routine open output file FILNAM, which can be a physical
C    or logical file name.  Returns FDEV to calling program.
C
C  PRECONDITIONS REQUIRED:
C    FILNAM is defined as a physical or logical file name.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Revised 7/2003 by A. Holland
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C...........   MODULES for public variables

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: RPT_

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)  CRLF
        INTEGER       GETEFILE
        INTEGER       JUNIT
        LOGICAL       SETENVVAR

        EXTERNAL    CRLF, GETEFILE, JUNIT, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FILNAM   ! physical or logical file name
        INTEGER     , INTENT(OUT) :: FDEV( RPT_%NUMFILES ) ! unit number of file
        INTEGER     , INTENT(OUT) :: MDEV                  ! unit number of file for src mapping

C...........   Local variables
        INTEGER         I, L, L2        ! counters and indices
        INTEGER         IOS             ! i/o status

        CHARACTER(3)    FILENO          ! tmp file number buffer
        CHARACTER(10)   FMT             ! tmp format buffer
        CHARACTER(16)   VARBUF          ! logical file name buffer
        CHARACTER(300)  MESG            ! message buffer
        CHARACTER(300)  PNAME           ! physical file name on ENVSTR test
        CHARACTER(300)  NNAME           ! new physical file name buffer

        CHARACTER(16) :: PROGNAME = 'OPENREPOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENREPOUT

        IOS = -1
C.........  If file name is less than 16 characters, check if file name is a
C           defined environment variable
        L = LEN_TRIM( FILNAM )
        IF( L .LE. 16 ) THEN
            MESG = 'Check of file name for logical status'
            VARBUF = FILNAM
            CALL ENVSTR( VARBUF, MESG, ' ', PNAME, IOS )
        ELSE
            MESG = 'Logical file name "'//TRIM(FILNAM)//
     &             '" is exceeding max 16 characters'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Open output file(s)
        DO I = 1, RPT_%NUMFILES

C.............  If it is a defined environment variable, open as logical file name
            IF( IOS .EQ. 0 ) THEN

C.................  If muliple files are being output then open them
                IF( RPT_%RPTMODE .EQ. 1 ) THEN

C.....................  Create new file name
                    IF( I .GE. 10 ) THEN
                        FMT = '( A, I2 )'
                    ELSE IF( I .GE. 100 ) THEN
                        FMT = '( A, I3 )'
                    ELSE
                        FMT = '( A, I1 )'
                    END IF

                    L = LEN_TRIM( PNAME )
                    WRITE( NNAME, FMT ) PNAME( 1:L )//'_', I

C.....................  Set logical file name to new file name
                    IF( .NOT. SETENVVAR( FILNAM, NNAME ) ) THEN
                        MESG = 'Could not set logical file '//
     &                         'name for file ' // TRIM( NNAME )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

                FDEV( I ) = GETEFILE( FILNAM, .FALSE.,
     &                                    .TRUE., PROGNAME )

C.............  If it is not a defined environment variable, open as physical filename
            ELSE

                IF( RPT_%RPTMODE .EQ. 1 ) THEN

                    WRITE( FILENO, '( I3 )' ) I
                    FILENO = ADJUSTL( FILENO )
                    NNAME = FILNAM // '_' // FILENO

                END IF

                FDEV( I ) = JUNIT()
                OPEN( FDEV( I ),ERR=1006,FILE=NNAME,
     &                    STATUS='UNKNOWN',RECL=2500 )

            END IF

        END DO

C.........  Open source mapping ancillary output file
        IF( RPT_%SRCMAP ) THEN

            NNAME = TRIM( PNAME )//'_src_crosswalk.txt'
            IF( .NOT. SETENVVAR( 'SRC_CROSSWALK', NNAME ) ) THEN
                 MESG = 'Could not set logical file '//
     &                  'name for file ' // TRIM( NNAME )
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            MDEV = GETEFILE( 'SRC_CROSSWALK', .FALSE.,
     &                            .TRUE., PROGNAME )

        END IF

        RETURN

C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) ' Could not open file:' //
     &         CRLF() // BLANK5 // FILNAM( 1:LEN_TRIM( FILNAM ) )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE OPENREPOUT
        
 
