
        SUBROUTINE INITEM( LDEV, NAMEVERS, INPROGNM )

***********************************************************************
C  program body starts at line 81
C
C  DESCRIPTION:
C       Writes out abridged copyright information, calling program version,
C       web address for documentation, general program and info, program-
C       specific info, and prompts to continue running the calling program
C
C  PRECONDITIONS REQUIRED:
C       Unit number LDEV defined
C       Version name set in calling program using SCCS
C       Calling program name defined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines, PROGDESC
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C       Copied from prototype version 1.2  in 10/98 by M Houyoux 
C
C***********************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C.........  Subroutine arguments
        INTEGER       LDEV      ! Log file unit number
        CHARACTER*50  NAMEVERS  ! ASCII field e/ version number appended to end
        CHARACTER*16  INPROGNM  ! Calling program name

C.........  Parameters
        INTEGER       STDOUT
        INTEGER       YEAR
        PARAMETER   ( STDOUT = 6,
     &                  YEAR = 1999 )

C.........  External functions
        LOGICAL       GETYN
        INTEGER       TRIMLEN

        EXTERNAL      GETYN, TRIMLEN

C.........  Local variables
        REAL          VERSION

        INTEGER       IOUT( 2 )  ! output unit numbers for stdout and logfile
        INTEGER       I, J
        INTEGER       NLOOP      ! Number of times to loop through output
         
        CHARACTER*16  VERCHAR
        CHARACTER*300 LINE0, LINE1, LINE2, LINE3, LINE4, LINE5

        CHARACTER*16 :: PROGNAME = 'INITEM'  ! program name

C***********************************************************************
C   begin body of program INITEM

        LINE0 = 'SMOKE ---------------' 
        WRITE( LINE1,94020 ) 'Copyright (c)', YEAR, 
     &                      'MCNC--North Carolina Supercomputing Center'

        LINE2 = 'All rights reserved'
        LINE3 = 'See file COPYRIGHT for conditions of use.'

        LINE4 = 'Online documentation available at:' 
        LINE5 = '    http://www.envpro.mcnc.org/products/smoke'

C.........  Set up program version information
        VERCHAR = ADJUSTL( NAMEVERS( 1:TRIMLEN( NAMEVERS ) ) )

        IF( VERCHAR( 1:1 ) .EQ. '%' ) THEN
            VERCHAR = 'Dev'
        ELSE
            J = INDEX( VERCHAR, ' ' )
            VERCHAR = ADJUSTL( VERCHAR( J+1:TRIMLEN( VERCHAR ) ) )
            J = INDEX( VERCHAR, '$' )
            IF( J .GT. 1 ) VERCHAR = VERCHAR( 1:J-1 )
        ENDIF

C.........  Set up writing loop
        IOUT( 1 ) = STDOUT
        IOUT( 2 ) = LDEV
        IF( LDEV .NE. STDOUT ) THEN 
            NLOOP = 2
        ELSE
            NLOOP = 1
        ENDIF

        DO I = 1, NLOOP

            LDEV = IOUT( I )

C.............  Write copyright information

            WRITE( LDEV,92000 ) LINE0( 1:TRIMLEN( LINE0 ) )
            WRITE( LDEV,92000 ) LINE1( 1:TRIMLEN( LINE1 ) )
            WRITE( LDEV,92000 ) LINE2( 1:TRIMLEN( LINE2 ) )
            WRITE( LDEV,92000 ) 
            WRITE( LDEV,92000 ) LINE3( 1:TRIMLEN( LINE3 ) )
            WRITE( LDEV,92000 ) 

C.............  Write program version information

            WRITE( LDEV,92010 ) INPROGNM( 1:TRIMLEN( INPROGNM ) ), 
     &                          VERCHAR
           
C.............  Write web site information
            WRITE( LDEV,92000 ) LINE4( 1:TRIMLEN( LINE0 ) )
            WRITE( LDEV,92000 ) LINE5( 1:TRIMLEN( LINE1 ) )
           
C.............  Write program-specific information
            CALL PROGDESC( LDEV, INPROGNM )

C.............  Write general information for SMOKE programs

            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program start,',
     &  'using "setenv <logicalname> <pathname>").',
     &      ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &      ' '

        ENDDO  ! End of write loop

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( INPROGNM, 0, 0, 'Ending program.', 2 )
        END IF

        RETURN 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )

92010   FORMAT( 5X, 'Program ', A, ', Version ', A )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A )

94020   FORMAT( 10( A, :, I4, :, 1X )  )

94030   FORMAT( F5.1 )

        END
