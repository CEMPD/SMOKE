
        SUBROUTINE RPLCTEMP( CTYPOS, NSCEN, NLINES, SCENARIO, STSCEN )

C***********************************************************************
C  subroutine body starts at line 75
C
C  DESCRIPTION:
C       Replaces temperatures in the MOBILE6 input file with values
C       from Premobl
C
C  PRECONDITIONS REQUIRED:
C       Temperature command must be present in the MOBILE6 scenario
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
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
        
C.........  MODULES for public variables
C.........  This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR
        
        IMPLICIT NONE

C...........   INCLUDES:

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)  CRLF    
        
        EXTERNAL  CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,        INTENT (IN)    :: CTYPOS             ! position of county in TKHOUR
        INTEGER,        INTENT (INOUT) :: NSCEN              ! no. actual lines in scenario
        INTEGER,        INTENT (IN)    :: NLINES             ! no. lines in scenario array
        CHARACTER(150), INTENT (INOUT) :: SCENARIO( NLINES ) ! scenario array
        INTEGER,        INTENT (INOUT) :: STSCEN             ! start of scenario level commands

C...........   Local allocatable arrays

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        CHARACTER(150) CURRLINE           ! current line from scenario
        CHARACTER(150) RPLCLINE           ! replacement line
        CHARACTER(19)  COMMAND            ! Mobile 6 command        
        CHARACTER(300) MESG               !  message buffer

        CHARACTER(16) :: PROGNAME = 'RPLCTEMP'   ! program name
        
C***********************************************************************
C   begin body of subroutine RPLCTEMP
        
        DO I = 1, NSCEN
        
            CURRLINE = SCENARIO( I )

C.............  Skip comment and blank lines
            IF( CURRLINE( 1:1 ) == '*' ) CYCLE
            IF( CURRLINE( 1:1 ) == '>' ) CYCLE
            IF( CURRLINE == ' ' ) CYCLE
            
C.............  Get Mobile6 command                   
            COMMAND = CURRLINE( 1:19 )            

C.............  Search command for min/max command
            IF( INDEX( COMMAND, 'MIN/MAX TEMPERATURE' ) > 0 ) THEN
            	
C.................  Move rest of scenario down one line to make room
C                   for second line of temperatures
                DO J = NSCEN, I + 1, -1
                    SCENARIO( J + 1 ) = SCENARIO( J )
                END DO

                NSCEN = NSCEN + 1

C.................  Update scenario level start if needed
                IF( I < STSCEN ) STSCEN = STSCEN + 1

            END IF

C.............  Replace either temperature command with hourly temperatures                	     
            IF( INDEX( COMMAND, 'TEMPERATURE' ) > 0 ) THEN

                WRITE( RPLCLINE,94020 ) 
     &               'HOURLY TEMPERATURES: ', TKHOUR( CTYPOS,1:12 )
                SCENARIO( I ) = RPLCLINE
                
                WRITE( RPLCLINE,94030 ) TKHOUR( CTYPOS,13:24 )
                SCENARIO( I + 1 ) = RPLCLINE
                
                EXIT

            END IF

        END DO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( A21, 12( 1X, F7.3 ) )
94030   FORMAT( 12( F7.3, :, 1X ) )
        
        END SUBROUTINE RPLCTEMP
        
