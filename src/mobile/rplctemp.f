
        SUBROUTINE RPLCTEMP( COUNTY, TEMPS, NCOUNTY, 
     &                       SCENARIO, NLINES, CTYPOS )

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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
        
        IMPLICIT NONE

C...........   INCLUDES:

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN)    :: COUNTY             ! current county being processed
        REAL,               INTENT (IN)    :: TEMPS( NCOUNTY, 24 ) ! hourly temps per county
        INTEGER,            INTENT (IN)    :: NCOUNTY  ! no. counties in temps array
        CHARACTER(LEN=150), INTENT (INOUT) :: SCENARIO( NLINES+1 ) ! scenario array
        INTEGER,            INTENT (IN)    :: NLINES             ! no. lines in scenario
        INTEGER,            INTENT (IN)    :: CTYPOS   ! position of county in TEMPS array

C...........   Local allocatable arrays

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        CHARACTER(LEN=150) CURRLINE           ! current line from scenario
        CHARACTER(LEN=150) RPLCLINE           ! replacement line
        CHARACTER(LEN=19)  COMMAND            ! Mobile 6 command        
        CHARACTER(LEN=300) MESG               !  message buffer

        CHARACTER*16 :: PROGNAME = 'RPLCTEMP'   ! program name
        
C***********************************************************************
C   begin body of subroutine RPLCTEMP
        
        DO I = 1, NLINES
        
            CURRLINE = SCENARIO( I )

C.............  Skip comment lines
            IF( CURRLINE( 1:1 ) == '*' ) CYCLE
            IF( CURRLINE( 1:1 ) == '>' ) CYCLE
            
C.............  Get Mobile6 command                   
            COMMAND = CURRLINE( 1:19 )            

C.............  Search command for min/max command
            IF( INDEX( COMMAND, 'MIN/MAX TEMPERATURE' ) > 0 ) THEN
            	
C.................  Move rest of scenario down one line to make room
C                   for second line of temperatures
                DO J = NLINES, I + 1, -1
                    SCENARIO( J + 1 ) = SCENARIO( J )
                END DO

            END IF

C.............  Replace either temperature command with hourly temperatures                	     
            IF( INDEX( COMMAND, 'TEMPERATURE' ) > 0 ) THEN

                WRITE( RPLCLINE,94020 ) 
     &               'HOURLY TEMPERATURES: ', TEMPS( CTYPOS,1:12 )
                SCENARIO( I ) = RPLCLINE
                
                WRITE( RPLCLINE,94030 ) TEMPS( CTYPOS,13:24 )
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
        