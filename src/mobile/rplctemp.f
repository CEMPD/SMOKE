
        SUBROUTINE RPLCTEMP( COUNTY, TEMPS, NCOUNTY, 
     &                       SCENARIO, NLINES, CTYPOS )
        
        IMPLICIT NONE

C...........   INCLUDES:

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN)    :: COUNTY             ! current county being processed
        REAL,               INTENT (IN)    :: TEMPS( NCOUNTY, 24 ) ! hourly temps per county
        INTEGER,            INTENT (IN)    :: NCOUNTY  ! no. counties in temps array
        CHARACTER(LEN=150), INTENT (INOUT) :: SCENARIO( NLINES ) ! scenario array
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

C.............  Search command for temperature commands
            IF( INDEX( COMMAND, 'HOURLY TEMPERATURES' ) > 0 ) THEN

                WRITE( RPLCLINE,94020 ) 
     &               'HOURLY TEMPERATURES: ', TEMPS( CTYPOS,1:12 )
                SCENARIO( I ) = RPLCLINE
                
                WRITE( RPLCLINE,94030 ) TEMPS( CTYPOS,13:24 )
                SCENARIO( I + 1 ) = RPLCLINE
             END IF

C NEED TO HANDLE MIN/MAX REPLACEMENT ALSO

        END DO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( A21, 12( 1X, F7.3 ) )
94030   FORMAT( 12( F7.3, :, 1X ) )
        
        END SUBROUTINE RPLCTEMP
        