
        SUBROUTINE RDHOURTEMP( TNAME, NCOUNTY, JDATE, JTIME, 
     &                         COUNTYTEMP )
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
                
C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=16), INTENT (IN) :: TNAME    ! logical name for temperature file
        INTEGER,           INTENT (IN) :: NCOUNTY  ! no. counties in temp file
        INTEGER,           INTENT (IN) :: JDATE    ! starting date (YYYYDDD)
        INTEGER,           INTENT (IN) :: JTIME    ! starting time (HHMMSS)
        REAL, INTENT(INOUT) :: COUNTYTEMP( NCOUNTY, 24 ) ! array of hourly temperatures by county

C...........   Local allocatable arrays

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDHOURTEMP'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDHOURTEMP
        
C.........  Initialize array
        COUNTYTEMP = 0.

C.........  Loop through the time steps        
        DO I = 1, 24
        
            IF( .NOT. READ3( TNAME, 'TKCOUNTY', 1, JDATE, JTIME, 
     &                       COUNTYTEMP( :,I ) ) ) THEN
                MESG = 'Could not read TKCOUNTY' //
     &                 ' from ' // TNAME 
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF
        
C.............  Convert temps from Kelvin to Fahrenheit 
            DO J = 1, NCOUNTY
                COUNTYTEMP( J,I ) = CTOF * 
     &                              ( COUNTYTEMP( J,I ) - CTOK ) + 32.
            END DO
                   
            CALL NEXTIME( JDATE, JTIME, 10000 )
        
        END DO

        END SUBROUTINE RDHOURTEMP
        