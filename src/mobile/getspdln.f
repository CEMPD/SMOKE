
        INTEGER FUNCTION GETSPDLN( SDEV, COUNTY, NLINES )
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   FUNCTION ARGUMENTS
        INTEGER,                INTENT (IN) :: SDEV       ! SPDSUM file unit no.
        CHARACTER(LEN=FIPLEN3), INTENT (IN) :: COUNTY     ! county to search for
        INTEGER,                INTENT (IN) :: NLINES     ! number of lines in SPDSUM file

C...........   Other local variables
        INTEGER I                          ! counters and indices                     

        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        
        LOGICAL :: FNDLINE = .FALSE.      ! true: found starting line for county
        
        CHARACTER(LEN=FIPLEN3) CURRCOUNTY  ! current county in SPDSUM file
        CHARACTER(LEN=100)     LINE     !  line buffer
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETSPDLN'   ! program name

C***********************************************************************
C   begin body of function GETSPDLN

        IREC = 0
        FNDLINE = .FALSE.

C.........  Loop through SPDSUM file
        DO I = 1, NLINES
        
C.........  Read line
            READ( SDEV, 93000, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS /= 0 ) THEN
                IF( IOS == -1 ) THEN
                    MESG = 'End of file reached unexpectedly. ' //
     &              'Check format of SPDSUM file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
                END IF
                
                WRITE( MESG, 94010 )
     &              'I/O error', IOS,
     &              'reading speed summary file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
                
            CURRCOUNTY = LINE( 1:6 )
            CALL PADZERO( CURRCOUNTY )
            
            IF( CURRCOUNTY == COUNTY ) THEN
                GETSPDLN = IREC
                FNDLINE = .TRUE.
                EXIT
            END IF

        END DO
        	
        IF( .NOT. FNDLINE ) THEN
            GETSPDLN = 0
        END IF
        	
        BACKSPACE( SDEV )
        	
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
        END FUNCTION GETSPDLN
        