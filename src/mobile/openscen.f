
        SUBROUTINE OPENSCEN( COUNTY, FDEV, FILENAME )

C.........  MODULES for public variables

        USE MODMBSET        

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER          JUNIT
        CHARACTER(LEN=2) CRLF 
        
        EXTERNAL JUNIT, CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=FIPLEN3), INTENT (IN)  :: COUNTY   ! ref. county
        INTEGER,                INTENT (OUT) :: FDEV     ! M6 scenario file unit no.
        CHARACTER(LEN=300),     INTENT (OUT) :: FILENAME ! name of M6 scenario file

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        LOGICAL :: FNDSCEN    = .FALSE.   ! true: found M6 scenario file

        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'OPENSCEN'   ! program name
        
C***********************************************************************
C   begin body of subroutine OPENSCEN
        
        FNDSCEN = .FALSE.
        
C.........  Find M6 scenario file in M6LIST
        DO J = 1, SIZE( M6LIST )
            K = INDEX( M6LIST( J ), COUNTY )  
           
            IF( K >= 1 ) THEN
                FILENAME = M6LIST( J )
                FNDSCEN = .TRUE.
                EXIT
            END IF                
        END DO
        
        IF( .NOT. FNDSCEN ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Could not find '
     &             // 'MOBILE6 scenario file for reference '
     &             // CRLF() // BLANK5 // 'county '
     &             // COUNTY // ' in M6LIST file.'
            CALL M3MESG( MESG )
        END IF
        
C.........  Open M6 scenario file
        IF( .NOT. EFLAG ) THEN
            FDEV = JUNIT()
     
            OPEN( FDEV, FILE=FILENAME( 1:LEN_TRIM( FILENAME ) ), 
     &            ACTION='READ', STATUS='OLD', IOSTAT=IOS )
     
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Could not open MOBILE6 '
     &                 // 'scenario file' // CRLF() // BLANK5
     &                 // FILENAME( 1:LEN_TRIM( FILENAME ) )
     &                 // CRLF() // BLANK5 // 'for county '
     &                 // COUNTY
                CALL M3MESG( MESG )
            END IF
        END IF

C.........  Abort if error found while opening scenario file
        IF( EFLAG ) THEN
            MESG = 'Problem opening MOBILE6 scenario file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE OPENSCEN