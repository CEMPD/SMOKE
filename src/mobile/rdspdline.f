
        SUBROUTINE RDSPDLINE( SDEV, SCENNUM, CURRCOUNTY, STLINE,
     &                        NLINES, ROADTYPE, SPEED )

C.........  MODULES for public variables
        USE MODMBSET
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER           STR2INT
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  STR2INT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,                INTENT (IN)  :: SDEV        ! SPDSUM file unit number
        INTEGER,                INTENT (IN)  :: SCENNUM     ! current scenario number
        CHARACTER(LEN=FIPLEN3), INTENT (IN)  :: CURRCOUNTY  ! current county being processed
        INTEGER,                INTENT (IN)  :: STLINE      ! starting line for current county
        INTEGER,                INTENT (IN)  :: NLINES      ! number of lines in SPDSUM file
        INTEGER,                INTENT (OUT) :: ROADTYPE    ! road type from SPDSUM file
        REAL,                   INTENT (OUT) :: SPEED       ! speed value from SPDSUM file

C...........   Local arrays
        INTEGER       SOURCES( 7 )          ! line of sources from SPDSUM file

C...........   Other local variables
        INTEGER I, J                      ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER, SAVE :: IREC = 0         ! record counter
        INTEGER COUNTY                    ! county from SPDSUM file

        LOGICAL   :: EFLAG     = .FALSE. !  true: error found
                
        CHARACTER   CONTCHAR              ! continuation character from SPDSUM file
        
        CHARACTER(LEN=300)     MESG       !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSPDLINE'   ! program name

C***********************************************************************
C   begin body of subroutine RDSPDLINE

        SOURCES = 0
        IREC = IREC + STLINE - 1

        DO

C.............  Make sure we don't try to read past the end of the file
            IF( IREC + 1 >= NLINES ) THEN
            	ROADTYPE = LOCAL
            	SPEED = 0.
            	EXIT
            END IF

C.............  Read line from SPDSUM file
            READ( SDEV, 93010, IOSTAT=IOS ) COUNTY, ROADTYPE, 
     &            SPEED, SOURCES( 1:7 ), CONTCHAR
            
            IREC = IREC + 1
            
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                
                IF( IOS == -1 ) THEN
                    MESG = 'End of file reached unexpectedly. ' //
     &                     'Check format of SPDSUM file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )   
                END IF
                
                WRITE( MESG, 94010 )
     &                 'I/O error', IOS,
     &                 'reading speed summary file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Make sure that we've reached the current county's starting line
            IF( IREC < STLINE ) CYCLE

C.............  Check that county from SPDSUM is still current county                
            IF( STR2INT( CURRCOUNTY ) /= COUNTY ) THEN
                BACKSPACE( SDEV )
                IREC = IREC - 1
            	ROADTYPE = LOCAL
            	SPEED = 0.
                EXIT
            END IF

C.............  Store scenario numbers with sources
            DO J = 1, 7
                IF( SOURCES( J ) /= 0 ) THEN
                    SCENLIST( SOURCES( J ) ) = SCENNUM
                END IF
            END DO

C.............  Check if there are more lines to read for this speed and road type
            IF( CONTCHAR /= '\' ) EXIT

        END DO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, 1X, I1, 1X, F6.2, 7( 1X, I6 ), 1X, 1A )  

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, F6.2, :, 1X ) )
        
        END SUBROUTINE RDSPDLINE
        