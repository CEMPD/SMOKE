
        SUBROUTINE RDSPDSRC( SDEV, NSRC, SRCARRAY )
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER           STR2INT
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  STR2INT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: SDEV               ! SPDSUM file unit number
        INTEGER, INTENT (IN)  :: NSRC               ! no. of sources
        INTEGER, INTENT (OUT) :: SRCARRAY( NSRC,2 ) ! array to hold county codes

C...........   Local arrays
        INTEGER       SOURCES( 7 )          ! line of sources from SPDSUM file

C...........   Other local variables
        INTEGER I, J                      ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER :: IREC = 0               ! record counter
        INTEGER COUNTY                    ! county from SPDSUM file
        INTEGER ROADTYPE                  ! roadtype from SPDSUM file
         
        REAL SPEED                        ! speed from SPDSUM file

        LOGICAL   :: EFLAG     = .FALSE.  !  true: error found
                
        CHARACTER   CONTCHAR              ! continuation character from SPDSUM file
        
        CHARACTER(LEN=300)     MESG       !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSPDSRC'   ! program name

C***********************************************************************
C   begin body of subroutine RDSPDSRC

        DO
        
            SOURCES = 0   ! array
        
C.............  Read line from SPDSUM file
            READ( SDEV, 93010, IOSTAT=IOS, END=10 ) COUNTY, ROADTYPE, 
     &            SPEED, SOURCES, CONTCHAR

C.............  Exit if we've reached the end of the file
            IF( IOS == -1 ) EXIT

            IREC = IREC + 1

C.............  Check for other I/O errors        
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.

                WRITE( MESG, 94010 )
     &                 'I/O error', IOS,
     &                 'reading speed summary file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Store county code for each source
            DO J = 1, 7
                IF( SOURCES( J ) /= 0 ) THEN
                    SRCARRAY( SOURCES( J ),1 ) = COUNTY
                END IF
            END DO
                        
        END DO

C.........  Abort if error found while reading SPDSUM file
        IF( EFLAG ) THEN
            MESG = 'Problem reading SPDSUM file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

10      RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, 1X, I1, 1X, F6.2, 7( 1X, I6 ), 1X, 1A )  

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, F6.2, :, 1X ) )
        
        END SUBROUTINE RDSPDSRC
        