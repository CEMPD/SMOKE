
        SUBROUTINE WRM6INPUT( GDEV, SDEV, MDEV )

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        USE MODMBSET

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER           GETFLINE
        INTEGER           ENVINT
        INTEGER           GETSPDLN
        CHARACTER(LEN=80) WRSPDVMT
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  GETFLINE, ENVINT, GETSPDLN, WRSPDVMT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,           INTENT (IN) :: GDEV     ! GROUP file unit no.
        INTEGER,           INTENT (IN) :: SDEV     ! SPDSUM file unit no.
        INTEGER,           INTENT (IN) :: MDEV     ! M6INPUT file unit no.

C...........   Local allocatable arrays
        CHARACTER(LEN=FIPLEN3), ALLOCATABLE :: GRPLIST( :,: ) ! contents of GROUP file
        CHARACTER(LEN=150),     ALLOCATABLE :: M6SCEN( : )    ! M6 scenario file

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER JYEAR                     ! emission factor year
        INTEGER NLINES                    ! number of lines in GROUP file
        INTEGER NLINESCEN                 ! number of lines in M6 scenario file
        INTEGER NLINESPD                  ! number of lines in speed summary file
        INTEGER FDEV                      ! unit no. for M6 scenario file
        INTEGER :: SCENNUM = 1            ! current scenario number
        INTEGER STLINE                    ! starting line in SPDSUM for current county
        INTEGER CURRROAD                  ! road type from SPDSUM file
        INTEGER PREVROAD                  ! previous road type
        
        REAL    CURRSPD                   ! speed value from SPDSUM file
        REAL    PREVSPD                   ! previous speed

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        LOGICAL :: NEWSCEN    = .FALSE.   ! true: print current and create a new scenario
        
        CHARACTER(LEN=FIPLEN3) CURRCOUNTY            ! current county FIPS code
        CHARACTER(LEN=FIPLEN3) REFCOUNTY             ! ref. county FIPS code for curr. county
        CHARACTER(LEN=6)       SCENARIO              ! scenario number
        
        CHARACTER(LEN=60)      SPDDIR   ! directory for creating speed vmt files
        CHARACTER(LEN=80)      SPDFILE  ! name of SPEED VMT file for M6 input file
        CHARACTER(LEN=300)     SCENFILE !  M6 scenario file name
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRM6INPUT'   ! program name
        
C***********************************************************************
C   begin body of subroutine WRM6INPUT

C.........  Get speed vmt directory information from the environment
        MESG = 'Path where speed vmt files for MOBILE6 will be written'
        CALL ENVSTR( 'SMK_SPDPATH', MESG, '.', SPDDIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: Speed vmt files being placed '//
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_SPDPATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF
        
C.........  Get the year for computing the emission factors
        MESG = 'Emission factors year'
        JYEAR = ENVINT( 'EF_YEAR', MESG, 1988, IOS )

C.........  Write message about which year emission factors will be for
        WRITE( MESG,94010 ) 
     &         'NOTE: Emission factors are for year', JYEAR
        CALL M3MSG2( MESG )

C.........  Get number of lines in speed summary file
        NLINESPD = GETFLINE( SDEV, 'Speed summary file' )

C.........  Get the number of lines in the GROUP file     
        NLINES = GETFLINE( GDEV, 'GROUP county list file' )
        
        ALLOCATE( GRPLIST( NLINES,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPLIST', PROGNAME )
                        
C.........  Read GROUP file
        CALL RDGRPLIST( GDEV, NLINES, GRPLIST )

C.........  Loop through counties in GROUP file
        DO I = 1, NLINES

            CURRCOUNTY = GRPLIST( I,1 )
            REFCOUNTY  = GRPLIST( I,2 )

C.............  Open M6 scenario file for current county
            CALL OPENSCEN( REFCOUNTY, FDEV, SCENFILE )

C.............  Get number of lines in M6 scenario file
            MESG = 'MOBILE6 scenario file for county' // REFCOUNTY
            NLINESCEN = GETFLINE( FDEV, MESG )
            
            ALLOCATE( M6SCEN( NLINESCEN ), STAT=IOS )
            CALL CHECKMEM( IOS, 'M6SCEN', PROGNAME )
            
C.............  Read M6 scenario file into array            
            CALL RDLINES( FDEV, MESG, NLINESCEN, M6SCEN )

C.............  Check calendar year of M6 scenario
            CALL CHKSCNYR( SCENFILE, M6SCEN, NLINESCEN, JYEAR )
            
C.............  Replace temperatures in M6 scenario

C.............  Find starting line for current county in SPDSUM file
            STLINE = GETSPDLN( SDEV, CURRCOUNTY, NLINESPD )

            IF( STLINE == 0 ) THEN
            	EFLAG = .TRUE.
            	
                WRITE( MESG, 93000 ) 'ERROR: Could not find county ' 
     &                 // CURRCOUNTY // ' in speed summary file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

            PREVROAD = LOCAL

C.............  Read speeds and sources
            DO
                CALL RDSPDLINE( SDEV, SCENNUM, CURRCOUNTY, STLINE,
     &                          NLINESPD, CURRROAD, CURRSPD )
                
                IF( CURRROAD == LOCAL ) THEN
                    CURRSPD = PREVSPD
                END IF
                 
                IF( CURRROAD == ARTERIAL ) THEN
                    NEWSCEN = .TRUE.
                END IF

                IF( PREVROAD == CURRROAD ) THEN
                    NEWSCEN = .TRUE.
                END IF

                IF( PREVROAD == FREEWAY .AND. CURRROAD == LOCAL ) THEN
                    NEWSCEN = .TRUE.
                END IF
              
                IF( NEWSCEN ) THEN
                	
C.....................  Create speed vmt file
                    SPDFILE = WRSPDVMT( PREVSPD, CURRSPD, SPDDIR )
                    
                    IF( SPDFILE == '' ) THEN
                    	EFLAG = .TRUE.
                    	
                        WRITE( MESG,93000 ) 'ERROR: Could not create '
     &                         // 'speed vmt file: ' // CRLF() 
     &                         // BLANK5 // SPDFILE
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF
                    
C.....................  Write M6 scenario to M6INPUT                   
                    WRITE( SCENARIO, '(I6)' ) SCENNUM
                    CALL PADZERO( SCENARIO )
                    WRITE( MDEV, 93000 ) 
     &                      'SCENARIO RECORD    : ' // SCENARIO
                    DO J = 1, NLINESCEN
                        WRITE( MDEV, 93000 ) 
     &                      M6SCEN( J )( 1:LEN_TRIM( M6SCEN( J ) ) )
                    END DO
                    WRITE( MDEV, 93000 )
     &                      'SPEED VMT          : ' 
     &                      // SPDFILE( 1:LEN_TRIM( SPDFILE ) )

C.....................  Increment scenario number
                    SCENNUM = SCENNUM + 1

                END IF

                NEWSCEN = .FALSE.

                IF( CURRROAD == LOCAL ) EXIT

                PREVROAD = CURRROAD
                PREVSPD  = CURRSPD

            END DO        
        
            DEALLOCATE( M6SCEN )
        
        END DO

C.........  Abort if error found while writing M6 input file
        IF( EFLAG ) THEN
            MESG = 'Problem writing MOBILE6 input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, 1X, I1, 1X, F6.2, 7( 1X, I6 ), 1X, 1A )  

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, F6.2, :, 1X ) )

        END SUBROUTINE WRM6INPUT
        