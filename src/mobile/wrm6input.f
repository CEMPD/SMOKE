
        SUBROUTINE WRM6INPUT( GRPLIST, NLINES, SDEV, MDEV, 
     &                        TEMPS, NCOUNTY, NSTEPS, VOLNAM, 
     &                        SCENNUM, SRCNUM )

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
        INTEGER           FIND1
        INTEGER           GETSPDLN
        INTEGER           STR2INT
        CHARACTER(LEN=80) WRSPDVMT
        CHARACTER(LEN=2)  CRLF    
        
        EXTERNAL  GETFLINE, ENVINT, FIND1, GETSPDLN, STR2INT, 
     &            WRSPDVMT, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN)  :: GRPLIST( NLINES,3 )          ! GROUP file contents
        INTEGER,      INTENT (IN)  :: NLINES                       ! no. lines in GROUP file
        INTEGER,      INTENT (IN)  :: SDEV                         ! SPDSUM file unit no.
        INTEGER,      INTENT (IN)  :: MDEV                         ! M6INPUT file unit no.
        REAL,         INTENT (IN)  :: TEMPS( NCOUNTY, 0:NSTEPS-1 ) ! temps per county
        INTEGER,      INTENT (IN)  :: NCOUNTY                      ! no. counties in temps array
        INTEGER,      INTENT (IN)  :: NSTEPS                       ! no. time steps in temps array
        CHARACTER(*), INTENT (IN)  :: VOLNAM                       ! volatile pollutant name
        INTEGER,      INTENT (OUT) :: SCENNUM                      ! total number of scenarios
        INTEGER,      INTENT (OUT) :: SRCNUM                       ! total number of sources

C...........   Local allocatable arrays
        CHARACTER(LEN=150),     ALLOCATABLE :: M6SCEN( : )    ! M6 scenario file
        INTEGER, ALLOCATABLE :: CTYLIST( : )   ! list of counties in temperature file

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER JYEAR                     ! emission factor year
        INTEGER NLINESCEN                 ! number of lines in M6 scenario file
        INTEGER NLINESPD                  ! number of lines in speed summary file
        INTEGER FDEV                      ! unit no. for M6 scenario file
        INTEGER PREVCTY                   ! previous county number in SPDSUM file
        INTEGER CURRCTY                   ! current county number in SPDSUM file
        INTEGER CTYNUM                    ! number of counties in SPDSUM file
        INTEGER CTYPOS                    ! position of current county in CTYLIST
        INTEGER STLINE                    ! starting line in SPDSUM for current county
        INTEGER CURRROAD                  ! road type from SPDSUM file
        INTEGER PREVROAD                  ! previous road type
        INTEGER STSCEN                    ! starting line of scenario data in scenario file
        INTEGER LASAFLAG                  ! local-as-arterial setting
        
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
        
C.........  Read county list from SPDSUM file
        ALLOCATE( CTYLIST( NCOUNTY ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTYLIST', PROGNAME )

        PREVCTY = 0
        CTYNUM  = 0

        DO I = 1, NLINESPD
            READ( SDEV,93010 ) CURRCTY
            
            IF( CURRCTY /= PREVCTY ) THEN
                CTYNUM = CTYNUM + 1
                CTYLIST( CTYNUM ) = CURRCTY
            END IF
            
            PREVCTY = CURRCTY
        END DO
        
        REWIND( SDEV )

C.........  Write Mobile6 header info to input file
        CALL WRM6HEADER( MDEV )
        
C.........  Loop through counties in GROUP file
        DO I = 1, NLINES

C.............  Get info from GROUP file
            WRITE( CURRCOUNTY, '(I6)' ) GRPLIST( I,1 )
            CALL PADZERO( CURRCOUNTY )
            
            WRITE( REFCOUNTY, '(I6)' ) GRPLIST( I,2 )
            CALL PADZERO( REFCOUNTY )
            
            LASAFLAG = GRPLIST( I,3 )

C.............  Open M6 scenario file for current county
            CALL OPENSCEN( REFCOUNTY, FDEV, SCENFILE )

C.............  Get number of lines in M6 scenario file
            MESG = 'MOBILE6 scenario file for county' // REFCOUNTY
            NLINESCEN = GETFLINE( FDEV, MESG )
            
            ALLOCATE( M6SCEN( NLINESCEN ), STAT=IOS )
            CALL CHECKMEM( IOS, 'M6SCEN', PROGNAME )
            
C.............  Read M6 scenario file into array            
            CALL RDLINES( FDEV, MESG, NLINESCEN, M6SCEN )

C.............  Close scenario file
            CLOSE( FDEV )

C.............  Find beginning of scenario commands 
            STSCEN = 0

            DO J = 1, NLINESCEN
                IF( M6SCEN( J )(1:20) == 'SCENARIO RECORD    :' ) THEN
                    STSCEN = J
                END IF
            END DO

            IF( STSCEN == 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 93000 ) 'ERROR: Scenario file for county '
     &                 // REFCOUNTY // ' does not include the '
     &                 // 'SCENARIO RECORD command.'
                CALL M3MESG( MESG )
                CYCLE
            END IF           

C.............  Check calendar year of M6 scenario
            CALL CHKSCNYR( SCENFILE, M6SCEN, NLINESCEN, JYEAR )

C.............  Find county in CTYLIST array
            CTYPOS = FIND1( STR2INT( CURRCOUNTY ), NCOUNTY, CTYLIST )
            
            IF( CTYPOS < 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 93000 ) 'ERROR: Could not find county '
     &                 // CURRCOUNTY // ' in speed summary file.'
                CALL M3MESG( MESG )
                CYCLE              
            END IF
            
C.............  Replace temperatures in M6 scenario
            CALL RPLCTEMP( CURRCOUNTY, TEMPS, NCOUNTY, NSTEPS, 
     &                     M6SCEN, NLINESCEN, CTYPOS )

C.............  Write run level commands to M6 input file
            DO J = 1, STSCEN - 1
                WRITE( MDEV, 93000 ) 
     &                      M6SCEN( J )( 1:LEN_TRIM( M6SCEN( J ) ) )
            END DO

C.............  Select M6 output based on volatile pollutant name            
            SELECT CASE( VOLNAM( 1:LEN_TRIM( VOLNAM ) )
            CASE ( 'THC' )
                WRITE( MDEV,93000 ) 'EXPRESS HC AS THC  :' // CRLF()
            CASE ( 'NMH' )
                WRITE( MDEV,93000 ) 'EXPRESS HC AS NMHC :' // CRLF()
            CASE ( 'TOG' )
                WRITE( MDEV,93000 ) 'EXPRESS HC AS TOG  :' // CRLF()
            CASE ( 'NMO' )
                WRITE( MDEV,93000 ) 'EXPRESS HC AS HMOG :' // CRLF()
            CASE DEFAULT
                WRITE( MDEV,93000 ) 'EXPRESS HC AS VOC  :' // CRLF()
            END SELECT

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
     &                          NLINESPD, LASAFLAG, CURRROAD, CURRSPD,
     &                          SRCNUM )
                
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
                    DO J = STSCEN, NLINESCEN
                        WRITE( MDEV, 93000 ) 
     &                      M6SCEN( J )( 1:LEN_TRIM( M6SCEN( J ) ) )
                    END DO
                    WRITE( MDEV, 93000 )
     &                      'SPEED VMT          : ' 
     &                      // SPDFILE( 1:LEN_TRIM( SPDFILE ) )
                    WRITE( MDEV, 93000 ) BLANK10
                    
C.....................  Increment scenario number
                    SCENNUM = SCENNUM + 1

                END IF

                NEWSCEN = .FALSE.

                IF( CURRROAD == LOCAL ) EXIT

                PREVROAD = CURRROAD
                PREVSPD  = CURRSPD

            END DO        
        
            DEALLOCATE( M6SCEN )
        
C.............  Write end of run command to M6 input file
            WRITE( MDEV,93000 ) 'END OF RUN         :' // CRLF()

        END DO

C.........  Abort if error found while writing M6 input file
        IF( EFLAG ) THEN
            MESG = 'Problem writing MOBILE6 input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, :, 1X, I1, 1X, F6.2, 7( 1X, I6 ), 1X, 1A )  

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 10( A, :, F6.2, :, 1X ) )

        END SUBROUTINE WRM6INPUT
        