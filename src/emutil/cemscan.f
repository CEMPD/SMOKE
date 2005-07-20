C TODO:
C   calculate replacement heat input values
C   report min and max emissions values

        PROGRAM CEMSCAN
        
        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters

C.........  EXTERNAL FUNCTIONS
        LOGICAL       BLKORCMT
        CHARACTER(2)  CRLF
        INTEGER       ENVINT
        INTEGER       FINDC
        INTEGER       GETFLINE
        INTEGER       JUNIT
        INTEGER       LOCATC
        INTEGER       PROMPTFFILE
        REAL          STR2REAL
        
        EXTERNAL      BLKORCMT, CRLF, ENVINT, FINDC, GETFLINE, JUNIT, 
     &                LOCATC, PROMPTFFILE, STR2REAL

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$'  ! CVS release tag
        INTEGER, PARAMETER :: MXSEG = 12        ! number of segments in line

C.........  LOCAL VARIABLES

C.........  Static arrays
        CHARACTER(20) SEGMENT(MXSEG)            ! parsed input line
        
C.........  Allocatable arrays
        CHARACTER(ORSLEN3+BLRLEN3), ALLOCATABLE :: UNITLIST( : ) ! list of units
        REAL,    ALLOCATABLE :: HEATINPUT( : )  ! heat input by unit
        REAL,    ALLOCATABLE :: FLOWRATE ( : )  ! flow rate by unit
        REAL,    ALLOCATABLE :: MAXHEAT  ( : )  ! max. heat input by unit
        REAL,    ALLOCATABLE :: MINHEAT  ( : )  ! min. heat input by unit
        REAL,    ALLOCATABLE :: MAXFLOW  ( : )  ! max. flow rate by unit
        REAL,    ALLOCATABLE :: MINFLOW  ( : )  ! min. flow rate by unit
        INTEGER, ALLOCATABLE :: NUMHOURS ( : )  ! total no. hours by unit
        INTEGER, ALLOCATABLE :: HEATHOURS( : )  ! no. hours with "good" heat input
        INTEGER, ALLOCATABLE :: FLOWHOURS( : )  ! no. hours for calculating ave. flow

C.........  File units
        INTEGER LDEV                            ! file unit for log file
        INTEGER IDEV                            ! file unit for input file list
        INTEGER TDEV                            ! tmp. file unit for input files
        INTEGER ODEV                            ! file unit output file
        INTEGER RDEV                            ! file unit for report file

C.........  Other local variables
        INTEGER I, J, K, L                      ! indices and counters
        INTEGER IOS                             ! I/O status
        INTEGER NUNITS                          ! total number of units
        INTEGER MXFILES                         ! max. number of input files
        INTEGER MXUNITS                         ! max. number of units
        INTEGER NLINES                          ! number of lines in input file
        
        REAL    HTINPUT                         ! tmp. heat input
        REAL    STKFL                           ! tmp. flow rate
        
        LOGICAL :: EFLAG = .FALSE.              ! true: an error has occurred
        LOGICAL :: FLOWCHECK = .TRUE.           ! true: check if input has flow data
        LOGICAL :: CALCFLOW                     ! true: CEM files have flow data
        LOGICAL :: FLOWBAD                      ! true: flow rate is bad for current line
        LOGICAL :: HEATBAD                      ! true: heat input is bad for current line
        
        CHARACTER(ORSLEN3+BLRLEN3) UNIT         ! tmp. unit string
        CHARACTER(ORSLEN3) ORIS                 ! tmp. ORIS ID
        CHARACTER(BLRLEN3) BLRID                ! tmp. boiler ID
        CHARACTER(300)     LINE                 ! line buffer
        CHARACTER(256)     MESG                 ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'CEMSCAN'   ! program name

C***********************************************************************
C   begin body of program METCOMBINE

        LDEV = INIT3()
        
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get program setting environment variables
        MESG = 'Maximum number of CEM units'
        MXUNITS = ENVINT( 'MAX_CEM_UNITS', MESG, 0, IOS )
        IF( IOS /= 0 ) THEN
            MESG = 'Environment variable MAX_CEM_UNITS must be set'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for storing CEM data          
        ALLOCATE( UNITLIST( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNITLIST', PROGNAME )
        ALLOCATE( HEATINPUT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HEATINPUT', PROGNAME )
        ALLOCATE( FLOWRATE( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FLOWRATE', PROGNAME )
        
        ALLOCATE( NUMHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NUMHOURS', PROGNAME )
        ALLOCATE( HEATHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HEATHOURS', PROGNAME )
        ALLOCATE( FLOWHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FLOWHOURS', PROGNAME )
        
        ALLOCATE( MAXHEAT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXHEAT', PROGNAME )
        ALLOCATE( MINHEAT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINHEAT', PROGNAME )
        ALLOCATE( MAXFLOW( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXFLOW', PROGNAME )
        ALLOCATE( MINFLOW( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINFLOW', PROGNAME )

C.........  Open list of input files
        MESG = 'Enter logical name for CEM inputs list'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )

C.........  Open output files
        MESG = 'Enter logical name for OUTPUT file'
        ODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'OUTFILE', PROGNAME )
        
        MESG = 'Enter logical name for REPORT file'
        RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'REPFILE', PROGNAME )

C.........  Determine maximum number of input files
        MXFILES = GETFLINE( IDEV, 'List of CEM input files' )
        
        TDEV = JUNIT()
        NUNITS = 0
        
        DO I = 1, MXFILES
            READ( IDEV, 93000, IOSTAT=IOS ) LINE
            
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: I/O error ', IOS,
     &              'reading list of input files at line ', I
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Open input file
            OPEN( UNIT=TDEV, FILE=TRIM( LINE ), STATUS='OLD',
     &            ACTION='READ', IOSTAT=IOS )
     
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not open input file ' //
     &                 CRLF() // BLANK10 // TRIM( LINE )
                CALL M3MESG( MESG )
                CYCLE
            END IF

            MESG = 'Successfully opened input file: ' // TRIM( LINE )
            CALL M3MSG2( MESG )

C.............  Get number of lines in current file            
            NLINES = GETFLINE( TDEV, 'CEM input file' )

C.............  Read through lines in file            
            DO J = 1, NLINES
                READ( TDEV, 93000, IOSTAT=IOS ) LINE
                
                IF( IOS /= 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: I/O error'
                    CALL M3MESG( MESG )
                    EXIT
                END IF

C.................  Skip blank or comment lines
                IF( BLKORCMT( LINE ) ) CYCLE
                
                ! temporarily skip 1st line until files are changed
                if( j == 1 ) cycle

C.................  Reset data checks and parse line
                FLOWBAD = .FALSE.
                HEATBAD = .FALSE.
                
                CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.................  After parsing first line, check if it has flow                
                IF( FLOWCHECK ) THEN
                    IF( SEGMENT( 12 ) == '' ) THEN
                        CALCFLOW = .FALSE.
                        MESG = 'NOTE: Program will not calculate ' //
     &                    'average hourly flow since CEM data' //
     &                    BLANK10 // 'does not contain flow rates.'
                        CALL M3MSG2( MESG )
                    ELSE
                        CALCFLOW = .TRUE.
                        MESG = 'NOTE: Program will calculate ' //
     &                    'average hourly flow.'
                        CALL M3MSG2( MESG )
                    END IF
                    
                    FLOWCHECK = .FALSE.
                END IF
                    
C.................  Store data from file
                UNIT = ''
                UNIT = SEGMENT( 1 )( 1:ORSLEN3 ) // 
     &                 SEGMENT( 2 )( 1:BLRLEN3 )
                
                HTINPUT = STR2REAL( SEGMENT( 11 ) )
                IF( HTINPUT < 0. ) HEATBAD = .TRUE.
                
                IF( CALCFLOW ) THEN
                    STKFL = STR2REAL( SEGMENT( 12 ) )
                    IF( STKFL < 0. ) FLOWBAD = .TRUE.
                END IF

C.................  Determine where current unit should go in sorted list of 
C                   units and check if it is already present
                IF( NUNITS == 0 ) THEN
                    K = 1
                ELSE
                    K = LOCATC( UNIT, NUNITS, UNITLIST )
                END IF

                IF( K == -1 ) THEN
                    K = FINDC( UNIT, NUNITS, UNITLIST )
                    
                    NUMHOURS ( K ) = NUMHOURS( K ) + 1
                    
                    IF( .NOT. HEATBAD ) THEN
                        HEATINPUT( K ) = HEATINPUT( K ) + HTINPUT
                        HEATHOURS( K ) = HEATHOURS( K ) + 1
                    END IF
                    
                    IF( MAXHEAT( K ) < HTINPUT ) THEN
                        MAXHEAT( K ) = HTINPUT
                    END IF
                    
                    IF( MINHEAT( K ) > HTINPUT ) THEN
                        MINHEAT( K ) = HTINPUT
                    END IF
                    
                    IF( CALCFLOW ) THEN
                        IF( .NOT. FLOWBAD ) THEN
                            FLOWRATE ( K ) = FLOWRATE ( K ) + STKFL
                            FLOWHOURS( K ) = FLOWHOURS( K ) + 1
                        END IF
                    
                        IF( MAXFLOW( K ) < STKFL ) THEN
                            MAXFLOW( K ) = STKFL
                        END IF
                    
                        IF( MINFLOW( K ) > STKFL ) THEN
                            MINFLOW( K ) = STKFL
                        END IF
                    END IF
                    
                    CYCLE
                END IF

C.................  Make sure we have space to add unit                
                IF( NUNITS + 1 > MXUNITS ) THEN
                    MESG = 'MAX_CEM_UNITS too small'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE

C.....................  Move existing entries down in arrays to make room
                    DO L = NUNITS, K, -1
                        UNITLIST ( L+1 ) = UNITLIST ( L )
                        HEATINPUT( L+1 ) = HEATINPUT( L )
                        FLOWRATE ( L+1 ) = FLOWRATE ( L )
                        
                        NUMHOURS ( L+1 ) = NUMHOURS ( L )
                        HEATHOURS( L+1 ) = HEATHOURS( L )
                        FLOWHOURS( L+1 ) = FLOWHOURS( L )
                        
                        MAXHEAT  ( L+1 ) = MAXHEAT  ( L )
                        MINHEAT  ( L+1 ) = MINHEAT  ( L )
                        MAXFLOW  ( L+1 ) = MAXFLOW  ( L )
                        MINFLOW  ( L+1 ) = MINFLOW  ( L )
                    END DO
                    
                    UNITLIST ( K ) = UNIT
                    NUMHOURS ( K ) = 1
                    
                    IF( .NOT. HEATBAD ) THEN
                        HEATINPUT( K ) = HTINPUT
                        HEATHOURS( K ) = 1
                    ELSE
                        HEATINPUT( K ) = 0.
                        HEATHOURS( K ) = 0
                    END IF
                    
                    MAXHEAT  ( K ) = HTINPUT
                    MINHEAT  ( K ) = HTINPUT
                    
                    IF( CALCFLOW ) THEN
                        IF( .NOT. FLOWBAD ) THEN
                            FLOWRATE ( K ) = STKFL
                            FLOWHOURS( K ) = 1
                        ELSE
                            FLOWRATE ( K ) = 0.
                            FLOWHOURS( K ) = 0
                        END IF
                        
                        MAXFLOW  ( K ) = STKFL
                        MINFLOW  ( K ) = STKFL
                    END IF
                    
                    NUNITS = NUNITS + 1
                END IF
            END DO
            
            CLOSE( TDEV )
        
        END DO
        
C.........  Write report header
        MESG = 'ORIS ID; Boiler ID; Tot Hrs; Heat Hrs; ' //
     &         'Ann Heat Input; Max Heat Input; Min Heat Input; ' //
     &         'Flow Hrs; Ave Flow Rate; Max Flow Rate; Min Flow Rate'
        WRITE( RDEV, '(A)' ) TRIM( MESG )

C.........  Write output and report file
        DO I = 1, NUNITS
            IF( FLOWRATE( I ) /= 0. ) THEN
                STKFL = FLOWRATE( I ) / FLOWHOURS( I )
            ELSE
                STKFL = 0.
            END IF
            
            ORIS  = UNITLIST( I )( 1:ORSLEN3 )
            BLRID = UNITLIST( I )( ORSLEN3+1:ORSLEN3+BLRLEN3 )
            
            WRITE( ODEV,93010 ) ORIS, BLRID, HEATINPUT( I ), STKFL
            
            WRITE( RDEV,93020 ) ORIS, BLRID, NUMHOURS( I ), 
     &          HEATHOURS( I ), HEATINPUT( I ), MAXHEAT( I ), 
     &          MINHEAT( I ), FLOWHOURS( I ), STKFL, MAXFLOW( I ), 
     &          MINFLOW( I )
        END DO

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( A6, 1X, A6, 1X, E12.5, 1X, E12.5 )
93020   FORMAT( 1X, A6, ';', 4X, A6, ';', 4X, I4.4, ';', 5X, I4.4, ';',
     &          3( 3X, E12.5, ';' ), 5X, I4.4, ';', 2( 2X, E12.5, ';' ),
     &          2X, E12.5 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END PROGRAM CEMSCAN
