
        PROGRAM CEMSCAN
        
C***********************************************************************
C  program body starts at line 39
C
C  DESCRIPTION: 
C      Program to read a year's worth of CEM data and calculate summed annual
C      NOx emissions, SO2 emissions, heat input, gross load, and steam load. 
C      The output from CEMScan is used by Smkinven to calculate hourly emissions
C      from annual inventory data when reading CEM data. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Written by C. Seppanen
C
C***********************************************************************

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
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv5.2.1_Sep2025$'  ! CVS release tag
        INTEGER, PARAMETER :: MXSEG = 16        ! number of segments in line

C.........  LOCAL VARIABLES

C.........  Static arrays
        CHARACTER(20) SEGMENT(MXSEG)            ! parsed input line
        
C.........  Allocatable arrays
        CHARACTER(OBRLEN3), ALLOCATABLE :: UNITLIST( : ) ! list of units
        INTEGER, ALLOCATABLE :: UNITIDX  ( : )  ! index to sorted unit list
        
        REAL,    ALLOCATABLE :: NOXEMIS  ( : )  ! NOx emissions by unit
        REAL,    ALLOCATABLE :: SO2EMIS  ( : )  ! SO2 emissions by unit
        REAL,    ALLOCATABLE :: OPTIMESUM( : )  ! operating time by unit
        REAL,    ALLOCATABLE :: GLOADSUM ( : )  ! gross load by unit
        REAL,    ALLOCATABLE :: SLOADSUM ( : )  ! s load by unit
        REAL,    ALLOCATABLE :: HEATINPUT( : )  ! heat input by unit
        
        REAL,    ALLOCATABLE :: MAXNOX   ( : )  ! max. NOx emissions by unit
        REAL,    ALLOCATABLE :: MAXSO2   ( : )  ! max. SO2 emissions by unit
        REAL,    ALLOCATABLE :: MAXOPTIME( : )  ! max. operating time by unit
        REAL,    ALLOCATABLE :: MAXGLOAD ( : )  ! max. gross load by unit
        REAL,    ALLOCATABLE :: MAXSLOAD ( : )  ! max. s load by unit
        REAL,    ALLOCATABLE :: MAXHEAT  ( : )  ! max. heat input by unit
        
        REAL,    ALLOCATABLE :: MINNOX   ( : )  ! min. NOx emissions by unit
        REAL,    ALLOCATABLE :: MINSO2   ( : )  ! min. SO2 emissions by unit
        REAL,    ALLOCATABLE :: MINOPTIME( : )  ! min. operating time by unit
        REAL,    ALLOCATABLE :: MINGLOAD ( : )  ! min. gross load by unit
        REAL,    ALLOCATABLE :: MINSLOAD ( : )  ! min. s load by unit
        REAL,    ALLOCATABLE :: MINHEAT  ( : )  ! min. heat input by unit
        
        INTEGER, ALLOCATABLE :: NUMHOURS ( : )  ! total no. hours by unit
        INTEGER, ALLOCATABLE :: NOXHOURS ( : )  ! no. hours with NOx emissions
        INTEGER, ALLOCATABLE :: SO2HOURS ( : )  ! no. hours with SO2 emissions
        INTEGER, ALLOCATABLE :: OPHOURS  ( : )  ! no. hours with operating time
        INTEGER, ALLOCATABLE :: GLDHOURS ( : )  ! no. hours with gross load
        INTEGER, ALLOCATABLE :: SLDHOURS ( : )  ! no. hours with s load
        INTEGER, ALLOCATABLE :: HEATHOURS( : )  ! no. hours with heat input

C.........  File units
        INTEGER LDEV                            ! file unit for log file
        INTEGER IDEV                            ! file unit for input file list
        INTEGER TDEV                            ! tmp. file unit for input files
        INTEGER ODEV                            ! file unit output file
        INTEGER RDEV                            ! file unit for report file

C.........  Other local variables
        INTEGER I, J, K, L                      ! indices and counters
        INTEGER IDX                             ! unit index
        INTEGER IOS                             ! I/O status
        INTEGER NUNITS                          ! total number of units
        INTEGER MXFILES                         ! max. number of input files
        INTEGER MXUNITS                         ! max. number of units
        
        REAL    NOX                             ! tmp. NOx emissions
        REAL    SO2                             ! tmp. SO2 emissions
        REAL    OPTIME                          ! tmp. operating time
        REAL    GLOAD                           ! tmp. gross load
        REAL    SLOAD                           ! tmp. s load
        REAL    HTINPUT                         ! tmp. heat input
        
        LOGICAL :: EFLAG = .FALSE.              ! true: an error has occurred
        
        CHARACTER          OUT                  ! Y: entry should be output
        CHARACTER(OBRLEN3) UNIT                 ! tmp. unit string
        CHARACTER(ORSLEN3) ORIS                 ! tmp. ORIS ID
        CHARACTER(BLRLEN3) BLRID                ! tmp. boiler ID
        CHARACTER(300)     LINE                 ! line buffer
        CHARACTER(300)     MESG                 ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'CEMSCAN'   ! program name

C***********************************************************************
C   begin body of program CEMSCAN

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
        ALLOCATE( UNITIDX( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNITIDX', PROGNAME )
        
        ALLOCATE( NOXEMIS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NOXEMIS', PROGNAME )
        ALLOCATE( SO2EMIS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SO2EMIS', PROGNAME )
        ALLOCATE( OPTIMESUM( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OPTIMESUM', PROGNAME )
        ALLOCATE( GLOADSUM( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GLOADSUM', PROGNAME )
        ALLOCATE( SLOADSUM( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLOADSUM', PROGNAME )
        ALLOCATE( HEATINPUT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HEATINPUT', PROGNAME )
        
        ALLOCATE( MAXNOX( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXNOX', PROGNAME )
        ALLOCATE( MAXSO2( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXSO2', PROGNAME )
        ALLOCATE( MAXOPTIME( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXOPTIME', PROGNAME )
        ALLOCATE( MAXGLOAD( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXGLOAD', PROGNAME )
        ALLOCATE( MAXSLOAD( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXSLOAD', PROGNAME )
        ALLOCATE( MAXHEAT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXHEAT', PROGNAME )
        
        ALLOCATE( MINNOX( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINNOX', PROGNAME )
        ALLOCATE( MINSO2( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINSO2', PROGNAME )
        ALLOCATE( MINOPTIME( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINOPTIME', PROGNAME )
        ALLOCATE( MINGLOAD( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINGLOAD', PROGNAME )
        ALLOCATE( MINSLOAD( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINSLOAD', PROGNAME )
        ALLOCATE( MINHEAT( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINHEAT', PROGNAME )
        
        ALLOCATE( NUMHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NUMHOURS', PROGNAME )
        ALLOCATE( NOXHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NOXHOURS', PROGNAME )
        ALLOCATE( SO2HOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SO2HOURS', PROGNAME )
        ALLOCATE( OPHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OPHOURS', PROGNAME )
        ALLOCATE( GLDHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GLDHOURS', PROGNAME )
        ALLOCATE( SLDHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLDHOURS', PROGNAME )
        ALLOCATE( HEATHOURS( MXUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HEATHOURS', PROGNAME )

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

            MESG = 'Successfully opened input file: ' // CRLF() //
     &             BLANK10 // TRIM( LINE )
            CALL M3MSG2( MESG )

C.............  Read through lines in file            
            DO
                READ( TDEV, 93000, IOSTAT=IOS ) LINE

C.................  Check for I/O errors                
                IF( IOS > 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: I/O error ', IOS,
     &                  'reading input file'
                    CALL M3MESG( MESG )
                    EXIT
                END IF

C.................  Check for end of file
                IF( IOS < 0 ) EXIT

C.................  Skip blank or comment lines
                IF( BLKORCMT( LINE ) ) CYCLE
                
                CALL PARSLINE( LINE, MXSEG, SEGMENT )
                    
C.................  Store data from file
                UNIT = ''
                UNIT = ADJUSTR( SEGMENT( 1 )( 1:ORSLEN3 ) ) // 
     &                 ADJUSTR( SEGMENT( 2 )( 1:BLRLEN3 ) )
                
                NOX     = STR2REAL( SEGMENT( 5 ) )
                SO2     = STR2REAL( SEGMENT( 6 ) )
                OPTIME  = STR2REAL( SEGMENT( 8 ) )
                GLOAD   = STR2REAL( SEGMENT( 9 ) )
                SLOAD   = STR2REAL( SEGMENT( 10 ) )
                HTINPUT = STR2REAL( SEGMENT( 11 ) )

C.................  Determine where current unit should go in sorted list of 
C                   units and check if it is already present
                IF( NUNITS == 0 ) THEN
                    K = 1
                ELSE
                    K = LOCATC( UNIT, NUNITS, UNITLIST )
                END IF

                IF( K == -1 ) THEN
                    K = FINDC( UNIT, NUNITS, UNITLIST )
                    
                    IDX = UNITIDX( K )
                    
                    NUMHOURS( IDX ) = NUMHOURS( IDX ) + 1
                    
                    CALL STORE_DATA( NOX, NOXEMIS( IDX ), 
     &                               NOXHOURS( IDX ), MAXNOX( IDX ), 
     &                               MINNOX( IDX ) )
                    CALL STORE_DATA( SO2, SO2EMIS( IDX ), 
     &                               SO2HOURS( IDX ), MAXSO2( IDX ), 
     &                               MINSO2( IDX ) )
                    CALL STORE_DATA( OPTIME, OPTIMESUM( IDX ), 
     &                               OPHOURS( IDX ), MAXOPTIME( IDX ), 
     &                               MINOPTIME( IDX ) )
                    CALL STORE_DATA( GLOAD, GLOADSUM( IDX ), 
     &                               GLDHOURS( IDX ), MAXGLOAD( IDX ), 
     &                               MINGLOAD( IDX ) )
                    CALL STORE_DATA( SLOAD, SLOADSUM( IDX ), 
     &                               SLDHOURS( IDX ), MAXSLOAD( IDX ), 
     &                               MINSLOAD( IDX ) )
                    CALL STORE_DATA( HTINPUT, HEATINPUT( IDX ), 
     &                               HEATHOURS( IDX ), MAXHEAT( IDX ), 
     &                               MINHEAT( IDX ) )
                    
                    CYCLE
                END IF

C.................  Make sure we have space to add unit                
                IF( NUNITS + 1 > MXUNITS ) THEN
                    MESG = 'MAX_CEM_UNITS value too small'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE

C.....................  Move existing entries down in unit list to make room
                    DO L = NUNITS, K, -1
                        UNITLIST( L+1 ) = UNITLIST( L )
                        UNITIDX ( L+1 ) = UNITIDX ( L )
                    END DO

C.....................  Store data                    
                    UNITLIST ( K ) = UNIT
                    
                    NUNITS = NUNITS + 1
                    IDX = NUNITS
                    
                    UNITIDX ( K ) = IDX
                    NUMHOURS( IDX ) = 1
                    
                    CALL INIT_DATA( NOX, NOXEMIS( IDX ), 
     &                              NOXHOURS( IDX ), MAXNOX( IDX ), 
     &                              MINNOX( IDX ) )
                    CALL INIT_DATA( SO2, SO2EMIS( IDX ), 
     &                              SO2HOURS( IDX ), MAXSO2( IDX ), 
     &                              MINSO2( IDX ) )
                    CALL INIT_DATA( OPTIME, OPTIMESUM( IDX ), 
     &                              OPHOURS( IDX ), MAXOPTIME( IDX ), 
     &                              MINOPTIME( IDX ) )
                    CALL INIT_DATA( GLOAD, GLOADSUM( IDX ), 
     &                              GLDHOURS( IDX ), MAXGLOAD( IDX ), 
     &                              MINGLOAD( IDX ) )
                    CALL INIT_DATA( SLOAD, SLOADSUM( IDX ), 
     &                              SLDHOURS( IDX ), MAXSLOAD( IDX ), 
     &                              MINSLOAD( IDX ) )
                    CALL INIT_DATA( HTINPUT, HEATINPUT( IDX ), 
     &                              HEATHOURS( IDX ), MAXHEAT( IDX ), 
     &                              MINHEAT( IDX ) )
                END IF
            END DO
            
            CLOSE( TDEV )
        
        END DO
        
C.........  Write output file header
        WRITE( ODEV, 93000 ) '#ORIS  BOILER  NOX               ' //
     &    'SO2               OPTIME            GLOAD' //
     &    '            SLOAD              HTINPUT'

C.........  Write report header
        WRITE( RDEV, 93000 ) 'ORIS ID; Boiler ID; Out; Tot Hrs; ' //
     &    ' NOx Hrs;      Ann NOx Emis;      Max NOx Emis; ' //
     &    '     Min NOx Emis;  SO2 Hrs;' //
     &    '      Ann SO2 Emis;      Max SO2 Emis;      ' //
     &    'Min SO2 Emis;   Op Hrs;       Ann Op Time;' //
     &    '       Max Op Time;       Min Op Time; ' //
     &    ' Gld Hrs;         Ann GLOAD;         Max GLOAD;    ' //
     &    '     Min GLOAD;  Sld Hrs;    ' //
     &    '     Ann SLOAD;         Max SLOAD;         Min SLOAD; ' //
     &    'Heat Hrs;      Ann Ht Input;      Max Ht Input; ' //
     &    '     Min Ht Input'

C.........  Write output and report file
        DO I = 1, NUNITS
            
            ORIS  = UNITLIST( I )( 1:ORSLEN3 )
            BLRID = UNITLIST( I )( ORSLEN3+1:OBRLEN3 )
            IDX   = UNITIDX ( I )
            
            OUT = 'Y'
            IF( NOXEMIS  ( IDX ) <= 0. .AND. 
     &          HEATINPUT( IDX ) <= 0. .AND.
     &          SLOADSUM ( IDX ) <= 0. .AND.
     &          GLOADSUM ( IDX ) <= 0.       ) OUT = 'N'
            
            IF( OUT == 'Y' ) THEN
                WRITE( ODEV,93010 ) ORIS, BLRID, NOXEMIS( IDX ), 
     &            SO2EMIS( IDX ), OPTIMESUM( IDX ), GLOADSUM( IDX ),
     &            SLOADSUM( IDX ), HEATINPUT( IDX )
            END IF
            
            WRITE( RDEV,93020 ) ORIS, BLRID, OUT, NUMHOURS( IDX ), 
     &        NOXHOURS( IDX ), NOXEMIS( IDX ), MAXNOX( IDX ), 
     &            MINNOX( IDX ),
     &        SO2HOURS( IDX ), SO2EMIS( IDX ), MAXSO2( IDX ), 
     &            MINSO2( IDX ),
     &        OPHOURS( IDX ), OPTIMESUM( IDX ), MAXOPTIME( IDX ), 
     &            MINOPTIME( IDX ),
     &        GLDHOURS( IDX ), GLOADSUM( IDX ), MAXGLOAD( IDX ), 
     &            MINGLOAD( IDX ),
     &        SLDHOURS( IDX ), SLOADSUM( IDX ), MAXSLOAD( IDX ), 
     &            MINSLOAD( IDX ),
     &        HEATHOURS( IDX ), HEATINPUT( IDX ), MAXHEAT( IDX ), 
     &            MINHEAT( IDX )
        END DO

        IF( EFLAG ) THEN
            MESG = 'See errors listed above.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
        END IF

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( A6, 1X, A6, 6( 1X, E17.10) )
93020   FORMAT( 1X, A6, ';', 4X, A6, ';', 3X, A1, ';', 4X, I4, ';', 
     &          6( 5X, I4, ';', 3( 1X, E17.10, ';' ) ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS
        
            SUBROUTINE INIT_DATA( VAL, SUM, HOURS, MAX, MIN )            
            
C.............  Subroutine arguments
            REAL,    INTENT(IN)     :: VAL    ! data value
            REAL,    INTENT(IN OUT) :: SUM    ! summed data
            INTEGER, INTENT(IN OUT) :: HOURS  ! no. valid hours
            REAL,    INTENT(IN OUT) :: MAX    ! max. data val
            REAL,    INTENT(IN OUT) :: MIN    ! min. data val

C----------------------------------------------------------------------

            HOURS = 0
            IF( VAL >= 0. ) HOURS = 1
            
            SUM = VAL
            MAX = VAL
            MIN = VAL

            END SUBROUTINE INIT_DATA

C----------------------------------------------------------------------

            SUBROUTINE STORE_DATA( VAL, SUM, HOURS, MAX, MIN )
            
C.............  Subroutine arguments
            REAL,    INTENT(IN)     :: VAL    ! data value
            REAL,    INTENT(IN OUT) :: SUM    ! summed data
            INTEGER, INTENT(IN OUT) :: HOURS  ! no. valid hours
            REAL,    INTENT(IN OUT) :: MAX    ! max. data val
            REAL,    INTENT(IN OUT) :: MIN    ! min. data val

C----------------------------------------------------------------------

            IF( VAL >= 0. ) THEN
                IF( SUM < 0. ) SUM = 0.
                SUM = SUM + VAL
                HOURS = HOURS + 1
            END IF
            
            IF( MAX < VAL ) MAX = VAL
            IF( MIN > VAL ) MIN = VAL

            RETURN

            END SUBROUTINE STORE_DATA
            
        END PROGRAM CEMSCAN
