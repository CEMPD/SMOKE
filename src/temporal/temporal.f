
        PROGRAM TEMPORAL

C***********************************************************************
C  program body starts at line 214
C
C  DESCRIPTION:
C    This program computes the hourly emissions data from inventory emissions 
C    and/or activity and emission factor data. It can read average-inventory,
C    day-specific and hour-specific emissions and activity data.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    copied by: M. Houyoux 01/99
C    origin: tmppoint.F 4.3
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC  
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the temporal cross-reference tables
        USE MODXREF

C.........  This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C.........  This module is used for MOBILE6 setup information 
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN

        IMPLICIT NONE
 
C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C..........  EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         RDTPROF
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        LOGICAL         SETENVVAR

        EXTERNAL    CHKINT, CRLF, ENVINT, ENVYN, FINDC, 
     &              GETDATE, GETFLINE, GETNUM, INDEX1, MMDDYY, RDTPROF, 
     &              PROMPTMFILE, SECSDIFF, STR2INT, SETENVVAR
                        
C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C.........  Emission arrays
        REAL   , ALLOCATABLE :: EMAC ( :,: ) !  inven emissions or activities
        REAL   , ALLOCATABLE :: EMACV( :,: ) !  day-adjst emis or activities
        REAL   , ALLOCATABLE :: EMIST( :,: ) !  timestepped output emssions

C.........  Emission factor arrays        
        CHARACTER(LEN=256), ALLOCATABLE :: EFLIST( : )  ! listing of emission factor file names
        CHARACTER(LEN=16) , ALLOCATABLE :: EFLOGS( : )  ! listing of ef logical file names
        INTEGER           , ALLOCATABLE :: EFDAYS( :,: )! ef file by day for each time period
        REAL              , ALLOCATABLE :: EMFAC ( :,: )! mobile emission factors by source
        REAL              , ALLOCATABLE :: TEMPEF( : )  ! temporary holding array for efs
        CHARACTER(LEN=1)  , ALLOCATABLE :: EFTYPE( : )  ! ef file type (day, week, etc.) for each src
        INTEGER           , ALLOCATABLE :: EFIDX ( : )  ! location of ef in file for each source
        INTEGER           , ALLOCATABLE :: SRCS  ( : )  ! temporary array for sources in each ef file

C.........  Temporal allocation Matrix.  
        REAL, ALLOCATABLE :: TMAT( :, :, : ) ! temporal allocation factors

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C.........  Actual-SCC  table
        INTEGER                                NSCC
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE :: SCCLIST( : )

C.........  Day-specific, hour-specific data, and elevated sources data. 
C.........  These need only to allow enough dimensions for one read per 
C           pollutant per time step.

        INTEGER                 NPELV        ! optional elevated source-count
        INTEGER, ALLOCATABLE :: INDXE( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISE( : )   ! elevated source emissions

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : )   ! contiguous ungridding matrix

C.........  Names of pollutants and activities associated with output variables
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: ALLIN( : ) 

C.........  Reshaped input variables and output variables
        INTEGER         NGRP                ! no. of pol/emis-types groups 
        INTEGER         NGSZ                ! no. of pols/emis-types per group 
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: ALLIN2D( :,: ) 
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EANAM2D( :,: ) 

C...........   Logical names and unit numbers

        INTEGER      :: CDEV = 0!  unit number for region codes file
        INTEGER      :: EDEV = 0!  unit number for ef file list
        INTEGER      :: HDEV = 0!  unit number for holidays file
        INTEGER         LDEV    !  unit number for log file
        INTEGER      :: MDEV = 0!  unit number for mobile codes file
        INTEGER         PDEV    !  unit number for supplemental tmprl file
        INTEGER         RDEV    !  unit number for temporal profile file
        INTEGER         SDEV    !  unit number for ASCII inventory file
        INTEGER         TDEV    !  unit number for emission processes file
        INTEGER         XDEV    !  unit no. for cross-reference file

        CHARACTER*16 :: ANAME = ' '    !  logical name for ASCII inven input 
        CHARACTER*16 :: GNAME = ' '    !  ungridding matrix
        CHARACTER*16 :: DNAME = 'NONE' !  day-specific  input file, or "NONE"
        CHARACTER*16 :: ENAME = ' '    !  logical name for I/O API inven input
        CHARACTER*16 :: FNAME = ' '    !  emission factors file
        CHARACTER*16 :: HNAME = 'NONE' !  hour-specific input file, or "NONE"
        CHARACTER*16 :: TNAME = ' '    !  timestepped (low-level) output file
        CHARACTER*16 :: TMPNAME = ' '  !  temporary inventory logical name

C...........   Other local variables

        INTEGER         I, J, K, L, L1, L2, N, S, T

        INTEGER         IOS, IOS1, IOS2, IOS3, IOS4 ! i/o status
        INTEGER         IOS6, IOS7, IOS8, IOS9      ! i/o status
        INTEGER         AVERTYPE            ! time period averaging type
        INTEGER         DYSTPOS, DYENDPOS   ! start and end position in file name string
        INTEGER         EARLYDATE           ! earliest starting date based on time zones
        INTEGER         EARLYTIME           ! earliest starting time
        INTEGER         EDATE, ETIME        ! ending Julian date and time
        INTEGER         EFSDATE, EFEDATE    ! start and end date of current ef file
        INTEGER         ENLEN               ! length of ENAME string
        INTEGER         ENDPOS              ! ending position in ef day array
        INTEGER         FIRSTPOS            ! temporary position in file name string
        INTEGER         FDATE, FTIME        ! emission factor date and time
        INTEGER         HYPPOS              ! position of hyphen in file name string
        INTEGER         JDATE, JTIME        ! Julian date and time
        INTEGER         LATEDATE            ! latest ending date based on time zones
        INTEGER         LATETIME            ! latest ending time
        INTEGER         NDAYS               ! no. days in episode
        INTEGER         NINVARR             ! no. inventory variables to read
        INTEGER         NLINES              ! no. lines in ef list file
        INTEGER         NMATX               ! size of ungridding matrix
        INTEGER         NMAJOR              ! no. major sources
        INTEGER         NPING               ! no. ping sources
        INTEGER         NSTEPS              ! number of output time steps
        INTEGER      :: PYEAR = 0           ! projected year
        INTEGER         SDATE, STIME        ! starting Julian date and time
        INTEGER         STPOS               ! starting position in ef day array
        INTEGER         TNLEN               ! length of TNAME string
        INTEGER         TSTEP               ! output time step
        INTEGER         TZONE               ! output-file time zone
        INTEGER         TZMIN               ! minimum time zone in inventory      
        INTEGER         TZMAX               ! maximum time zone in inventory      

        REAL            RTMP                ! tmp float

        LOGICAL         DFLAG   !  true: day-specific  file available
        LOGICAL      :: EFLAG = .FALSE.  !  error-flag
        LOGICAL      :: EFLAG2= .FALSE.  !  error-flag (2)
        LOGICAL         ENDFLAG !  true: couldn't find file end date
        LOGICAL         HFLAG   !  true: hour-specific file available
        LOGICAL         MFLAG   !  true: mobile codes file available
        LOGICAL         NFLAG   !  true: use all uniform temporal profiles
        LOGICAL         WFLAG   !  true: write QA on current time step
        LOGICAL      :: USETIME( 4 ) = .FALSE. ! true: time period is used

        CHARACTER*8              TREFFMT   ! tmprl x-ref format (SOURCE|STANDARD)
        CHARACTER*14             DTBUF     ! buffer for MMDDYY
        CHARACTER*3              INTBUF    ! buffer for integer
        CHARACTER*20             MODELNAM  ! emission factor model name
        CHARACTER(LEN=256)       CURFNM    ! current emission factor file name
        CHARACTER(LEN=16)        CURLNM    ! current ef logical file name
        CHARACTER(LEN=IOVLEN3)   VOLNAM    ! volatile pollutant name
        CHARACTER*300            MESG      ! buffer for M3EXIT() messages
        CHARACTER(LEN=IOVLEN3)   CBUF      ! pollutant name temporary buffer 
        CHARACTER(LEN=IOVLEN3)   EBUF      ! pollutant name temporary buffer 
        CHARACTER(LEN=20)        SEARCHSTR ! string used in search
        CHARACTER(LEN=MXDLEN3)   TEMPLINE  ! line from file description

        CHARACTER*16 :: PROGNAME = 'TEMPORAL' ! program name

C***********************************************************************
C   begin body of program TEMPORAL

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Obtain settings from the environment...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Get environment variable that overrides temporal profiles and 
C               uses only uniform profiles.
        NFLAG = ENVYN( 'UNIFORM_TPROF_YN', MESG, .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get the name of the emission factor model to use for one run
        IF ( CATEGORY .EQ. 'MOBILE' ) THEN
            MESG = 'Emission factor model'
            CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE6', MODELNAM, IOS)
        ELSE
            MODELNAM = ' '
        END IF

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input files
C.........  Also, store source-category specific information in the MODINFO 
C           module.
        CALL OPENTMPIN( MODELNAM, NFLAG, ENAME, ANAME, DNAME, HNAME, 
     &                  GNAME, SDEV, XDEV, RDEV, CDEV, HDEV, TDEV, 
     &                  MDEV, EDEV, PYEAR )

C.........  Determine status of some files for program control purposes
        DFLAG = ( DNAME .NE. 'NONE' )  ! Day-specific emissions
        HFLAG = ( HNAME .NE. 'NONE' )  ! Hour-specific emissions
        MFLAG = ( MDEV .NE. 0 )        ! Use mobile codes file

C.........  Get length of inventory file name
        ENLEN = LEN_TRIM( ENAME )

C.........  Get episode settings from the Models-3 environment variables
        SDATE  = 0
        STIME  = 0
        NSTEPS = 1
        TSTEP  = 10000  
        CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C.........  Compare base year with episode and warn if not consistent
        IF( SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR, 
     &             'is inconsistent with year ' // CRLF() // BLANK10 //
     &             'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

C.........  Give a note if running for a projected year
        IF( PYEAR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//
     &             'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

C.........  Calculate the ending date and time
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

C.........  Check requested episode against available emission factors

C.........  For day-specific data input...
        IF( DFLAG ) THEN

C.............  Get header description of day-specific input file
            IF( .NOT. DESC3( DNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &                       'Could not get description of file "' 
     &                       // DNAME( 1:LEN_TRIM( DNAME ) ) // '"', 2 )
            END IF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( DYPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPNAM', PROGNAME )
            ALLOCATE( DYPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPDSC', PROGNAME )
            DYPNAM = ' '  ! array
            DYPDSC = ' '  ! array

C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( DNAME, SDATE, STIME, EDATE, ETIME, TZONE,  
     &                    NIPPA, EANAM, NDYPOA, NDYSRC, EFLAG, DYPNAM,
     &                    DYPDSC )

        ENDIF

C.........  Allocate memory for reading day-specific emissions data
C.........  NDYSRC is initialized to zero in case DFLAG is false
        ALLOCATE( INDXD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXD', PROGNAME )
        ALLOCATE( EMACD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACD', PROGNAME )

C.........  For hour-specific data input...
        IF( HFLAG ) THEN

C............. Get header description of hour-specific input file
            IF( .NOT. DESC3( HNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &                       'Could not get description of file "' 
     &                       // HNAME( 1:LEN_TRIM( HNAME ) ) // '"', 2 )
            ENDIF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( HRPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPNAM', PROGNAME )
            ALLOCATE( HRPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPDSC', PROGNAME )
            HRPNAM = ' '  ! array
            HRPDSC = ' '  ! array

C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( HNAME, SDATE, STIME, EDATE, ETIME, TZONE,  
     &                    NIPPA, EANAM, NHRPOA, NHRSRC, EFLAG2, HRPNAM,
     &                    HRPDSC )

        ENDIF

C.........  Allocate memory for reading hour-specific emissions data
C.........  NHRSRC is initialized to 0 in case HFLAG is false
        ALLOCATE( INDXH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXH', PROGNAME )
        ALLOCATE( EMACH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACH', PROGNAME )

        IF( EFLAG .OR. EFLAG2 ) THEN
            MESG = 'Problem with day- or hour-specific inputs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'TZONES'
        IVARNAMS( 3 ) = 'TPFLAG'
        IVARNAMS( 4 ) = 'CSCC'
        IVARNAMS( 5 ) = 'CSOURC'

C.........  Set inventory variables to read for specific source categories
        IF( CATEGORY .EQ. 'AREA' ) THEN
            NINVARR = 5

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
            NINVARR = 9
            IVARNAMS( 6 ) = 'IRCLAS'
            IVARNAMS( 7 ) = 'IVTYPE'
            IVARNAMS( 8 ) = 'CLINK'
            IVARNAMS( 9 ) = 'CVTYPE'

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
            NINVARR = 5

        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Reset TPFLAG if ozone-season emissions are being used since
C           we don't want to apply the monthly adjustment factors in this case.
        IF ( INVPIDX .EQ. 1 ) THEN
            DO S = 1, NSRC
                IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN
                    TPFLAG( S ) = TPFLAG( S ) / MTPRFAC
                END IF
            END DO
        END IF

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

C.........  Adjust TZMIN for possibility of daylight savings
        TZMIN = MAX( TZMIN - 1, 0 )

C.........  Read special files...

C.........  Read region codes file
        CALL RDSTCY( CDEV, NINVIFIP, INVIFIP )

C.........  Populate filter for sources that use daylight time
        CALL SETDAYLT

C.........  Read holidays file
        CALL RDHDAYS( HDEV, SDATE, EDATE )

C.........  When mobile codes file is being used read mobile codes file
        IF( MFLAG ) CALL RDMVINFO( MDEV )

C.........  Perform steps needed for using activities and emission factors

        IF( NIACT .GT. 0 ) THEN

C.............  Allocate memory for emission factor arrays
            ALLOCATE( TEMPEF( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TEMPEF', PROGNAME )
            ALLOCATE( EFTYPE( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFFILE', PROGNAME )
            ALLOCATE( EFIDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFIDX', PROGNAME )

            TEMPEF = 0.
            EFTYPE = ' '
            EFIDX  = 0

C.............  Determine number of days in episode

C.............  Earliest day is start time in maximum time zone
            EARLYDATE = SDATE
            EARLYTIME = STIME
            CALL NEXTIME( EARLYDATE, EARLYTIME, 
     &                   -( TZMAX - TZONE )*10000 )
            
C.............  If time is before 6 am, need previous day also
            IF( EARLYTIME < 60000 ) THEN
                EARLYDATE = EARLYDATE - 1
            END IF
            
C.............  Latest day is end time in minimum time zone
            LATEDATE = EDATE
            LATETIME = ETIME
            CALL NEXTIME( LATEDATE, LATETIME, 
     &                   -( TZMIN - TZONE )*10000 )

C.............  If time is before 6 am, don't need last day
            IF( LATETIME < 60000 ) THEN
                LATEDATE = LATEDATE - 1
            END IF

            NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
            NDAYS = NDAYS + 1
            ALLOCATE( EFDAYS( NDAYS,4 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFDAYS', PROGNAME )
            EFDAYS = 0

C.............  Read header of ungridding matrix
            IF( .NOT. DESC3( GNAME ) ) THEN
            	MESG = 'Could not get description for file ' // GNAME
            	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Store number of ungridding factors
            NMATX = NCOLS3D

C.............  Allocate memory for ungridding matrix
            ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UMAT', PROGNAME )

C.............  Read ungridding matrix
            CALL RDUMAT( GNAME, NSRC, NMATX, NMATX, UMAT( 1 ),
     &                   UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 )  )

C.............  Store sources that are outside the grid
            DO S = 1, NSRC
                IF( UMAT( S ) == 0 ) THEN
                    EFIDX( S ) = -9
                END IF
            END DO

C.............  Read emission processes file.  Populate array in MODEMFAC.
            CALL RDEPROC( TDEV )

C.............  Loop through activities and...
C.............  NOTE - this is not fully implemented for multiple activities. 
C               To do this, the data structures and RDEFACS will need to be 
C               updated. Also, the variable names in the emission factor file
C               are not truly supporting 16-character pollutant and 
C               emission process names, because it is only set up for MOBILE5
            DO I = 1, NIACT

C.................  Skip activities that do not have emissions types
                IF( NETYPE( I ) .LE. 0 ) CYCLE            

C.................  Set up emission process variable names
                CALL EFSETUP( 'NONE', MODELNAM, MXVARS3, NEFS, VNAME3D, 
     &                         UNITS3D, VDESC3D, VOLNAM )

            END DO

C.............  Read list of emission factor files
            NLINES = GETFLINE( EDEV, 'Emission factor file list' )

            ALLOCATE( EFLIST( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFLIST', PROGNAME )
            ALLOCATE( EFLOGS( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFLOGS', PROGNAME )
        
            EFLIST = ' '
            EFLOGS = ' '
            CALL RDLINES( EDEV, 'Emission factor file list', NLINES, 
     &                    EFLIST )

C.............  Loop through EF files
            MESG = 'Checking emission factor files...'
            CALL M3MSG2( MESG )

            DO N = 1, NLINES

                CURFNM = EFLIST( N )
                
C.................  Skip any blank lines
                IF( CURFNM == ' ' ) CYCLE

C.................  Determine file type
                IF( INDEX( CURFNM, 'daily' ) > 0 ) THEN
                    AVERTYPE = DAILY
                ELSE IF( INDEX( CURFNM, 'weekly' ) > 0 ) THEN
                    AVERTYPE = WEEKLY
                ELSE IF( INDEX( CURFNM, 'monthly' ) > 0 ) THEN
                    AVERTYPE = MONTHLY
                ELSE IF( INDEX( CURFNM, 'episode' ) > 0 ) THEN
                    AVERTYPE = EPISLEN
                ELSE
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not determine time period ' //
     &                     'of file ' // CURFNM( 1:LEN_TRIM( CURFNM ) )
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
C.................  Assign and store logical file name
                WRITE( INTBUF,94030 ) N
                CURLNM = 'EMISFAC_' // ADJUSTL( INTBUF )
                EFLOGS( N ) = CURLNM

C.................  Set logical file name
                IF( .NOT. SETENVVAR( CURLNM, CURFNM ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not set logical file name ' //
     &                     'for file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                USETIME( AVERTYPE ) = .TRUE.

C.................  Try to open file   
                IF( .NOT. OPEN3( CURLNM, FSREAD3, PROGNAME ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open emission factors ' //
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Read file description
                IF( .NOT. DESC3( CURLNM ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not get description for ' // 
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                EFSDATE = SDATE3D

C.................  Find end date in file description
                SEARCHSTR = '/END DATE/ '
                L = LEN_TRIM( SEARCHSTR ) + 1
                ENDFLAG = .FALSE.
                
                DO I = 1, MXDESC3
                   IF( INDEX( FDESC3D( I ), 
     &                        SEARCHSTR( 1:L ) ) > 0 ) THEN
     	                   TEMPLINE = FDESC3D( I )
                   	   IF( CHKINT( TEMPLINE( L+1:L+8 ) ) ) THEN
                               EFEDATE = STR2INT( TEMPLINE( L+1:L+8 ) )
                           EXIT
                       ELSE
                           ENDFLAG = .TRUE.
                           EXIT
                       END IF
                   END IF
                   
                   IF( I == MXDESC3 ) THEN
                       ENDFLAG = .TRUE.
                   END IF
                END DO
                	
                IF( ENDFLAG ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not get ending date of ' //
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Determine starting and ending positions in array
                STPOS  = SECSDIFF( EARLYDATE, 0, EFSDATE, 0 )/(24*3600)
                STPOS  = STPOS + 1
                ENDPOS = SECSDIFF( EARLYDATE, 0, EFEDATE, 0 )/(24*3600)
                ENDPOS = ENDPOS + 1

C.................  Make sure starting and ending positions are valid
                IF( STPOS < 1 ) THEN
                    IF( ENDPOS > 0 ) THEN
                        STPOS = 1
                    ELSE
                    	CYCLE
                    END IF
                END IF 
                
                IF( ENDPOS > NDAYS ) THEN
                    IF( STPOS <= NDAYS ) THEN
                        ENDPOS = NDAYS
                    ELSE
                        CYCLE
                    END IF
                END IF

C.................  Store day info
                DO I = STPOS, ENDPOS
                    EFDAYS( I, AVERTYPE ) = N
                END DO

C.................  Allocate memory for temporary source info
                IF( ALLOCATED( SRCS ) ) DEALLOCATE( SRCS )
                ALLOCATE( SRCS( NROWS3D ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SRCS', PROGNAME )
            
                SRCS = 0
            
C.................  Read source information
                IF( .NOT. READ3( CURLNM, 'SOURCES', ALLAYS3, 
     &                           SDATE3D, STIME3D, SRCS ) ) THEN
     	            EFLAG = .TRUE.
                    MESG = 'ERROR: Could not read SOURCES ' // 
     &                     'from file ' // CRLF() // BLANK10 // '"' //
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Store source information
                DO S = 1, NROWS3D
                                    
                    IF( SRCS( S ) /= 0 ) THEN

C.........................  Make sure source number is valid
                        IF( SRCS( S ) < 1 .OR. SRCS( S ) > NSRC ) CYCLE
                        
C.........................  Skip sources that are outside the grid
                        IF( EFIDX( SRCS( S ) ) == -9 ) CYCLE
                    
                        EFIDX( SRCS( S ) ) = S
                        WRITE( EFTYPE( SRCS( S ) ),'(I1)' ) AVERTYPE
                    END IF
                END DO
                
C.................  Close current file
                IF( .NOT. CLOSE3( CURLNM ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not close file ' // 
     &                     CURFNM( 1:LEN_TRIM( CURFNM ) )
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
            END DO
            	
C.............  Exit if there was a problem with the emission factor files
            IF( EFLAG ) THEN
                MESG = 'Problem checking emission factor files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Make sure all days are covered
            DO I = DAILY, EPISLEN
                IF( USETIME( I ) .EQV. .TRUE. ) THEN 
                    IF( MINVAL( EFDAYS( 1:NDAYS,I ) ) == 0 ) THEN
                        MESG = 'ERROR: Emission factor files do not ' //
     &                         'cover requested time period.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF
            END DO

C............  Print warning for sources that don't have emission factors
           DO I = 1, NSRC
               IF( EFIDX( I ) == 0 ) THEN
                   WRITE( MESG,94070 ) 'WARNING: No VMT or emission ' //
     &                    'factors available for' // CRLF() // 
     &                    BLANK10 // 'Region: ', IFIP( I ),
     &                    ' SCC: ' // CSCC( I )
                   CALL M3MESG( MESG )
                   EFIDX( I ) = -1
               END IF
           END DO

        END IF

C.........  Determine all of the variables to be output based on the 
C           activities and input pollutants.  
C.........  NOTE - Uses NETYPE, EMTACTV, and EMTNAM, and sets units, and units
C           conversion factors for all pollutants and activities
        CALL TMNAMUNT
        
C.........  Reset the number of all output variables as the number of pollutants
C           and emission types, instead of the number of pollutants and 
C           activities
        NIPPA = NIPOL
        DO I = 1, NIACT
            NIPPA = NIPPA + NETYPE( I )
        END DO

C.........  Allocate memory for I/O pol names, activities, & emission types
C.........  Will be resetting EANAM to include the emission types instead
C           of the activities
        DEALLOCATE( EANAM )
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( ALLIN( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN', PROGNAME )

C.........  Create 1-d arrays of I/O pol names
C.........  If pollutant is created by one or more activities, then give a
C           warning.
        N = 0
        DO I = 1, NIPOL

C.............  Look for pollutant in list of pollutants created by activities
            J = INDEX1( EINAM( I ), NEPOL, EMTPOL )

C.............  If pollutant created by activity, skip from this list, unless
C               pollutant is also part of the inventory pollutants
            IF( J .GT. 0 ) THEN
                L1 = LEN_TRIM( EINAM ( I ) )
                L2 = LEN_TRIM( EAREAD( I ) )
                MESG = 'WARNING: Pollutant "' // EINAM(I)( 1:L1 ) //
     &                 '" is explicitly in the inventory and' //
     &                 CRLF() // BLANK10 // 'it is also generated by '
     &                 // 'activity data.'
                CALL M3MSG2( MESG )
            END IF

            N = N + 1
            ALLIN( N ) = EAREAD( I )
            EANAM( N ) = EINAM ( I )

        END DO

C.........  Add activities, & emission types to read and output lists
        J = NIPOL
        DO I = 1, NIACT

            K = NETYPE( I )  ! Number of emission types

C.............  If any emissions types associated with this activity, store them
            IF ( K .GT. 0 ) THEN
                ALLIN( J+1:J+K ) = ACTVTY( I )
                EANAM( N+1:N+K ) = EMTNAM( 1:K, I )
                N = N + K
            END IF
            J = J + K

        END DO

C.........  Reset number of pollutants and emission types based on those used
        NIPPA = N

C.........  Set up and open I/O API output file(s) ...
        CALL OPENTMP( ENAME, SDATE, STIME, TSTEP, TZONE, NPELV,
     &                TNAME, PDEV )

        TNLEN = LEN_TRIM( TNAME )

C.........  Read temporal-profile cross-reference file and put into tables
C.........  Only read entries for pollutants that are in the inventory.
C.........  Only read if not using uniform temporal profiles
        IF( .NOT. NFLAG ) CALL RDTREF( XDEV, TREFFMT )

C.........  Read temporal-profiles file:  4 parts (monthly, weekly, 
C           weekday diurnal, and weekend diurnal)
        CALL M3MSG2( 'Reading temporal profiles file...' )

        NMON = RDTPROF( RDEV, 'MONTHLY', NFLAG )
        NWEK = RDTPROF( RDEV, 'WEEKLY' , NFLAG )
        NHRL = RDTPROF( RDEV, 'DIURNAL', NFLAG )

C.........  Adjust temporal profiles for use in generating temporal emissions
C.........  NOTE - All variables are passed by modules.
        CALL NORMTPRO

C.........  It is important that all major arrays must be allocated by this 
C           point because the next memory allocation step is going to pick a
C           data structure that will fit within the limits of the host.

C.........  Allocate memory, but allow flexibility in memory allocation
C           for second dimension.
C.........  The second dimension (the number of pollutants and emission types)
C            can be different depending on the memory available.
C.........  To determine the approproate size, first attempt to allocate memory
C           for all pollutants & emission types to start, and if this fails,
C           divide pollutants into even groups and try again.

        NGSZ = NIPPA            ! No. of pollutant & emis types in each group
        NGRP = 1               ! Number of groups
        DO

            ALLOCATE( TMAT ( NSRC, NGSZ, 24 ), STAT=IOS1 )
            ALLOCATE( MDEX ( NSRC, NGSZ )    , STAT=IOS2 )
            ALLOCATE( WDEX ( NSRC, NGSZ )    , STAT=IOS3 )
            ALLOCATE( DDEX ( NSRC, NGSZ )    , STAT=IOS4 )
            ALLOCATE( EMAC ( NSRC, NGSZ )    , STAT=IOS6 )
            ALLOCATE( EMACV( NSRC, NGSZ )    , STAT=IOS7 )
            ALLOCATE( EMIST( NSRC, NGSZ )    , STAT=IOS8 )
            
C.............  Only need to allocate if we have activities            
            IF( NIACT .GT. 0 ) THEN
                ALLOCATE( EMFAC( NSRC, NGSZ )    , STAT=IOS9 )
            ELSE
            	IOS9 = 0
            END IF
            
            IF( IOS1 .GT. 0 .OR. IOS2 .GT. 0 .OR. IOS3 .GT. 0 .OR.
     &          IOS4 .GT. 0 .OR. IOS6 .GT. 0 .OR.
     &          IOS7 .GT. 0 .OR. IOS8 .GT. 0 .OR. IOS9 .GT. 0 ) THEN

                IF( NGSZ .EQ. 1 ) THEN
                    J = 8 * NSRC * 31    ! Assume 8-byte reals
                    WRITE( MESG,94010 ) 
     &                'Insufficient memory to run program.' //
     &                CRLF() // BLANK5 // 'Could not allocate ' // 
     &                'pollutant-dependent block of', J, 'bytes.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NGRP = NGRP + 1
                NGSZ = NGSZ / NGRP 
                NGSZ = NGSZ + ( NIPPA - NGSZ * NGRP )

                IF( ALLOCATED( TMAT  ) ) DEALLOCATE( TMAT )
                IF( ALLOCATED( MDEX  ) ) DEALLOCATE( MDEX )
                IF( ALLOCATED( WDEX  ) ) DEALLOCATE( WDEX )
                IF( ALLOCATED( DDEX  ) ) DEALLOCATE( DDEX )
                IF( ALLOCATED( EMAC  ) ) DEALLOCATE( EMAC )
                IF( ALLOCATED( EMACV ) ) DEALLOCATE( EMACV )
                IF( ALLOCATED( EMIST ) ) DEALLOCATE( EMIST )
                IF( ALLOCATED( EMFAC ) ) DEALLOCATE( EMFAC )

            ELSE
                EXIT

            END IF

        END DO
        
C.........  Allocate a few small arrays based on the size of the groups
C.........  NOTE that this has a small potential for a problem if these little
C           arrays exceed the total memory limit.
        ALLOCATE( ALLIN2D( NGSZ, NGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN2D', PROGNAME )
        ALLOCATE( EANAM2D( NGSZ, NGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM2D', PROGNAME )
        ALLOCATE( LDSPOA( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDSPOA', PROGNAME )
        ALLOCATE( LHSPOA( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LHSPOA', PROGNAME )
        ALLOCATE( LHPROF( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LHPROF', PROGNAME )

C.........  Create 2-d arrays of I/O pol names, activities, & emission types 
        ALLIN2D = ' '
        EANAM2D = ' '
        I = 0
        DO N = 1, NGRP 
            DO J = 1, NGSZ
                I = I + 1
                IF ( I .LE. NIPPA ) THEN
                    EANAM2D( J,N ) = EANAM( I )
                    ALLIN2D( J,N ) = ALLIN( I )
                END IF
            END DO
        END DO

C.........  Loop through pollutant/emission-type groups
        DO N = 1, NGRP

C.............  Skip group if the first pollutant in group is blank (this
C               shouldn't happen, but it is happening, and it's easier to
C               make this fix).
            IF ( EANAM2D( 1,N ) .EQ. ' ' ) CYCLE

C.............  Write message stating the pols/emission-types being processed
            CALL POLMESG( NGSZ, EANAM2D( 1,N ) )

C.............  Set up logical arrays that indicate which pollutants/activities
C               are day-specific and which are hour-specific.
C.............  Also set flag for which hour-specific pollutants/activities
C               are actually diurnal profiles instead of emissions
            LDSPOA = .FALSE.   ! array
            DO I = 1, NDYPOA
                J = INDEX1( DYPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                LDSPOA( J ) = .TRUE.
            END DO

            LHSPOA = .FALSE.   ! array
            LHPROF = .FALSE.   ! array
            DO I = 1, NHRPOA
                J = INDEX1( HRPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                LHSPOA( J ) = .TRUE.

                CALL UPCASE( HRPDSC( I ) )
                K = INDEX( HRPDSC( I ), 'PROFILE' )
                IF( K .GT. 0 ) LHPROF( J ) = .TRUE.
            END DO

C.............  Initialize emissions, activities, and other arrays for this
C               pollutant/emission-type group
            TMAT  = 0.
            MDEX  = IMISS3
            WDEX  = IMISS3
            DDEX  = IMISS3
            EMAC  = 0.
            EMACV = 0.
            EMIST = 0.
            IF( NIACT .GT. 0 ) EMFAC = IMISS3

C.............  Assign temporal profiles by source and pollutant
            CALL M3MSG2( 'Assigning temporal profiles to sources...' )

C.............  If using uniform profiles, set all temporal profile number
C               to 1; otherwise, assign profiles with cross-reference info
            IF( NFLAG ) THEN
                MDEX = 1
        	WDEX = 1
        	DDEX = 1

            ELSE
                CALL ASGNTPRO( NGSZ, EANAM2D( 1,N ), TREFFMT )

            END IF

C.............  Read in pollutant emissions or activities from inventory for 
C               current group
            DO I = 1, NGSZ

                EBUF = EANAM2D( I,N )
                CBUF = ALLIN2D( I,N )
                L1   = LEN_TRIM( CBUF )

C.................  Skip blanks that can occur when NGRP > 1
                IF ( CBUF .EQ. ' ' ) CYCLE

C.................  Read the emissions data in either map format
C                   or the old format.
                CALL RDMAPPOL( ENAME, NMAP, MAPNAM, MAPFIL, NSRC,
     &                         1, 1, EBUF, CBUF, 1, EMAC( 1,I )   )

C...............  If there are any missing values in the data, give an
C                 error to avoid problems in genhemis routine
                RTMP = MINVAL( EMAC( 1:NSRC,I ) )
                IF( RTMP .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Missing or negative emission '//
     &                     'value(s) in inventory for "' // 
     &                     CBUF( 1:L1 ) // '".'
                    CALL M3MSG2( MESG )
                END IF

C.................  If pollutant name is ozone-season-based, remove the
C                   prefix from the input pollutant name
                K = INDEX1( CBUF, NIPPA, EAREAD )
                J = INDEX( CBUF, OZNSEART )
                IF( J .GT. 0 ) THEN
                    CBUF = CBUF( CPRTLEN3+1:L1 )
                    ALLIN2D( I,N ) = CBUF
                    EAREAD ( K )   = CBUF
                END IF

            END DO

C.............  Abort if error found
            IF( EFLAG ) THEN
                MESG = 'Problem with input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  For each time step and pollutant or emission type in current 
C               group, generate hourly emissions, and write layer-1 emissions
C               file (or all data).
            JDATE = SDATE
            JTIME = STIME

C.............  Write supplemental temporal profiles file
            CALL WRTSUP( PDEV, NSRC, NGSZ, EANAM2D( 1,N ) )

C.............  Loop through time steps for current pollutant group
            DO T = 1, NSTEPS

                IF( NIACT .GT. 0 ) THEN

C.....................  Create array of emission factors by source

C.....................  Loop through pollutants/emission-types in this group
                    DO I = 1, NGSZ
                        CBUF = EANAM2D( I,N )
                        L1   = LEN_TRIM( CBUF )

C.........................  Skip blanks that can occur when NGRP > 1
                        IF ( CBUF .EQ. ' ' ) CYCLE

C.........................  Check that this pollutant uses emission factors
C                           Look for double underscore in pollutant name
                        K = INDEX( CBUF, ETJOIN )
                        IF( K == 0 ) THEN
                            EMFAC( :,I ) = -1
                            CYCLE
                        END IF 

C.........................  Loop through time zones
                        DO J = TZMIN, TZMAX
                   
C.............................  Adjust time zone J based on output time zone and account
C                               for 6 AM starting time in files
                            K = J - TZONE + 6
                   
                            FDATE = JDATE
                            FTIME = JTIME
                            CALL NEXTIME( FDATE, FTIME, -K * 10000 )

C.............................  Use date and time to find appropriate ef file
                            STPOS = SECSDIFF( EARLYDATE, 0, FDATE, 0 )
                            STPOS = STPOS / ( 24*3600 )
                            STPOS = STPOS + 1
                            
                            DO L = DAILY, EPISLEN
                                IF( USETIME( L ) .EQV. .FALSE. ) CYCLE
                                
                                IF( STPOS <= 0 .OR. STPOS > NDAYS ) THEN
                                    WRITE( *,* ) NDAYS
                                    MESG = 'ERROR: Invalid position'
                                    CALL M3EXIT( PROGNAME, FDATE, FTIME,
     &                                           MESG, 2 )
                                END IF
                                
                                CURFNM = EFLIST( EFDAYS( STPOS,L ) )
                                CURLNM = EFLOGS( EFDAYS( STPOS,L ) )

C.................................  Set logical file name
                                IF( .NOT. SETENVVAR( CURLNM, 
     &                                               CURFNM ) ) THEN
                                    EFLAG = .TRUE.
                                    MESG = 'ERROR: Could not set ' //
     &                                     'logical file name for ' //
     &                                     'file ' // CRLF() // BLANK10
     &                                     // '"' // TRIM( CURFNM ) // 
     &                                     '".'
                                    CALL M3EXIT( PROGNAME, FDATE, FTIME,
     &                                           MESG, 2 )
                                END IF

C.................................  Open current file
                                IF( .NOT. OPEN3( CURLNM, FSREAD3, 
     &                                           PROGNAME ) ) THEN
     	                            EFLAG = .TRUE.
                                    MESG = 'ERROR: Could not open ' //
     &                                     'emission factors file ' //
     &                                     CRLF() // BLANK10 // '"' //
     &                                     CURFNM( 1:LEN_TRIM(CURFNM) )
     &                                     // '".'
                                    CALL M3EXIT( PROGNAME, FDATE, FTIME,
     &                                           MESG, 2 )
                                END IF

C.................................  Read file description
                                IF( .NOT. DESC3( CURLNM ) ) THEN
                                    MESG = 'ERROR: Could not get ' //
     &                                     'description for file ' //
     &                                     CRLF() // BLANK10 // '"' // 
     &                                     CURFNM( 1:LEN_TRIM(CURFNM) )
     &                                     // '".'
                                    CALL M3EXIT( PROGNAME, FDATE, FTIME, 
     &                                           MESG, 2 )
                                END IF

C.................................  Read emission factors from current file
                                IF( .NOT. READ3( CURLNM, CBUF, ALLAYS3, 
     &                              SDATE3D, FTIME, TEMPEF ) ) THEN
                                    EFLAG = .TRUE.
                                    MESG = 'Error reading "'// 
     &                                     CBUF(1:L1) //
     &                                     '" from file ' // 
     &                                     CRLF() // BLANK10 // '"' // 
     &                                     CURFNM( 1:LEN_TRIM(CURFNM) ) 
     &                                     // '."'
                                    CALL M3EXIT( PROGNAME, FDATE, FTIME,
     &                                           MESG, 2 )
                                END IF
                
C.................................  Store emission factors by source                            
                                DO S = 1, NSRC
                            
C.....................................  Skip sources that are outside the grid or
C                                       don't use emission factors
                                    IF( EFIDX( S ) == -9 .OR. 
     &                                  EFIDX( S ) == -1 ) THEN
                                        EMFAC( S,I ) = 0
                                        CYCLE
                                    END IF

                                    IF( TZONES( S ) == J .AND. 
     &                                  STR2INT(EFTYPE( S )) == L ) THEN
                                        EMFAC( S,I ) = 
     &                                       TEMPEF( EFIDX( S ) )
                                    END IF
                                END DO   ! end source loop

                            END DO   ! end time period loop

                        END DO   ! end time zone loop

C.........................  If there are any missing values in the data, give an
C                           error to avoid problems in genhemis routine
                        RTMP = MINVAL( EMFAC( 1:NSRC,I ) )
                        IF( RTMP == IMISS3 ) THEN
                            EFLAG = .TRUE.
                            MESG = 'ERROR: Missing emission ' //
     &                         'factors(s) for "'// CBUF( 1:L1 ) // '".'
                            CALL M3MSG2( MESG )
                        END IF
                    END DO  ! End loop on pollutants/emission-types I in this group

C.....................  Abort if error found
                    IF( EFLAG ) THEN
                        MESG = 'Problem with emission factors.'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
 
                END IF

C.................  Generate hourly emissions for current hour
                CALL GENHEMIS( NGSZ, JDATE, JTIME, TZONE, DNAME, HNAME, 
     &                         ALLIN2D( 1,N ), EANAM2D( 1,N ), 
     &                         EMAC, EMFAC, EMACV, TMAT, EMIST )

C.................  Loop through pollutants/emission-types in this group
                DO I = 1, NGSZ

                    CBUF = EANAM2D( I,N )

C.....................  Skip blanks that can occur when NGRP > 1
                    IF ( CBUF .EQ. ' ' ) CYCLE

C.....................  Write hourly emissions to I/O API NetCDF file
                    IF( .NOT. WRITE3( TNAME, CBUF, JDATE, JTIME, 
     &                                EMIST( 1,I )              ) ) THEN

                        L = LEN_TRIM( CBUF )
                        MESG = 'Could not write "' // CBUF( 1:L ) // 
     &                         '" to file "' // TNAME( 1:TNLEN ) // '."'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF

                END DO  ! End loop on pollutants/emission-types I in this group

C.................  Advance the output date/time by one time step
                CALL NEXTIME( JDATE, JTIME, TSTEP )

C.................  Call QA report routine
c               WFLAG = ( T .EQ. NSTEPS )
c               CALL QATMPR( LDEV, NGSZ, T, JDATE, JTIME, WFLAG, 
c    &                       EANAM2D( 1,N ), EMAC )

            END DO      ! End loop on time steps T

        END DO          ! End loop on pollutant groups N

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( I3 )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )
     
94070   FORMAT( A, I5, A )

        END PROGRAM TEMPORAL

