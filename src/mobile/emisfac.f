
        PROGRAM EMISFAC

C***********************************************************************
C  subroutine body starts at line 153
C
C  DESCRIPTION:
C       Reads information from SPDSUM and hourly meteorology files to
C       creates MOBILE6 input files. Calls modified version of MOBILE6 
C       to generate emission factors. Matches factors to sources and 
C       writes to output files. Loops over all days for a given meteorology
C       averaging group, stopping when no more meteorology files are available.
C
C  PRECONDITIONS REQUIRED:
C       MBSETUP and PREMOBL must have been run.
C
C  SUBROUTINES AND FUNCTIONS CALLED: 
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************        
        
C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NSRC, NIACT

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDVTS
        
C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, QVHOUR, BPHOUR, BPDAY, RHHOUR

C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS, NUMSCEN, SCENLIST, EMISSIONS,
     &                      INPUTHC, OUTPUTHC, NSUBPOL, SUBPOLS,
     &                      NEPOL, EMTPOL, EMTNAM

C.........This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         CVTRDTYPE
        INTEGER         CVTVEHTYPE
        CHARACTER(2)    CRLF
        LOGICAL         SETENVVAR
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        LOGICAL         CHKINT
        INTEGER         JUNIT

        EXTERNAL        PROMPTFFILE, PROMPTMFILE, ENVYN, GETFLINE, 
     &                  CVTRDTYPE, CVTVEHTYPE, CRLF, SETENVVAR
     &                  SECSDIFF, STR2INT, CHKINT, JUNIT

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: GRPLIST( :,: ) ! contents of GROUP file
        INTEGER, ALLOCATABLE :: TEMPCTY( : )   ! list of counties from temperature file
        CHARACTER(IOVLEN3), ALLOCATABLE :: RAWSUBS( : )  ! list of pollutants to be
                                                             ! subtracted from HC

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )
        
C.........  Unit numbers and logical file names
        INTEGER         CDEV     ! unit number for county MOBILE6 scenarios file (M6LIST)
        INTEGER         GDEV     ! unit number for time period group file (GROUP)
        INTEGER         IDEV     ! tmp unit number if ENAME is map file
        INTEGER         LDEV     ! unit number for log file
        INTEGER         MDEV     ! unit number for concatenated MOBILE6 input file (M6INPUT)
        INTEGER         PDEV     ! unit number for speeds summary file (SPDSUM)
        INTEGER         SDEV     ! unit number for ASCII inventory file
        INTEGER         TDEV     ! unit number for emission processes file
        INTEGER         VDEV     ! unit number for inventory data table
        INTEGER         ZDEV     ! unit number for speed profiles file
        
        CHARACTER(16)   ANAME   !  logical name for ASCII inventory file
        CHARACTER(16)   ENAME   !  logical name for I/O API inventory file   
        CHARACTER(16)   FNAME   !  logical name for I/O API emission factors file
        CHARACTER(16)   INAME   !  tmp name for inven file of unknown fmt
        CHARACTER(16)   TNAME   !  logical name for I/O API temperature file

C.........   Other local variables
        INTEGER    I, J, K, L, L2! counters and indices
        INTEGER    IOS           ! i/o status
        INTEGER    NINVARR       ! number inventory variables to input
        INTEGER    SDATE         ! current start date
        INTEGER    STIME         ! current start time
        INTEGER    EDATE         ! current end date
        INTEGER    TZONE         ! time zone (not used)
        INTEGER    TSTEP         ! time step of input temperature data (HHMMSS)
        INTEGER    TIMEDIFF      ! difference between file date and current sdate
        INTEGER    TEMPDATE      ! temporary date to pass to read routine
        INTEGER    TEMPTIME      ! temporary time to pass to read routine
        INTEGER    NSTEPS        ! no. time steps in temperature data
        INTEGER    NROWS         ! no. grid rows   
        INTEGER    NGRPLINES     ! no. lines in GROUP file
        INTEGER    NUMSRC        ! total number of sources
        INTEGER :: MAXPOL = 0    ! max. number of pollutants
        INTEGER :: MAXVEH = 0    ! max. number of vehicle types
        INTEGER :: MAXFAC = 0    ! max. number of facility types

        LOGICAL :: EFLAG    = .FALSE.    ! error flag
        LOGICAL :: TEMPFLAG = .TRUE.     ! true: replace temperatures in M6 scenarios
        LOGICAL :: RHFLAG   = .TRUE.     ! true: replace humidity values in M6 scenarios
        LOGICAL :: INITIAL  = .TRUE.     ! true: first time through loop
        LOGICAL :: FEXIST   = .FALSE.    ! true: file exists
        LOGICAL :: FILEOPEN = .FALSE.    ! true: emisfacs file is open
        LOGICAL :: FNDINPUT = .FALSE.    ! true: found input hydrocarbon
        LOGICAL :: FNDOUTPUT = .FALSE.   ! true: found output hydrocarbon
        LOGICAL :: SPDFLAG  = .TRUE.     ! true: use speed profiles
        
        CHARACTER(20)      MODELNAM  ! emission factor model name
        CHARACTER(20)      GRP_NAME  ! temperature aggregation group
        CHARACTER(IOVLEN3) VOLNAM    ! volatile pollutant name
        CHARACTER(280)     M6INPUT   ! Mobile6 input file name
        CHARACTER(200)     TEMPDIR   ! location of hourly temperature files
        CHARACTER(256)     TEMPNAME  ! full temperature file name
        CHARACTER(200)     M6DIR     ! location of MOBILE6 files
        CHARACTER(200)     EMISDIR   ! directory for output EF files
        CHARACTER(20)      SEARCHSTR ! string used in search
        CHARACTER(MXDLEN3) TEMPLINE  ! line from file description
        CHARACTER(300)     MESG      ! message buffer 

        CHARACTER(16) :: PROGNAME = 'EMISFAC' ! program name
        
C***********************************************************************
C   begin body of program EMISFAC

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  End program if source category is not mobile sources
        IF( CATEGORY /= 'MOBILE' ) THEN
            L = LEN_TRIM( PROGNAME )
            MESG = 'Program ' // PROGNAME( 1:L ) // ' does not ' //
     &             'support ' // TRIM( CATEGORY ) // ' sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Obtain settings from the environment...

C.........  Get the name of the emission factor model to use
ccs        MESG = 'Emission factor model'
ccs        CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE6', MODELNAM, IOS )
        MODELNAM = 'MOBILE6'

C.........  Get environment variables that control program behavior
        TEMPFLAG = ENVYN( 'REPLACE_TEMPERATURES', 
     &                 'Replace temperatures in MOBILE6 scenarios',
     &                 .TRUE., IOS )
        RHFLAG = ENVYN( 'REPLACE_HUMIDITY', 'Replace humidity ' //
     &                  'values in MOBILE6 scenarios', .TRUE., IOS )
     
C.........  Get meteorology aggregation length
        MESG = 'Meteorology aggregation group'
        CALL ENVSTR( 'GROUP_TYPE', MESG, 'daily', GRP_NAME, IOS )
     
C.........  Check if speed profiles are to be used
        SPDFLAG = ENVYN( 'USE_SPEED_PROFILES', 
     &            'Use speed profiles instead of inventory speeds', 
     &            .FALSE., IOS )
     
C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open inventory file 
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

        PDEV = PROMPTFFILE(
     &           'Enter logical name for SPDSUM speed summary file',
     &           .TRUE., .TRUE., 'SPDSUM', PROGNAME )
        
        CDEV = PROMPTFFILE(
     &           'Enter logical name for M6LIST scenarios file',
     &           .TRUE., .TRUE., 'M6LIST', PROGNAME )
        
        TDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION PROCESSES file',
     &           .TRUE., .TRUE., 'MEPROC', PROGNAME )

        VDEV = PROMPTFFILE( 
     &           'Enter logical name for INVENTORY DATA TABLE file',
     &           .TRUE., .TRUE., 'INVTABLE', PROGNAME )

        IF( SPDFLAG ) THEN
            ZDEV = PROMPTFFILE(
     &           'Enter logical name for SPDPRO speed profiles file',
     &           .TRUE., .TRUE., 'SPDPRO', PROGNAME )
        END IF
        
        SELECT CASE( GRP_NAME )
        CASE( 'daily' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for DAILYGROUP file',
     &           .TRUE., .TRUE., 'DAILYGROUP', PROGNAME )
        CASE( 'weekly' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for WEEKLYGROUP file',
     &           .TRUE., .TRUE., 'WEEKLYGROUP', PROGNAME )
        CASE( 'monthly' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for MONTHLYGROUP file',
     &           .TRUE., .TRUE., 'MONTHLYGROUP', PROGNAME )
        CASE( 'episode' )
           GDEV = PROMPTFFILE(
     &           'Enter logical name for EPISODEGROUP file',
     &           .TRUE., .TRUE., 'EPISODEGROUP', PROGNAME )
        CASE DEFAULT
            MESG = 'ERROR: Unrecognized value for GROUP_TYPE.' //
     &             CRLF() // BLANK10 // 'Valid values are daily, ' //
     &             'weekly, monthly, or episode.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT
     
C.........  Get meteorology file directory from the environment
        MESG = 'Location of hourly meteorology files'
        CALL ENVSTR( 'SMK_METPATH', MESG, '.', TEMPDIR, IOS )
     
        TNAME = PROMPTMFILE(
     &          'Enter logical name for first hourly meteorology file',
     &          FSREAD3, 'HOURLYT', PROGNAME )
        
C.........  Get MOBILE6 directory from the environment
        MESG = 'Location of MOBILE6 input and output files ' //
     &         '(50 characters or less)'
        CALL ENVSTR( 'SMK_M6PATH', MESG, '.', M6DIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: MOBILE6 files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_M6PATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF        

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

C.........  Ensure that there is at least one activity in the inventory 
C           file, or else this program does not need to be run
        IF( NIACT == 0 ) THEN
            MESG = 'No activities are found in the ' //
     &             'inventory file!  Program cannot be used.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Set up emission process variable names
        CALL EFSETUP( 'NONE', MODELNAM, NEFS, VOLNAM )

C.........  Read emission processes file.  Populate array in MODEMFAC.
        CALL RDEPROC( TDEV )
        
C.........  Read inventory table
        CALL RDCODNAM( VDEV )

C.........  Set up arrays for processing NONHAP values

C.........  Set input and output hydrocarbon names
        INPUTHC = TRIM( VOLNAM )
        OUTPUTHC = 'NONHAP' // TRIM( INPUTHC )

        ALLOCATE( RAWSUBS( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RAWSUBS', PROGNAME )
        
        FNDINPUT = .FALSE.
        FNDOUTPUT = .FALSE.
        K = 0

C.........  Loop through all pollutants        
        DO I = 1, MXIDAT
            IF( INVDNAM( I ) == INPUTHC ) THEN
                FNDINPUT = .TRUE.
                CYCLE
            END IF
            
            IF( INVDNAM( I ) == OUTPUTHC ) THEN
                FNDOUTPUT = .TRUE.
                CYCLE
            END IF

C.............  If requested hydrocarbon is not TOG or VOC, skip rest of loop
            IF( INPUTHC /= 'TOG' .AND. INPUTHC /= 'VOC' ) CYCLE
         
            IF( INVDVTS( I ) /= 'N' ) THEN
            
C.................  Check that pollutant is generated by MOBILE6   
                DO J = 1, NEPOL
                    IF( INVDNAM( I ) == EMTPOL( J ) ) THEN
                        IF( INVDVTS( I ) == 'V' ) THEN
                            K = K + 1
                            RAWSUBS( K ) = INVDNAM( I )
                        ELSE IF( INPUTHC == 'TOG' ) THEN
                            K = K + 1
                            RAWSUBS( K ) = INVDNAM( I )
                        END IF
                        EXIT
                    END IF
                END DO
            END IF
        END DO

C.........  Check that the input hydrocarbon was found
        IF( .NOT. FNDINPUT ) THEN
            MESG = 'Requested MOBILE6 hydrocarbon ' //
     &             TRIM( VOLNAM ) // 'is not in the inventory ' //
     &             'data table.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
        END IF

C.........  If output was not found, set name to blank and set no. polls to zero        
        IF( .NOT. FNDOUTPUT .OR. K == 0 ) THEN
            OUTPUTHC = ' '
            NSUBPOL = 0
        ELSE
            NSUBPOL = K
        END IF

        IF( NSUBPOL > 0 ) THEN
            ALLOCATE( SUBPOLS( NSUBPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SUBPOLS', PROGNAME )
            
            SUBPOLS = RAWSUBS( 1:NSUBPOL )
            
C.............  Write message to log file with pollutants to be subtracted
            MESG = 'NOTE: Emissions from the following pollutants ' //
     &             'will be subtracted from ' // TRIM( INPUTHC ) // 
     &             ' to create ' // TRIM( OUTPUTHC ) // ':'
            CALL M3MESG( MESG )
            
            DO I = 1, NSUBPOL
                MESG = BLANK10 // TRIM( SUBPOLS( I ) ) 
                CALL M3MESG( MESG )
            END DO

        END IF
        
        DEALLOCATE( RAWSUBS )

C.........  Rename emission factors if necessary
        IF( OUTPUTHC /= ' ' ) THEN
            DO I = 1, SIZE( EMTNAM,1 )
                L = INDEX( EMTNAM( I,1 ), ETJOIN )
                L2 = LEN_TRIM( ETJOIN )
                
                IF( EMTNAM( I,1 )( L+L2:IOVLEN3 ) == INPUTHC ) THEN
                    EMTNAM( I,1 )( L+L2:IOVLEN3 ) = OUTPUTHC
                    CYCLE
                END IF
            END DO
        END IF

C.........  Get output directory information from the environment
        MESG = 'Path where emission factors files will be written'
        CALL ENVSTR( 'SMK_EMISPATH', MESG, '.', EMISDIR, IOS )

        IF( IOS /= 0 ) THEN
            MESG = 'WARNING: Emission factors files being placed in ' //
     &             'executable directory because ' // CRLF() //
     &             BLANK10 // 'environment variable SMK_EMISPATH '//
     &             'is not set properly'
            CALL M3MSG2( MESG )
        END IF
        
C.........  Read the GROUP list file into an array
        MESG = 'Reading GROUP list file...'
        CALL M3MSG2( MESG )

C.........  Get the number of lines in the GROUP file     
        NGRPLINES = GETFLINE( GDEV, 'GROUP county list file' )

C.........  If no counties in group file, exit program        
        IF( NGRPLINES == 0 ) THEN
            MESG = 'ERROR: No counties listed in GROUP file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 
        
        ALLOCATE( GRPLIST( NGRPLINES,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPLIST', PROGNAME )

C.........  Read GROUP file
        CALL RDGRPLIST( GDEV, NGRPLINES, GRPLIST )
        
C.........  Read the M6LIST file into an array
        MESG = 'Reading M6LIST file...'
        CALL M3MSG2( MESG )

        CALL RDM6LIST( CDEV )

C.........  Read speed profiles file
        IF( SPDFLAG ) CALL RDSPDPROF( ZDEV )
        
C.........  Allocate memory for the source/scenario number array
        ALLOCATE( SCENLIST( NSRC,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCENLIST', PROGNAME )
        
        SCENLIST = 0

C.........  Start loop over meteorology files
        DO

c            IF( TEMPFLAG ) THEN
C.................  Read header of meteorology file
                IF ( .NOT. DESC3( TNAME ) ) THEN
                
                    MESG = 'Could not get description of file ' // TNAME
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                
C.................  Save header information that will be needed later
                ELSE
                    SDATE  = SDATE3D
                    STIME  = STIME3D
                    TSTEP  = TSTEP3D
                    NSTEPS = MXREC3D
                    NROWS  = NROWS3D
                    
C.....................  Find end date in file description
                    SEARCHSTR = '/END DATE/ '
                    L = LEN_TRIM( SEARCHSTR ) + 1
                    
                    DO I = 1, MXDESC3
                       IF( INDEX( FDESC3D( I ), 
     &                            SEARCHSTR( 1:L ) ) > 0 ) THEN
                           TEMPLINE = FDESC3D( I )
                           IF( CHKINT( TEMPLINE( L+1:L+8 ) ) ) THEN
                               EDATE = STR2INT( TEMPLINE( L+1:L+8 ) )
                               EXIT
                           ELSE
                               EFLAG = .TRUE.
                               EXIT
                           END IF
                       END IF
                       
                       IF( I == MXDESC3 ) THEN
                           EFLAG = .TRUE.
                       END IF
                    END DO

                    IF( EFLAG ) THEN
                        MESG = 'Could not get ending date of file ' //
     &                         TNAME
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF

C.................  Allocate space for meteorology info            
                IF( INITIAL ) THEN
                    ALLOCATE( TEMPCTY( NROWS ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'TEMPCTY', PROGNAME )
                    ALLOCATE( TKHOUR( NROWS, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
                    ALLOCATE( QVHOUR( NROWS, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'QVHOUR', PROGNAME )
                    ALLOCATE( BPHOUR( NROWS, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'BPHOUR', PROGNAME )
                    ALLOCATE( BPDAY( NROWS ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'BPDAY', PROGNAME )
                    ALLOCATE( RHHOUR( NROWS, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'RHHOUR', PROGNAME )
                END IF
                
                TEMPCTY = 0

C.................  Read contents of hourly meteorology file
                MESG = 'Reading hourly meteorology file...'
                CALL M3MSG2( MESG )

C.................  Read county list from file
                IF( .NOT. READ3( TNAME, 'COUNTIES', 1, SDATE, STIME, 
     &                           TEMPCTY ) ) THEN
                    MESG = 'Could not read COUNTIES from ' // TNAME
                    CALL M3EXIT( PROGNAME, SDATE, STIME, MESG, 2 )
                END IF
              
C.................  Read temperature data from file            
                TEMPDATE = SDATE
                TEMPTIME = STIME
                CALL RDHOURTEMP( TNAME, NROWS, TEMPDATE, TEMPTIME )

C.............  If not using temperature files            
c            ELSE
                
C...................  Get start and end dates from enviroment
C                     Have to adjust for 6 am local time start
c            END IF
            
C.............  Create the concatenated MOBILE6 input file
            MESG = 'Writing MOBILE6 input file...'
            CALL M3MSG2( MESG )

C.............  Open new input file
            MDEV = JUNIT()
            WRITE( M6INPUT,94010 ) M6DIR( 1:LEN_TRIM( M6DIR ) ) // 
     &                             '/m6input.' // 
     &                             GRP_NAME( 1:LEN_TRIM( GRP_NAME ) ) //
     &                             '.', SDATE, '.txt'
            OPEN( UNIT=MDEV, FILE=M6INPUT, STATUS='REPLACE', 
     &            ACTION='WRITE', IOSTAT=IOS )
            
            IF( IOS /= 0 ) THEN
                MESG = 'Could not create file ' // M6INPUT
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF     
            
            NUMSCEN = 1
            NUMSRC  = 0
            CALL WRM6INPUT( GRPLIST, NGRPLINES, PDEV, MDEV, 
     &                      TEMPCTY, NROWS, VOLNAM, NUMSCEN, 
     &                      NUMSRC, TEMPFLAG, RHFLAG, SPDFLAG )
            
            CLOSE( MDEV )
            
            NUMSCEN = NUMSCEN - 1
            
C.............  Open file for storing emission factors (check this now rather than
C               waste time running Mobile6)
            FNAME = 'EMISFACS'
            
            IF( FILEOPEN ) THEN
                IF( .NOT. CLOSESET( FNAME ) ) THEN
                    MESG = 'Could not close file set ' // 
     &                      FNAME( 1:LEN_TRIM( FNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    FILEOPEN = .FALSE.
                END IF
            END IF
            
            CALL OPENSEF( NUMSRC, GRP_NAME, SDATE, EDATE, 
     &                    EMISDIR, FNAME )
            FILEOPEN = .TRUE.
            
C.............  Allocate space for storing emission factors
            IF( INITIAL ) THEN
                ALLOCATE( EMISSIONS( MXM6EPR ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMISSIONS', PROGNAME )

                DO I = 1,MXM6EPR
                
C.....................  Calculate maximum values for this emission process
C                       Can't use MAXVAL on Linux
                    DO J = 1, MXM6POLS
                        IF( M6POL2EF( I,J ) > MAXPOL ) THEN
                            MAXPOL = M6POL2EF( I,J )
                        END IF
                    END DO
                    
                    DO J = 1, MXM6VTYP
                        IF( M6VEH2EF( I,J ) > MAXVEH ) THEN
                            MAXVEH = M6VEH2EF( I,J )
                        END IF
                    END DO
                    
                    DO J = 1, MXM6FACS
                        IF( M6FAC2EF( I,J ) > MAXFAC ) THEN
                            MAXFAC = M6FAC2EF( I,J )
                        END IF
                    END DO
                    
                    ALLOCATE( 
     &                  EMISSIONS( I )%PTR( NUMSCEN, MAXPOL, MAXVEH,
     &                                      MAXFAC, 24 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'EMISSIONS%PTR', PROGNAME )
                END DO

                INITIAL = .FALSE.
            END IF
            
            DO I = 1,MXM6EPR
                EMISSIONS( I )%PTR = 0.
            END DO

C.............  Call Mobile6 with M6 input file name
C               Custom driver will use SMOKE log file for any screen output
            MESG = 'Running MOBILE6...' // CRLF()
            CALL M3MSG2( MESG )
            
            CALL SMKDRIVER( M6INPUT, LDEV, EFLAG )
            
            IF( EFLAG ) THEN
                MESG = 'Problem running MOBILE6'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Match up emission factors with sources and write to file
            MESG = 'Writing emission factors to file...' // CRLF()
            CALL M3MSG2( MESG )
            
            CALL WREMFACS( FNAME, NUMSRC, SDATE, VOLNAM )

            IF( TEMPFLAG ) THEN

C.................  Close current temperature file
                IF( .NOT. CLOSE3( TNAME ) ) THEN
                    MESG = 'Could not close file ' // 
     &                     TNAME( 1:LEN_TRIM( TNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
C.................  Open next temperature file if available
             
C.................  Set start date of next file based on end date of current file
                SDATE = EDATE
                CALL NEXTIME( SDATE, STIME, 24*TSTEP )
             
C.................  Construct next file name
                WRITE( TEMPNAME,94010 ) 
     &                 TEMPDIR( 1:LEN_TRIM( TEMPDIR ) ) //
     &                 '/' // GRP_NAME ( 1:LEN_TRIM( GRP_NAME ) ) //
     &                 '.', SDATE, '.ncf'
             
C.................  Check if file exists
                INQUIRE( FILE=TEMPNAME, EXIST=FEXIST )
             
C.................  If file does not exist, we're done            
                IF( .NOT. FEXIST ) EXIT
                
C.................  Set logical file name
                IF( .NOT. SETENVVAR( TNAME, TEMPNAME ) ) THEN
                    MESG = 'Could not set logical file name for ' //
     &                     'file ' // TRIM( TEMPNAME )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                                    
C.................  Open temperature file
                IF( .NOT. OPEN3( TNAME, FSREAD3, PROGNAME ) ) THEN
                    MESG = 'Could not open temperature file ' // 
     &                     TEMPNAME( 1:LEN_TRIM( TEMPNAME ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            ELSE
                EXIT
            END IF

        END DO

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94020   FORMAT( A, I7, A1, I7, A )

        END PROGRAM EMISFAC
