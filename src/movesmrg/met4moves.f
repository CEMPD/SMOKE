
        PROGRAM MET4MOVES

C***********************************************************************
C  program body starts at line 205
C
C  DESCRIPTION:
C       Creates county-based 24-hour temperature profiles and relative  
C       humidity based on gridded meteorology data. Temperatures and 
C       RH can be averaged across counties and different time periods. 
C       Also requires hourly barometric pressure values to compute RH.
C       using barometric pressure and dew point. Tries to 
C       account for missing meteorology data as much as possible.
C
C  PRECONDITIONS REQUIRED:
C       Program Mbsetup has been run
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     2/10: Created by B.H. Baek
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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

C...........   MODULES for public variables
C.........  This module is used for MOBILE6 setup information                
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX,
     &                      DAILY, WEEKLY, MONTHLY, EPISLEN,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NGRID, NCOLS, NROWS, COORD, GDTYP,
     &                     P_ALP, P_BET, P_GAM, XCENT, YCENT, NGRID,
     &                     OFFLAG

        USE MODSURG, ONLY: NSRGREC, IDXSRGA, SFIPSA, SCELLA, SSRGIDA,
     &                     IDXSRGB, NSRGFIPS, NSRGS, SRGFIPS, SRGLIST,
     &                     FIPCELL, SRGFRAC, SRGCSUM, NCELLS, SFRACA,
     &                     SRGFMT, NTSRGDSC, SRGFNAM, SRGFCOD, SRGFREG

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: LSTCYPOP, STCYPOPYR, USEDAYLT,
     &                     NCOUNTRY, NSTATE,   NCOUNTY,
     &                     CTRYCOD,  STATCOD,  CNTYCOD,
     &                     CTRYNAM,  STATNAM,  CNTYNAM,
     &                     CTRYPOPL, STATPOPL, CNTYPOPL,
     &                     CNTYTZON, CNTYTZNM

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: MINTEMP, MAXTEMP, TASRC, QVSRC, PRESSRC,
     &                    TKHOUR, RHHOUR, NDAYSRC, MAXTSRC, MINTSRC,
     &                    MAXTFUEL, MINTFUEL, RHFUEL, TKFUEL, FUELIDX,
     &                    NFUEL
     
        IMPLICIT NONE
        
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
c        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables
c        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

        CHARACTER(2)    CRLF
        LOGICAL         DSCM3GRD
        INTEGER         GETIFDSC
        INTEGER         GETFLINE
        INTEGER         GETEFILE
        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         INDEX1
        INTEGER         INTLIST
        INTEGER         FIND1
        INTEGER         FIND1FIRST
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         SECSDIFF
        LOGICAL         SETENVVAR
        INTEGER         WKDAY
        
        EXTERNAL     CRLF, DSCM3GRD, GETIFDSC, GETFLINE, ENVINT, FIND1
     &               ENVREAL, INDEX1, MMDDYY, PROMPTFFILE, PROMPTMFILE, 
     &               SECSDIFF, SETENVVAR, WKDAY, GETEFILE, INTLIST,
     &               FIND1FIRST
        
C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$' ! CVS release tag
        INTEGER, PARAMETER :: MXVAR = 20

C...........   LOCAL VARIABLES and their descriptions:
        
C...........   Gridded meteorology data (dim: NGRID)
        REAL   , ALLOCATABLE :: TA( : )   !  one layer of temperature
        REAL   , ALLOCATABLE :: QV( : )   !  water vapor mixing ratio
        REAL   , ALLOCATABLE :: PRES( : ) !  pressure
        REAL   , ALLOCATABLE :: TKREFHR( : ) !  ref county temp profile
 
C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C...........  Allocatable per-source arrays
        INTEGER, ALLOCATABLE :: DAYBEGT  ( : ) ! daily start time HHMMSS
        INTEGER, ALLOCATABLE :: DAYENDT  ( : ) ! daily end time HHMMSS
        LOGICAL, ALLOCATABLE :: LDAYSAV  ( : ) ! true: src uses daylight time
        INTEGER, ALLOCATABLE :: CNTYSRC  ( : ) ! FIPS code and averaging values
        INTEGER, ALLOCATABLE :: TZONES   ( : ) ! county-specific time zones
        INTEGER                 SRGIDS  ( MXVAR ) ! listing of surrogates

C...........  Allocatable arrays for met data
        INTEGER, ALLOCATABLE :: METDAYS( : ) ! dimension: nsteps in episode, 
                                             ! value indicates which met data file covers that hour
        CHARACTER(256), ALLOCATABLE :: METLIST( : ) ! listing of met file names

C...........   File units and logical names:
        INTEGER      ADEV  ! unit number for county fuelmonth file
        INTEGER      CDEV  ! unit number for co/st/cy file
        INTEGER      DDEV  ! unit number for daily group file
        INTEGER      EDEV  ! unit number for episode group file
        INTEGER      GDEV  ! tmp unit number for individual surrogate file
        INTEGER      IDEV  ! tmp unit number if ENAME is map file
        INTEGER      LDEV  ! unit number for log file
        INTEGER      MDEV  ! unit number for monthly group file
        INTEGER      PDEV  ! unit number for speeds summary file (SPDSUM)
        INTEGER      QDEV  ! unit number for surrogate(s) file
        INTEGER      SDEV  ! unit number for surrogate(s) file
        INTEGER      TDEV  ! unit number for meteorology list file
        INTEGER      WDEV  ! unit number for weekly group file
        INTEGER      ODEV  ! unit number for tmp surrogate file
        INTEGER      ODEV1 ! unit number for SMOKE-ready output file
        INTEGER      ODEV2 ! unit number for MOVES-ready output file
        INTEGER      XDEV  ! unit number for mobile x-ref file

        CHARACTER(16) DNAME   ! logical name for daily output ungridded hourly temps
        CHARACTER(16) INAME   ! tmp name for inven file of unknown fmt
        CHARACTER(16) METNAME ! logical name for meteorology files
        CHARACTER(16) MNAME   ! logical name for monthly output hourly temps
        CHARACTER(16) PNAME   ! logical name for episode output hourly temps
        CHARACTER(16) UNAME   ! logical name for ungridding-matrix input file
        CHARACTER(16) WNAME   ! logical name for weekly output hourly temps
                
C...........   Other local variables:
        INTEGER    I, IC, J, K, L, N, NRH, NR, S, T, TT, T2, V  ! Counters and pointers

        INTEGER    EPI_SDATE      ! episode start date from E.V. (YYYYDDD)
        INTEGER    EPI_STIME      ! episode start time from E.V. (HHMMSS)
        INTEGER    EPI_RUNLEN     ! episode duration from E.V. (HHMMSS)
        INTEGER    EPI_NSTEPS     ! episode number of time steps
        INTEGER    EPI_EDATE      ! episode ending date based on ERUNLEN
        INTEGER    EPI_ETIME      ! episode ending time based on ERUNLEN
        
        INTEGER    ARRAYPOS    ! position in 24-hour arrays
        INTEGER    CURCNTY     ! tmp current processing county
        INTEGER    DAY         ! tmp day of week number
        INTEGER    DDATE       ! output date for daily counties or date difference
        INTEGER    DMONTH      ! month difference 
        INTEGER    DUMMYTIME   ! dummy time variable to use in calls to NEXTIME
        INTEGER    EDATE       ! ending input date counter (YYYYDDD) in GMT
        INTEGER    ETIME       ! ending input time counter (HHMMSS)  in GMT
        INTEGER    FIP         ! tmp inventory county
        INTEGER    FILENUM     ! file number of current meteorology file
        INTEGER    IOS         ! temporary I/O status
        INTEGER    JDATE       ! input date counter (YYYYDDD) in GMT
        INTEGER    JTIME       ! input time counter (HHMMSS)  in GMT
        INTEGER    LDATE       ! date from previous loop iteration
        INTEGER    MDATE       ! output date for monthly counties
        INTEGER    SMONTH      ! month of start date
        INTEGER    SDAY        ! date of month (start) 
        INTEGER    EMONTH      ! month of end date
        INTEGER    EDAY        ! date of month (end) 
        INTEGER    TMONTH      ! tmp month of end date
        INTEGER    TDAY        ! tmp date of month (end) 
        INTEGER    MONTH       ! tmp month
        INTEGER    CURMONTH    ! current month
        INTEGER    TMPMONTH    ! previous month
        INTEGER    FUELMONTH   ! current fuelmonth
        INTEGER    PRVFMONTH   ! previous fuelmonth
        INTEGER    METNGRID    ! no. grid cells in met data
        INTEGER    NLINES      ! no. lines in met list file
        INTEGER    NSRC        ! no. source (=counties)
        INTEGER    NVARS       ! no. surrogates
        INTEGER    NF          ! current and prvious county fuelmonths
        INTEGER    PRNF        ! prvious county fuelmonths
        INTEGER    NMON        ! no of fuelmonth per refcounty
        INTEGER    NFMON       ! tmp no of fuelmonths per refcounty
        INTEGER    NSTEPS      ! number of time steps to process temperature data
        INTEGER    OTIME       ! output time in local time
        INTEGER    POS         ! position in time step loop
        INTEGER    PDTEMP      ! temp increment for rateperdistance lookup table
        INTEGER    PVTEMP      ! temp increment for ratepervehicle lookup table
        INTEGER    PPTEMP      ! temp increment for rateperprofile lookup table
        INTEGER    RH_STRHR    ! start hour to read RH from met file
        INTEGER    RH_ENDHR    ! ending hour to read RH from met file
        INTEGER    REFCOUNTY   ! ref. county FIPS code
        INTEGER    INVCOUNTY   ! inv. county FIPS code
        INTEGER    PRCOUNTY    ! previous ref. county
        INTEGER    RDATE       ! date to read met file
        INTEGER    SDATE       ! output start date
        INTEGER    STIME       ! output start time
        INTEGER    TDATE       ! temporary julian date
        INTEGER    TTIME       ! temporary time
        INTEGER    SRGNCOLS    ! surrogate no of cols
        INTEGER    SRGNROWS    ! surrogate no of rows
        INTEGER    TMPMNTH     ! temporary month
        INTEGER    TSPREAD     ! time spread: difference between TZMAX and TZMIN
        INTEGER    TZONE       ! zone to determine output days
        INTEGER    TZMIN       ! minimum time zone in inventory
        INTEGER    TZMAX       ! maximum time zone in inventory
        INTEGER    WDATE       ! output date for weekly counties

        REAL       RHSUM               ! sum of RH
        REAL       RHAVG               ! avg of RH
        REAL       RHREFSUM            ! sum of RH for ref county
        REAL       RHREFAVG            ! avg of RH for ref county
        REAL       MAXTREF, MINTREF    ! ref county min/max temperatures

        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: COMPLETE = .FALSE.  !  true: program successful complete
        LOGICAL :: GRID_ERR = .FALSE.  !  true: error found in grid settings
        LOGICAL :: DAYAVER  = .FALSE.  !  true: daily averaging
        LOGICAL :: MONAVER  = .FALSE.  !  true: monthly averaging
        LOGICAL :: EPIAVER  = .FALSE.  !  true: episode averaging
        LOGICAL :: FUELAVER = .FALSE.  !  true: fuelmonth averaging
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs
        LOGICAL :: MONOPEN  = .FALSE.  !  true: monthly fuelmonth is processed
        LOGICAL :: FILEOPEN = .FALSE.  !  true: met file is open
        LOGICAL :: FND_DATA = .FALSE.  !  true: found met data for this hour
        LOGICAL :: ALT_DATA = .FALSE.  !  true: using alternate data for this hour

        CHARACTER(16)      COORUNIT    !  coordinate system projection units
        CHARACTER(80)      GDESC       !  grid description
        CHARACTER(256)     CFNAME      !  fuelmonth input file
        CHARACTER(16)      SRG_CNTRY   !  surrogate country
        CHARACTER(256)     LINE        !  line buffer
        CHARACTER(256)     FULLNAME    !  full file name
        CHARACTER(512)     METFILE     !  tmp physical file name
        CHARACTER(512)     PREVFILE    !  previous physical file name
        CHARACTER(IOVLEN3) AVG_TYPE    !  averaging method name
        CHARACTER(IOVLEN3) TVARNAME    !  temperature variable name
        CHARACTER(IOVLEN3) PRESNAME    !  pressure variable name
        CHARACTER(IOVLEN3) MIXNAME     !  mixing ratio name
        CHARACTER(200)     TEMPDIR     !  directory for output files
        CHARACTER(300)     MESG        !  message buffer

        CHARACTER(16) :: PROGNAME = 'MET4MOVES'  !  program name

C***********************************************************************
C   begin body of program MET4MOVES

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get file name for country, state, and county file, with time zones
        CDEV = PROMPTFFILE(
     &         'Enter logical name for COUNTRY, STATE, AND ' //
     &         'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  Open input surrogates file
        QDEV = PROMPTFFILE(
     &         'Enter logical name for Surrogate Description file',
     &         .TRUE., .TRUE., 'SRGDESC', PROGNAME )

C.........  Open met list file
        TDEV= PROMPTFFILE(
     &         'Enter logical name for METEOROLOGY LIST file'
     &         ,.TRUE., .TRUE., 'METLIST', PROGNAME )

C.........  Open mobile county x-ref file to determine representative counties
        XDEV = PROMPTFFILE(
     &           'Enter logical name for MCXREF cross-reference file',
     &           .TRUE., .TRUE., 'MCXREF', PROGNAME )

C.........  Obtain episode settings from the environment...
C.........  Define averaging method.
        CALL ENVSTR( 'AVERAGING_METHOD', MESG,' ',AVG_TYPE, IOS )

        MESG ='Define averaging method for meteorology processing'
        CALL M3MSG2( MESG )

        CALL UPCASE( AVG_TYPE )

        IF( IOS .NE. 0 ) THEN
            MESG = 'ERROR: AVERAGE_METHOD environment variable ' //
     &             ' is not defined for meteorology processing'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check group file unit numbers to see which types of averaging are needed
        IF( AVG_TYPE == 'DAILY'   ) THEN
            DAYAVER  = .TRUE.
        ELSE IF( AVG_TYPE == 'MONTHLY' ) THEN
            MONAVER  = .TRUE.
        ELSE IF( AVG_TYPE == 'EPISODE' ) THEN
            EPIAVER  = .TRUE.
        ELSE
            MESG = 'ERROR: AVERAGE_METHOD environment variable ' //
     &             TRIM( AVG_TYPE ) // 'is not recognized.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Get temperature increments for rateperdistance lookup table
        MESG = 'Temperature increment for rateperdistance'
        PDTEMP = ENVINT( 'PD_TEMP_INCREMENT', MESG, 5, IOS )

C.........  Get temperature increments for ratepervehicle lookup table
        MESG = 'Temperature increment for ratepervehicle'
        PVTEMP = ENVINT( 'PV_TEMP_INCREMENT', MESG, 5, IOS )

C.........  Get temperature increments for ratepervehicle lookup table
        MESG = 'Temperature increment for rateperprofile'
        PPTEMP = ENVINT( 'PP_TEMP_INCREMENT', MESG, 10, IOS )

C.........  Get episode starting date and time and ending date
        MESG = 'Episode start date (YYYYDDD)'
        EPI_SDATE = ENVINT( 'EPI_STDATE', MESG, 0, IOS )
        
        MESG = 'Episode start time (HHMMSS)'
        EPI_STIME = ENVINT( 'EPI_STTIME', MESG, 0, IOS )
         
        MESG = 'Episode end date (YYYYDDD)'
        EPI_EDATE = ENVINT( 'EPI_ENDATE', MESG, 0, IOS )

        MESG = 'Episode end time (HHMMSS)'
        EPI_ETIME = ENVINT( 'EPI_ENDTIME', MESG, 230000, IOS )

C.........  Find the total number of time steps
        EPI_NSTEPS = 1 + 
     &       SECSDIFF( EPI_SDATE,EPI_STIME,EPI_EDATE,EPI_ETIME ) / 3600

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Get RH averaging start and ending hours
        MESG = 'RH start averaging hour (HHMMSS)'
        RH_STRHR = ENVINT( 'RH_STR_HOUR', MESG, 60000, IOS )
        RH_STRHR = RH_STRHR / 10000
        
        MESG = 'RH ending average hour (HHMMSS)'
        RH_ENDHR = ENVINT( 'RH_END_HOUR', MESG, 180000, IOS )
        RH_ENDHR = RH_ENDHR / 10000
                
C.........  Get the name of the temperature variable
        MESG = 'Temperature variable name'
        CALL ENVSTR( 'TVARNAME', MESG, 'TEMP2', TVARNAME, IOS )

C.........  Set default names for additional variables
        PRESNAME = 'PRES'
        MIXNAME = 'QV'

C.........  Read the surrogates header and initialize the grid description
C.........  Open/read SRGDESC file
        CALL RDSRGDESC( QDEV )

C.........  Open output file for temporary combined surrogate file
        FULLNAME = 'TMP_COMBINED_SRG.txt'
        IF( .NOT. SETENVVAR( 'TMP_SRG_FILE', FULLNAME )) THEN
             MESG = 'Could not set logical file name of file ' 
     &              // TRIM( FULLNAME )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ODEV = GETEFILE( 'TMP_SRG_FILE', .FALSE., .TRUE., PROGNAME )

C.........  Get name of surrogate IDs to use
        MESG = 'Country of surrogate files'
        CALL ENVSTR( 'SRG_COUNTRY', MESG, 'USA', SRG_CNTRY, IOS )
        
        MESG = 'Spatial surrogate IDs'
        IF( .NOT. INTLIST( 'SRG_LIST', MESG, 20, NVARS, SRGIDS ) ) THEN
            MESG = 'Could not read list of surrogate IDs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Prompt for and open I/O API output file(s)...
        DO I = 1, NVARS
 
            L = 0
            DO K = 1,NTSRGDSC
                CALL UPCASE( SRGFREG( K ) )
                IF( SRG_CNTRY    == SRGFREG( K ) ) THEN
                    IF( SRGIDS( I )  == SRGFCOD( K ) ) L = K
                END IF
            END DO
 
            IF( L < 1 ) CYCLE
            
            CALL GETENV( 'SRGPRO_PATH', TEMPDIR )
            WRITE( FULLNAME, '(3A)' ) TRIM( TEMPDIR ),'/',TRIM( SRGFNAM( L ) )

C.........  Set logical file name
            IF( .NOT. SETENVVAR( 'TMP_SRG', FULLNAME )) THEN
                MESG = 'Could not set logical file ' //
     &                 'name of file ' // TRIM( FULLNAME )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Get the number of lines in the surrogate description file desription file
            GDEV = GETEFILE( 'TMP_SRG',.TRUE., .TRUE., PROGNAME )
            NLINES = GETFLINE( GDEV, 'Surrogate file' )
            DO J = 1, NLINES
                READ(  GDEV,'(A)', END=999, IOSTAT=IOS ) LINE
                WRITE( ODEV,'(A)' ) TRIM( LINE )
999         END DO

            CLOSE( GDEV )

        END DO
        
        CLOSE( ODEV )

        SDEV = GETEFILE( 'TMP_SRG_FILE',.TRUE., .TRUE., PROGNAME )

C.........  Also get the format of the surrogates file.
        CALL RDSRGHDR( .FALSE., SDEV, SRGFMT )
        SRGNCOLS = NCOLS
        SRGNROWS = NROWS

        L = LEN_TRIM( SRGFMT )
        MESG = 'NOTE: Input surrogates are ' // SRGFMT( 1:L ) //
     &         ' format.'
        CALL M3MSG2( MESG )

C.........  Obtain settings from the environment...
C.........  Get grid name from the environment and read grid parameters
        IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D,
     &       COORUNIT, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,
     &       YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &       YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )
        IF ( EFLAG ) THEN
            MESG = 'Problem with gridded input data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for and read the gridding surrogates file,
C           extracting data for a subgrid, if necessary
        CALL M3MSG2( 'Processing gridding surrogate(s)...' )

        CALL RDSRG( .FALSE., SDEV, SRGFMT, SRGNROWS, SRGNCOLS )

        MESG = 'NOTE: A list of surrogates used in the modeling'
        CALL M3MESG( MESG )

        DO I = 1,NSRGS
            WRITE( MESG,94010 ) 'Surrogate ID ::', SRGLIST( I )
            CALL M3MESG( MESG )
        END DO

C.........  Read region codes file
        CALL RDSTCY( CDEV, NSRGFIPS, SRGFIPS )

C.........  Read county x-ref file to determine representative counties
        CALL RDMXREF( XDEV, NSRGFIPS, SRGFIPS  )

C.........  open county-specific fuelmonth input files in file
        MESG = 'Enter logical name for a reference county-' //
     &         'specific FuelMonth input file'
        ADEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'MFMREF', PROGNAME )

        CALL RDFMREF( ADEV ) 

C.........  define a number of target sources (=inventory counties)
        NSRC = NINVC

C.........  Allocate arrays for county time zone
        ALLOCATE( TZONES( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONES', PROGNAME )
        TZONES = 0

C.........  Assign time zone to inventory counties
        DO I = 1, NSRC
            FIP = MCREFSORT( I,1 )
            J = FIND1( FIP, NCOUNTY, CNTYCOD )
            IF( J < 1 ) THEN
                WRITE( MESG,94010 ) 'ERROR: Could not find time zone '//
     &               'for county :', FIP, ' from COSTCY file' 
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                TZONES( I ) = CNTYTZON( J )
            END IF
        END DO

C.........  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

C.........  Adjust TZMIN for possibility of daylight savings
        TZMIN = MAX( TZMIN - 1, 0 )

C.........  Calculate time spread based on min and max time zone
        TSPREAD = TZMAX - TZMIN

C.........  Calculate required starting date and time based on episode settings
        SDATE = EPI_SDATE

C.........  For MOVES model, earliest time required will be 12 a.m. in time zone closest to GMT
        STIME = TZMIN - TZONE          ! starting time in output time zone
        
C.........  Make sure the starting time is between 0 and 23
        IF( STIME < 0 ) THEN
            STIME = STIME + 24
        ELSE
            STIME = MOD( STIME, 24 )
        END IF
        STIME = STIME*10000

C.........  If the episode start time is earlier than our calculated start time,
C           we need to set the starting date back one day
        IF( EPI_STIME < STIME ) THEN
            CALL NEXTIME( SDATE, STIME, -24*10000 )
        END IF
        
C.........  Calculate required ending date and time based on episode settings
        EDATE = EPI_EDATE
        ETIME = EPI_ETIME
        
C.........  For MOVES model, latest time required will be 11 p.m. in time zone farthest from GMT
        ETIME = TZMAX - TZONE - 1        ! ending time in output time zone
        
C.........  Make sure the ending time is between 0 and 23
        IF( ETIME < 0 ) THEN
            ETIME = ETIME + 24
        ELSE
            ETIME = MOD( ETIME, 24 )
        END IF
        ETIME = ETIME*10000

C.........  If the episode ending time is later than calculated end time,
C           set the ending date forward one day
        IF( EPI_ETIME > ETIME ) THEN
            CALL NEXTIME( EDATE, ETIME, 24*10000 )
        END IF

C.........  Convert start and end dates and times back to GMT
        CALL NEXTIME( SDATE, STIME, TZONE*10000 )
        CALL NEXTIME( EDATE, ETIME, TZONE*10000 )

C.........  Find the total number of time steps
        NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600

C.........  Configure the month(s) of modeling period
        CALL DAYMON( SDATE, SMONTH, SDAY )
        CALL DAYMON( EDATE, EMONTH, EDAY )
        DDATE  = EDATE  - SDATE
        DMONTH = EMONTH - SMONTH

        IF( MONAVER ) THEN

            IF( DDATE < 15 ) THEN
                MESG = 'ERROR: MONTHLY averaging method '//
     &              'is not applicable (Episode period must > 15 days).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( DMONTH == 0 ) THEN
                MESG = 'ERROR: MONTHLY averaging method '//
     &              'is not applicable (Must cross the month).'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( DMONTH == 1 .AND. EDAY > 15 ) THEN
                WRITE( MESG,94010 ) 'ERROR: MONTHLY averaging method '//
     &              'is not applicable (MUST cross the month).'//CRLF() 
     &               //BLANK5//'End date', EDAY,'on month of',  
     &               EMONTH,' must cross the month'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF
        END IF

        IF( EPIAVER ) THEN
            IF( DDATE < 2 ) THEN
                MESG = 'ERROR: EPISODE averaging method requires '//
     &            'more than 2 days for the episode period'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Get number of lines in met list file
        NLINES = GETFLINE( TDEV, 'METLIST file' )

C.........  Allocate memory for storing the met file information
        ALLOCATE( METLIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METLIST', PROGNAME )
        ALLOCATE( METDAYS( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METDAYS', PROGNAME )
        
        METLIST = ' '    ! array
        METDAYS = 0
  
C.........  Store lines of METLIST file
        CALL RDLINES( TDEV, 'METLIST file', NLINES, METLIST )

        MESG = 'Checking meteorology files...'
        CALL M3MSG2( MESG )

        METNAME = 'METFILE'

C.........  Loop through all meteorology files
        DO N = 1, NLINES

C.............  Close previous file if needed
            IF( FILEOPEN ) THEN
                IF( .NOT. CLOSE3( METNAME ) ) THEN
                    MESG = 'Could not close meteorology file ' //
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    FILEOPEN = .FALSE.
                END IF
            END IF

C.............  Get physical file name for current iteration        
            METFILE = METLIST( N )  

C.............  Set logical file name
            IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not set logical file name ' //
     &                 'for file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Try to open file   
            IF( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not open meteorology file ' //
     &                 TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            ELSE
                FILEOPEN = .TRUE.
            END IF

C.............  Try to get description from file            
            IF( .NOT. DESC3( METNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not get description of ' //
     &                 'meteorology file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Check that requested variables are in the file
            J = INDEX1( TVARNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( TVARNAME ) // 
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF
            
            J = INDEX1( PRESNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( PRESNAME ) //
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF
            
            J = INDEX1( MIXNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TRIM( MIXNAME ) //
     &                 ' in file ' // TRIM( METFILE )
                CALL M3MESG( MESG )
            END IF

            IF( EFLAG ) CYCLE

C.............  Initialize (or check) reference grid with meteorology data
            CALL CHKGRID( METNAME, 'GRID', 0, GRID_ERR )
            METNGRID = NGRID

C............. If the dimensions were in error, print message and cycle
            IF( GRID_ERR ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Grid in meteorology file ' //
     &                 TRIM( METFILE ) // ' is inconsistent ' //
     &                 'with previous files.' 
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Build METDAYS array linking met files to a specific date and time

C............. Find starting position in METDAYS array
            POS = 1 + SECSDIFF( SDATE, STIME, SDATE3D, STIME3D ) / 3600

C.............  Make sure the met data is within the episode            
            IF( POS + MXREC3D - 1 <= 0 .OR. POS > NSTEPS ) CYCLE

C.............  Fill in array for number of steps in current met file            
            DO L = POS, POS + MXREC3D - 1
                IF( L <= 0 ) CYCLE
                IF( L > NSTEPS ) EXIT
                METDAYS( L ) = N
            END DO

        END DO

C.........  Exit if there was a problem with the meteorology files
        IF( EFLAG ) THEN
            MESG = 'Problem checking meteorology files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check for missing meteorology data
        MESG = 'Checking for missing meteorology data...'
        CALL M3MSG2( MESG )

C.........  Check that all days are covered
        IF( DAYAVER ) THEN
            T = 1
            
C.............  Loop over all days in episode
            DO

C.................  Make sure we're within episode bounds
                IF( T > NSTEPS ) EXIT

C.................  Check all hours in current day, accounting for time zone spread            
                DO I = 0, 23 + TSPREAD
                    K = T + I

C.....................  Double check episode bounds
                    IF( K > NSTEPS ) EXIT

C.....................  If no met data for current step, try to find data
                    IF( METDAYS( K ) == 0 ) THEN

C.........................  Try 24 hours back
                        IF( K - 24 > 0 ) THEN
                            IF( METDAYS( K - 24 ) > 0 ) THEN
                                METDAYS( K ) = - METDAYS( K - 24 )
                                CYCLE
                            END IF

C.........................  Try 12 hours forward
                        ELSE IF( K + 24 < NSTEPS ) THEN
                            IF( METDAYS( K + 24 ) > 0 ) THEN
                                METDAYS( K ) = - METDAYS( K + 24 )
                                CYCLE
                            END IF
                        END IF

C.........................  No data available, exit with error                    
                        MESG = 'Meteorology data does not cover ' //
     &                         'requested episode.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                END DO   ! end loop over hours in day

C.................  Skip to start of next day
                T = T + 24
            END DO   ! end loop over days in episode

        END IF   ! end day averaging check

C.........  Check that all months are covered
        IF( MONAVER ) THEN
            T = 1
            JDATE = SDATE
            JTIME = STIME
            CALL DAYMON( JDATE, MONTH, DAY )

C.............  Loop over all months in episode     
            DO

C.................  Check episode bounds
                IF( T > NSTEPS ) EXIT

C.................  Check hours in day            
                DO I = 0, 23 + TSPREAD
                    K = T + I
                    
C.....................  Check episode bounds
                    IF( K > NSTEPS ) EXIT
                    
C.....................  Check for met data at this step                
                    IF( METDAYS( K ) <= 0 ) THEN
                        FND_DATA = .FALSE.

C.........................  Loop through previous days in month
                        DO J = 1, 15
                            TDATE = JDATE
                            TTIME = JTIME

C.............................  Check episode bounds                        
                            IF( K - J*24 > 0 ) THEN

C.................................  Make sure it's still the same month
                                CALL NEXTIME( TDATE, TTIME, -J*240000 )
                                CALL DAYMON( TDATE, TMPMNTH, DAY )

                                IF( TMPMNTH == MONTH ) THEN
                                    IF( METDAYS( K - J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF
                                    
C.................................  Otherwise, it's the previous month, so exit
                                ELSE
                                    EXIT
                                END IF
                                
C.............................  Otherwise, it's too far back, exit
                            ELSE
                                EXIT
                            END IF
                        END DO   ! end loop over previous days
 
C.........................  Skip rest of loop if we've found data                        
                        IF( FND_DATA ) CYCLE

C.........................  Loop through remaining days in month                    
                        DO J = 1, 15
                            TDATE = JDATE
                            TTIME = JTIME
                        
C.............................  Check episode bounds
                            IF( K + J*24 < NSTEPS ) THEN

C.................................  Make sure it's still the same month
                                CALL NEXTIME( TDATE, TTIME, J*240000 )
                                CALL DAYMON( TDATE, TMPMNTH, DAY )
                                IF( TMPMNTH == MONTH ) THEN
                                    IF( METDAYS( K + J*24 ) > 0 ) THEN
                                        FND_DATA = .TRUE.
                                        EXIT
                                    END IF
                                    
C.................................  Otherwise, it's the next month, so exit
                                ELSE
                                    EXIT
                                END IF    
                                
C.............................  Otherwise, it's too far forward, exit
                            ELSE
                                EXIT        
                            END IF
                        END DO   ! end loop over remaining days
                    
                        IF( FND_DATA ) CYCLE

C.....................  Still no data, exit with error                         
                        MESG = 'Meteorology data does not cover ' //
     &                         'requested episode.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                
                END DO   ! end loop over hours in day
            
C.................  Advance T to next first day of next month at midnight            
                T = T + 1
                CALL NEXTIME( JDATE, JTIME, 10000 )
                CALL DAYMON( JDATE, MONTH, DAY )
                IF( DAY == 1 .AND. JTIME == 0 ) EXIT

            END DO   ! end loop over months in episode

        END IF   ! end month averaging check
        
C.........  Check that the episode is covered
        IF( EPIAVER ) THEN

C.............  Loop over hours in day accounting for time zones
            DO T = 1, 24 + TSPREAD
                K = T

C.................  Check episode bounds
                IF( K > NSTEPS ) EXIT
                
C.................  Check for met data at this step
                IF( METDAYS( K ) <= 0 ) THEN
                    FND_DATA = .FALSE.

C.....................  Loop through additional days in episode
                    DO
                        K = K + 24

C.........................  Double check episode bounds
                        IF( K > NSTEPS ) EXIT

C.........................  If found data, go on to next hour                        
                        IF( METDAYS( K ) > 0 ) THEN
                            FND_DATA = .TRUE.
                            EXIT
                        END IF
                    END DO   ! end loop over days in episode
                
                    IF( FND_DATA ) CYCLE

C.....................  Still no data, exit with error                
                    MESG = 'Meteorology data does not cover ' //
     &                     'requested episode.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

            END DO   ! end loop over hours in day

        END IF   ! end episode averaging check

C.........  Get the info of counties
        ALLOCATE( CNTYSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CNTYSRC', PROGNAME )
        CNTYSRC = 0
        
C.........  Assign inventory counties to new array
        CNTYSRC( : ) = MCREFSORT( :, 1 )


C.........  Open and write out the header 
C           information for SMOKE and MOVES-ready output files

        ODEV1 = PROMPTFFILE(
     &       'Enter logical name for SMOKE-ready output files',
     &       .FALSE., .TRUE., 'SMOKE_OUTFILE', PROGNAME )
        
C.........  Define temporal resolution header
        WRITE( ODEV1,'(A)' )'#DESC SMOKE-ready input file'
        WRITE( ODEV1,94010 )'#AVERAGING_METHOD: ' // TRIM( AVG_TYPE ) //
     &      '  (', EPI_SDATE, '-', EPI_EDATE, ')'
        WRITE( ODEV1,'(A)' )'#DATA FIPS,fuelmonthID,monthID,'//
     &       'JulianDate,avgRH,minimum_Temp,maximum_Temp'
        
C.........  Open output file
        ODEV2 = PROMPTFFILE(
     &       'Enter logical name for MOVES-ready output files',
     &       .FALSE., .TRUE., 'MOVES_OUTFILE', PROGNAME )

        WRITE( ODEV2,'(A)' )'#DESC MOVES-ready input file'
        WRITE( ODEV2,94010 )'#AVERAGING_METHOD: ' // TRIM( AVG_TYPE ) //
     &      '  (', EPI_SDATE, '-', EPI_EDATE, ')'
        WRITE( ODEV2,'(A)' )'#DATA FIPS,MonthID,Temperature'//
     &          'ProfileID,RH,temp1(min),temp2(max),,,,,,,,,,,,temp24'
        WRITE( ODEV2,'(A,I5)' ) 'PD_TEMP_INCREMENT ' , PDTEMP
        WRITE( ODEV2,'(A,I5)' ) 'PV_TEMP_INCREMENT ' , PVTEMP
        WRITE( ODEV2,'(A,I5)' ) 'PP_TEMP_INCREMENT ' , PPTEMP

C.......................................................................
C.........  Allocate met variable arrays
        ALLOCATE( TA( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( QV( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QV', PROGNAME )
        ALLOCATE( PRES( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )

C.........  dates/daylight saving arrays        
        ALLOCATE( DAYBEGT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( LDAYSAV( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        
C.........  Create array of which sources are affected by daylight savings
        CALL GETDYSAV( NSRC, CNTYSRC, LDAYSAV )

C.........  Source met variable arrays
        ALLOCATE( MAXTSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXTSRC', PROGNAME )
        ALLOCATE( MINTSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINTSRC', PROGNAME )
        ALLOCATE( TASRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
        ALLOCATE( QVSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QVSRC', PROGNAME )
        ALLOCATE( PRESSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRESSRC', PROGNAME )
        MAXTSRC = BADVAL3
        MINTSRC = -1*BADVAL3

C.........  Allocate memory for storing meteorology profiles
        ALLOCATE( TKHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKHOUR', PROGNAME )
        ALLOCATE( RHHOUR( NSRC, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RHHOUR', PROGNAME )
        ALLOCATE( NDAYSRC( NSRC,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NDAYSRC', PROGNAME )
        TKHOUR = 0.0
        RHHOUR = 0.0
        NDAYSRC = 0

C.........  Process meteorology data...
        MESG = 'Processing meteorology data using variables ' //
     &         TRIM( TVARNAME ) // ', ' // TRIM( MIXNAME ) // 
     &         ', ' // TRIM( PRESNAME ) // '...'
        CALL M3MSG2( MESG )

C.........  Loop through days/hours of meteorology files
        DDATE = SDATE
        OTIME = 0
        WDATE = SDATE
        MDATE = SDATE
        
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        
        NF = 0

C.........  Process max/min temperatures and avg RH for county-specific
C           fuelmonth

C...........  Define max no of fuelmonth
        NFUEL = ( EMONTH - SMONTH ) + 1

C...........  Allocate fuelmonth arrays
        ALLOCATE( FUELIDX( NREFC,NFUEL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FUELIDX', PROGNAME )
        ALLOCATE( RHFUEL( NREFC,NFUEL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RHFUEL', PROGNAME )
        ALLOCATE( TKFUEL( NREFC,NFUEL,24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKFUEL', PROGNAME )
        ALLOCATE( TKREFHR( 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKREFHR', PROGNAME )
        ALLOCATE( MAXTFUEL( NREFC,NFUEL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAXTFUEL', PROGNAME )
        ALLOCATE( MINTFUEL( NREFC,NFUEL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINTFUEL', PROGNAME )
        RHFUEL   = 0.0
        TKFUEL   = 0.0
        TKREFHR  = 0.0
        MAXTFUEL = BADVAL3
        MINTFUEL = -1*BADVAL3

        N = 0
        DO I = 1,NFUEL
            FUELIDX( :,I ) = SMONTH + N
            N = N + 1
        END DO

C...........  loop over hours
        DO T = 1, NSTEPS

C.............  When new day...
            IF ( JDATE /= LDATE ) THEN
C.................  Set start and end hours of day for all sources
                CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, .TRUE.,
     &                         DAYBEGT, DAYENDT )
            END IF

C.............  Determine input file for this hour
            POS = T
            
C.............  Get file number for current iteration
            FILENUM = METDAYS( POS )
            IF( FILENUM <= 0 ) THEN
                ALT_DATA = .TRUE.
            ELSE
                ALT_DATA = .FALSE.
            END IF

C.............  Skip file opening when not doing day averaging and using alternate data
            IF( .NOT. ALT_DATA ) THEN
            
C.................  Get file name
                METFILE = METLIST( ABS( FILENUM ) )

C.................  Close previous file if needed
                IF( METFILE .NE. PREVFILE ) THEN
                    IF( FILEOPEN ) THEN
                        IF( .NOT. CLOSE3( METNAME ) ) THEN
                            MESG = 'Could not close meteorology ' //
     &                             'file ' // TRIM( PREVFILE )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ELSE
                            FILEOPEN = .FALSE.
                        END IF
                    END IF

                    PREVFILE = METFILE

                END IF

C.................  Set logical file name
                IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                    MESG = 'Could not set logical file name for ' //
     &                     'file ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Open the meteorology data file
                IF ( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                    MESG = 'Could not open meteorology file ' // 
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                ELSE
                    FILEOPEN = .TRUE.
                END IF

C.................  Read current meteorology file
                IF ( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                            JDATE, JTIME, TA ) ) THEN
                    MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                     ' from ' // TRIM( METFILE ) 
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
                IF ( .NOT. READ3( METNAME, MIXNAME, 1,
     &                            JDATE, JTIME, QV ) ) THEN
                    MESG = 'Could not read ' // TRIM( MIXNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

                IF ( .NOT. READ3( METNAME, PRESNAME, 1,
     &                            JDATE, JTIME, PRES ) ) THEN
                    MESG = 'Could not read ' // TRIM( PRESNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Apply ungridding matrix 
                CALL GRDFIPS( NSRC, CNTYSRC, TA, TASRC, .TRUE. )
                CALL GRDFIPS( NSRC, CNTYSRC, QV, QVSRC, .FALSE. )
                CALL GRDFIPS( NSRC, CNTYSRC, PRES, PRESSRC, .FALSE. )

C.................  Create hourly meteorology arrays by source
                CALL HOURMET( NSRC, AVG_TYPE, JDATE, JTIME, DAYBEGT,
     &             ALT_DATA, LDAYSAV, RH_STRHR, RH_ENDHR )

            END IF  ! check for using alternate data or day averaging

C.............  Make sure we've waited long enough to catch all time zones
            IF( POS > TSPREAD ) THEN

C.................  Adjust time step for 24-hour arrays
                ARRAYPOS = MOD( POS - TSPREAD, 24 )
                IF( ARRAYPOS == 0 ) ARRAYPOS = 24

C.................  If last day of month, process monthly averages
                CALL DAYMON( DDATE, MONTH, DAY )
                CALL DAYMON( DDATE + 1, TMPMNTH, DAY )

C.................  Estimate fuelmonth averaged monthly ref county temp and RH
                IF( TMPMNTH /= MONTH .AND. OTIME == 230000 ) THEN

C.....................  Averaging met data over no of days
                    DO K = 1,24
                        CALL AVGMET( NSRC,K )
                    ENDDO
                        
                    CALL AVG_REF_COUNTY_RH_TEMP( MONTH )
                    
                    MONOPEN = .TRUE.

                END IF
                
            END IF    ! time zone check

C.............  Estimate fuelmonth averaged episodic ref county temp and RH
            IF( T == NSTEPS .AND. .NOT. MONOPEN ) THEN

C.................  Averaging met data over no of days
                DO K = 1,24
                    CALL AVGMET( NSRC,K )
                ENDDO

                CALL AVG_REF_COUNTY_RH_TEMP( MONTH )

            END IF

C.............  Increment output time
            IF( POS > TSPREAD ) THEN
                CALL NEXTIME( DDATE, OTIME, 10000 )
            END IF

C.............  Increment loop time
            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, 10000 )

        END DO   !  End loop on hours of temperature files

C...........  Compute ref county min/max temperatures and averaged RH per
C             fuel month
        DO NR = 1, NREFC
          
            REFCOUNTY = MCREFIDX( NR,1 )

C.............  Choose month-specific fulemonth county
            L = FIND1FIRST( REFCOUNTY, NREFF,  FMREFSORT( :,1 ) )
            K = FIND1FIRST( REFCOUNTY, NFUELC, FMREFLIST( :,1 ) )

            IF( L < 0 .OR. K < 0 ) THEN
                WRITE( MESG,94010 ) 'ERROR: FUELMONTH input' //
     &              'file MUST contain reference county ', REFCOUNTY
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
           
            NMON = FMREFLIST( K, 2 )   ! no month available per ref county

C.............  Loop over months per ref. county
            N = 0
            NFMON = 0
            TMPMONTH = 0    ! prv tmp processing month 
            PRVFMONTH = 0   ! prv tmp fuelmonth 
            DO J = L, L + NMON - 1

                FUELMONTH = FMREFSORT( J,2 )    ! processing fuelmonth/county
                CURMONTH  = FMREFSORT( J,3 )    ! processing  current month per ref. county

C.................  Skip other months
                IF( CURMONTH <  SMONTH ) CYCLE
                IF( CURMONTH >  EMONTH ) CYCLE

                NFMON = NFMON + 1    ! county no of fuelmonth per refcouty for QA check

C.................  Store fuelmonth specific values into arrays
                IF( FUELMONTH /= PRVFMONTH ) THEN
                    IF( N > 0 ) THEN
                        PRNF = NF - N + 1
                        RHFUEL  ( NR,PRNF:NF ) = RHSUM / N
                        TKFUEL  ( NR,NF  ,:  ) = TKREFHR( : ) / N
                        MAXTFUEL( NR,PRNF:NF ) = MAXTEMP
                        MINTFUEL( NR,PRNF:NF ) = MINTEMP

c              print*,PRVFMONTH,NR,NF,PRNF,MAXTEMP,RHSUM/N,N,'REF1,,,'
c              print*,TKFUEL(NR,NF,:)

                        CALL WRTEMPROF( ODEV2, SDATE, AVG_TYPE,
     &                                  REFCOUNTY, TMPMONTH, PPTEMP,
     &                                  TKFUEL( NR,NF,: ) ) 

                    END IF

C.....................  initialize local variables
                    N     = 0
                    RHSUM = 0
                    TKREFHR = 0.0
                    MAXTEMP = BADVAL3
                    MINTEMP = -1*BADVAL3

                END IF

                NF = FIND1( CURMONTH, NFUEL, FUELIDX( NR,: ) )

                IF( RHFUEL( NR,NF ) > 0.0 ) N  = N + 1
                RHSUM = RHSUM + RHFUEL( NR,NF )
                MAXTEMP = MAX( MAXTEMP, MAXTFUEL( NR,NF ) )
                MINTEMP = MIN( MINTEMP, MINTFUEL( NR,NF ) )

                DO TT = 1,24
                    TKREFHR( TT ) = TKREFHR( TT ) + TKFUEL( NR,NF,TT )
                END DO
                
                TMPMONTH  = CURMONTH
                PRVFMONTH = FUELMONTH

c             print*,PRVFMONTH,NR,NF,PRNF,MAXTEMP,RHFUEL(NR,NF),N,'INV,,,'

            END DO

C...............  Store last fuelmonth specific values into arrays
            PRNF = NF - N + 1
            RHFUEL  ( NR,PRNF:NF ) = RHSUM / N
            TKFUEL  ( NR,NF  ,:  ) = TKREFHR( : ) / N  
            MAXTFUEL( NR,PRNF:NF ) = MAXTEMP
            MINTFUEL( NR,PRNF:NF ) = MINTEMP

c            print*,PRVFMONTH,NR,NF,PRNF,MAXTEMP,RHSUM/N,N,'REF2,,,'
c            print*,TKFUEL(NR,NF,:)

            CALL WRTEMPROF( ODEV2, SDATE, AVG_TYPE, REFCOUNTY,
     &                      TMPMONTH, PPTEMP, TKFUEL( NR,NF,: ) )

C..............  Check processing month is listed in the fuelmonth for each ref. county
            IF( NFMON < NFUEL ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Could not find the '//
     &              'modeling month(s) for reference county ',REFCOUNTY, 
     &              'in fuel month input file'
                CALL M3MSG2( MESG )
            END IF

        END DO   ! end of loop of reference couties
      
C.........  Exit if there was a problem with the meteorology files
        IF( EFLAG ) THEN
            MESG = 'Problem checking meteorology files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Initializing arrays
C.........  Allocate memory for storing meteorology profiles
        TKHOUR = 0.0
        RHHOUR = 0.0
        NDAYSRC = 0
        MAXTSRC = BADVAL3
        MINTSRC = -1*BADVAL3

        PREVFILE = ' '
        FILEOPEN = .TRUE.

C.........  Loop through days/hours of meteorology files
        DDATE = SDATE
        OTIME = 0
        WDATE = SDATE
        MDATE = SDATE
        
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9
        
        DO T = 1, NSTEPS

C.............  When new day...
            IF ( JDATE /= LDATE ) THEN

C.................  Write message for day of week and date
                DAY = WKDAY( JDATE )
                MESG = 'Processing ' // DAYS( DAY ) // ' ' // 
     &                 MMDDYY( JDATE )
                CALL M3MSG2( MESG )

C.................  Set start and end hours of day for all sources for MOVES model
                CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, .TRUE.,
     &                         DAYBEGT, DAYENDT )

            END IF

C.............  Determine input file for this hour
            POS = T
            
C.............  Get file number for current iteration
            FILENUM = METDAYS( POS )

            IF( FILENUM <= 0 ) THEN
                ALT_DATA = .TRUE.
            ELSE
                ALT_DATA = .FALSE.
            END IF

C.............  Skip file opening when not doing day averaging and using alternate data
            IF( .NOT. ALT_DATA .OR. DAYAVER ) THEN
            
C.................  Get file name
                METFILE = METLIST( ABS( FILENUM ) )

C.................  Close previous file if needed
                IF( METFILE .NE. PREVFILE ) THEN
                    IF( FILEOPEN ) THEN
                        IF( .NOT. CLOSE3( METNAME ) ) THEN
                            MESG = 'Could not close meteorology ' //
     &                             'file ' // TRIM( PREVFILE )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ELSE
                            FILEOPEN = .FALSE.
                        END IF
                    END IF

                    PREVFILE = METFILE

                END IF

C.................  Set logical file name
                IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                    MESG = 'Could not set logical file name for ' //
     &                     'file ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Open the meteorology data file
                IF ( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                    MESG = 'Could not open meteorology file ' // 
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                ELSE
                    FILEOPEN = .TRUE.
                END IF

C.................  Reset read date when using alternate data while processing daily (no averaging)
                IF( ALT_DATA ) THEN
                    IF( .NOT. DESC3( METNAME ) ) THEN
                        MESG = 'Could not get description of ' //
     &                         'meteorology file ' // TRIM( METFILE )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    ELSE
                        IF( SDATE3D < JDATE ) THEN
                            RDATE = JDATE - 1
                        ELSE
                            RDATE = JDATE + 1
                        END IF
                    END IF

                ELSE
                    RDATE = JDATE
                END IF
            
C.................  Read current meteorology file
                IF ( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                            RDATE, JTIME, TA ) ) THEN
                    MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                     ' from ' // TRIM( METFILE ) 
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF
                IF ( .NOT. READ3( METNAME, MIXNAME, 1,
     &                            RDATE, JTIME, QV ) ) THEN
                    MESG = 'Could not read ' // TRIM( MIXNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF

                IF ( .NOT. READ3( METNAME, PRESNAME, 1,
     &                            RDATE, JTIME, PRES ) ) THEN
                    MESG = 'Could not read ' // TRIM( PRESNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                END IF

C.................  Apply ungridding matrix 
                CALL GRDFIPS( NSRC, CNTYSRC, TA, TASRC, .TRUE. )
                CALL GRDFIPS( NSRC, CNTYSRC, QV, QVSRC, .FALSE. )
                CALL GRDFIPS( NSRC, CNTYSRC, PRES, PRESSRC, .FALSE. )

C.................  Create hourly meteorology arrays by source
                CALL HOURMET( NSRC, AVG_TYPE, JDATE, JTIME, DAYBEGT,
     &             ALT_DATA, LDAYSAV, RH_STRHR, RH_ENDHR )

            ELSE
                IF( OTIME == 230000 ) THEN
                    MESG = 'NOTE: Missing meteorology file on '//
     &                      MMDDYY( JDATE )
                    CALL M3MSG2( MESG )
                END IF

            END IF  ! check for using alternate data or day averaging

C.............  Make sure we've waited long enough to catch all time zones
            IF( POS > TSPREAD ) THEN

C.................  Adjust time step for 24-hour arrays
                ARRAYPOS = MOD( POS - TSPREAD, 24 )
                IF( ARRAYPOS == 0 ) ARRAYPOS = 24

C.................  Process daily averages
                IF( DAYAVER ) THEN

C.....................  Average temperatures across county group                
                    CALL AVGMET( NSRC, ARRAYPOS )

C.....................  Write averaged daily county temp and RH to file
                    IF( OTIME == 230000 ) THEN
                        CALL WRAVGMET( NSRC, ODEV1, DDATE )
                    END IF

                END IF

C.................  If last day of month, process monthly averages
                CALL DAYMON( DDATE, MONTH, DAY )
                CALL DAYMON( DDATE + 1, TMPMNTH, DAY )
            
                IF( MONAVER .AND. TMPMNTH /= MONTH ) THEN

C.....................  Average temperatures across county group 
                    CALL AVGMET( NSRC, ARRAYPOS )

C.....................  Write averaged monthly county temp and RH to file
C                       only output when more than 15 days are processed for monthly avg
                    IF( OTIME == 230000 .AND. POS > 15*24 ) THEN
                        CALL WRAVGMET( NSRC, ODEV1, DDATE )
                    END IF
                  
                END IF

            END IF    ! time zone check

C.............  Output episode averaged temperatures
C.............  Write averaged episodic county temp and RH to file
            IF( T == NSTEPS .AND. EPIAVER ) THEN
C.................  Average temperatures across county group 
                DO K = 1, 24
                    CALL AVGMET( NSRC, K )
                END DO

                CALL WRAVGMET( NSRC, ODEV1, DDATE )
            END IF

C.............  Increment output time
            IF( POS > TSPREAD ) THEN
                CALL NEXTIME( DDATE, OTIME, 10000 )
            END IF

C.............  Increment loop time
            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, 10000 )

        END DO   !  End loop on hours of temperature files
 
C......... End program successfully

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
 
C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Formatted file I/O formats............ 93xxx
 
93010   FORMAT( I8, 1X, F13.5, 1X, F13.5, 1X, I8, 1X, A )
 
C...........   Internal buffering fosrmats............ 94xxx

94000   FORMAT( A )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, :, 2( I7.7, :, A ) )

94030   FORMAT( A, I5 )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )
     
94070   FORMAT( A, F5.1, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram estimates ref. county level 
C               averaged RH and min/max Temperatures over fuelmonth

            SUBROUTINE AVG_REF_COUNTY_RH_TEMP( MONTH )

C.............  local argument
            INTEGER, INTENT( IN ) :: MONTH

C---------------------------------------------------------------------- 
C.............  Loop over sources
            IC = 0
            PRCOUNTY = 0
            DO S = 1, NSRC
                        
                INVCOUNTY = MCREFSORT( S,1 )
                REFCOUNTY = MCREFSORT( S,2 )

                N = 0
                RHSUM = 0.0
C.................  averaging RH per inventory county
                DO TT = 1,24
                    IF( RHHOUR( S,TT ) == 0.0 ) CYCLE
                    N = N + 1
                    RHSUM = RHSUM + RHHOUR( S,TT )
                END DO

                RHAVG   = RHSUM / N
                MAXTEMP = MAXTSRC( S )
                MINTEMP = MINTSRC( S )

c           print*,INVCOUNTY,REFCOUNTY,MAXTEMP,MINTEMP,RHAVG,N,RHSUM,'RHAVG'
C.................  Calculation monthly max/min temp and avg RH
C                   per ref. county

C.................  Averaging RH for ref. county
                IF( PRCOUNTY /= REFCOUNTY ) THEN
                     IF( IC > 0 ) THEN
                          NR = FIND1( PRCOUNTY,NREFC, MCREFIDX( :,1 ) )
                          NF = FIND1( MONTH, NFUEL, FUELIDX( NR,: ) )

                          RHFUEL( NR,NF )   = RHREFSUM / IC
                          TKFUEL( NR,NF,: ) = TKREFHR( : ) / IC
                          MAXTFUEL( NR,NF ) = MAXTREF
                          MINTFUEL( NR,NF ) = MINTREF

c           print*,prcounty,NR,NF,MAXTREF,MINTREF,IC,RHREFSUM/IC,'REF1,,'
                     END IF
                            
                     IC = 0
                     MAXTREF = BADVAL3
                     MINTREF = -1*BADVAL3
                     RHREFSUM = 0.0
                     TKREFHR  = 0.0

                END IF
                        
                IC = IC + 1
                RHREFSUM = RHREFSUM + RHAVG

                DO TT = 1, 24
                    TKREFHR( TT ) = TKREFHR( TT ) + TKHOUR( S,TT )
                END DO

                MAXTREF = MAX( MAXTREF, MAXTEMP )
                MINTREF = MIN( MINTREF, MINTEMP )

c         print*,refcounty,NR,NF,MAXTREF,MINTREF,IC,RHAVG,'INV,,.'
                
                PRCOUNTY = REFCOUNTY
                        
            END DO

            NR = FIND1( PRCOUNTY,NREFC, MCREFIDX( :,1 ) )
            NF = FIND1( MONTH, NFUEL, FUELIDX( NR,: ) )

            RHFUEL( NR,NF )   = RHREFSUM / IC
            TKFUEL( NR,NF,: ) = TKREFHR( : ) / IC
            MAXTFUEL( NR,NF ) = MAXTREF
            MINTFUEL( NR,NF ) = MINTREF

C.............  reinitializing local arrays
            NDAYSRC = 0
            TKHOUR = 0.0
            RHHOUR = 0.0
            MAXTSRC = BADVAL3
            MINTSRC = -1*BADVAL3
       
c         print*,prcounty,NR,NF,MAXTREF,MINTREF,IC,RHREFSUM/IC,'REF2,,'

            END SUBROUTINE AVG_REF_COUNTY_RH_TEMP
               
        END PROGRAM MET4MOVES
