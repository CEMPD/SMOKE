
        PROGRAM GENTPRO

C***********************************************************************
C
C  DESCRIPTION:
C       Temporal profile generator.
C
C       The temporal profile generator (GenTPRO) calculates county
C       average meteorology data to estimate temporal profiles of
C       emissions using algorithms that relate meteorology to emissions
C       fluxes. GenTPRO reads in hourly meteorology data from the
C       meteorology-chemistry interface processor (MCIP) and a gridded
C       spatial surrogate (commonly used in SMOKE) to produce temporal
C       profiles and cross-reference data in a comma delimited (CSV)
C       format and a temporal data binary netCDF file.
C
C       More specifically, GenTPRO creates county-based temporal
C       profiles using RWC (residential wood combustion) algorithms,
C       NH3 (agricultural ammonia) equations, and MET (meteorological)
C       variables based on meteorological variables from gridded
C       meteorology data. Temperatures and wind speed (or aerodynamic resistance)
C       can be averaged across counties and over different time periods.
C
C       GenTPRO produces two CSV files for input to SMOKE. The temporal
C       profile file contains scalars that allocate emission inventory
C       data to specific temporal periods. A temporal cross-reference
C       file relates FIPS codes for the inventory administrative units
C       (i.e. counties) and surrogate numbers.  Every administrative
C       unit with non-zero inventory data contained in the modeling grid
C       will have an entry in the cross-reference file. If the hour-of-
C       year option for the AGNH3 and MET profile types are selected,
C       GenTPRO outputs a netCDF file of temporal data.
C
C  PRECONDITIONS REQUIRED: none
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     2/11: Created  by B.H. Baek
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
        USE MODINFO,  ONLY: CRL, CATEGORY, NSRC

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NGRID, NCOLS, NROWS, COORD, GDTYP,
     &                     P_ALP, P_BET, P_GAM, XCENT, YCENT, NGRID

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

C.........  Force explicit declaration of all variables
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        CHARACTER(14)   MMDDYY
        CHARACTER(16)   PROMPTMFILE
        CHARACTER(16)   VERCHAR
        INTEGER         FIND1, FINDC
        INTEGER         FIND1FIRST
        INTEGER         GETIFDSC
        INTEGER         GETFLINE
        INTEGER         GETEFILE
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        INTEGER         WKDAY
        INTEGER         ENVINT
        REAL            ENVREAL
        REAL            YR2DAY
        REAL            STR2REAL
        LOGICAL         BLKORCMT
        LOGICAL         ENVYN
        LOGICAL         DSCM3GRD
        LOGICAL         INTLIST
        LOGICAL         SETENVVAR
        LOGICAL         STRLIST

        EXTERNAL     CRLF, DSCM3GRD, GETIFDSC, GETFLINE, ENVINT, FIND1
     &               ENVREAL, INDEX1, MMDDYY, PROMPTFFILE, PROMPTMFILE,
     &               SECSDIFF, SETENVVAR, WKDAY, GETEFILE, INTLIST, 
     &               FIND1FIRST, STRLIST, STR2INT, BLKORCMT, VERCHAR,
     &               YR2DAY, ENVYN, STR2REAL, FINDC

C.....  Define temporal profile type constants for enumeration
        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv4.8_Jun2020$' ! CVS release tag

        INTEGER, PARAMETER :: MXVAR  = 100
        INTEGER, PARAMETER :: MXSEG  = 16
        INTEGER, PARAMETER :: NMONTH = 12
        INTEGER, PARAMETER :: MXDAYS = 366
        INTEGER, PARAMETER :: NTPRO  = 14

C...........   real arrays
        REAL   , ALLOCATABLE :: TA( : )            !  one layer of temperature
        REAL   , ALLOCATABLE :: VARTMP( : )        !  one layer of 2nd variable for NH3 profile method
        REAL   , ALLOCATABLE :: TASRC( : )         !  averaged gridded temperature by FIPS(source)
        REAL   , ALLOCATABLE :: VARSRC( : )        !  averaged gridded 2nd variable by FIPS(source)
        REAL   , ALLOCATABLE :: TOTSRC( :,:,: )    !  hourly daily/month total
        REAL   , ALLOCATABLE :: HRLSRC( :,: )      !  hourly values
        REAL   , ALLOCATABLE :: DAYSRC( :,: )      !  daily values
        REAL   , ALLOCATABLE :: MONSRC( :,: )      !  monthly total
        REAL   , ALLOCATABLE :: ANNSRC( : )        !  annual total value
        REAL   , ALLOCATABLE :: TMPDSRC( : )       !  tmp daily total
        REAL   , ALLOCATABLE :: TMPMSRC( : )       !  tmp monthly total
        REAL   , ALLOCATABLE :: MINTEMP ( : )      !  min temp values by source (FIPS)
        REAL   , ALLOCATABLE :: RWC_TEMP( : )      !  county-min temp values for RWC
        REAL   , ALLOCATABLE :: PROF_MON(:,:)      ! Array of monthly profile fractional data values
        REAL   , ALLOCATABLE :: PROF_DAY(:,:)      ! Array of daily profile fractional data values

C...........   integer arrays
        INTEGER, ALLOCATABLE :: DT      ( : )      ! daily total starting hour
        INTEGER, ALLOCATABLE :: MT      ( : )      ! monthly total starting hour
        INTEGER, ALLOCATABLE :: DAYBEGT ( : )      ! current daily start time HHMMSS
        INTEGER, ALLOCATABLE :: DAYENDT ( : )      ! current daily end time HHMMSS
        INTEGER, ALLOCATABLE :: PRVBEGT ( : )      ! previous day start time HHMMSS
        INTEGER, ALLOCATABLE :: PRVENDT ( : )      ! previous day end time HHMMSS
        INTEGER, ALLOCATABLE :: TZONES  ( : )      ! county-specific time zones
        INTEGER, ALLOCATABLE :: METDAYS ( : )      ! dimension: nsteps in episode,
        INTEGER, ALLOCATABLE :: PROCDAYS( : )      ! no of processing hours 
        INTEGER, ALLOCATABLE :: SRGIDS  ( : )      ! list of surrogates
        INTEGER, ALLOCATABLE :: SRGSTA  ( : )      ! list of state in surrogates
        INTEGER, ALLOCATABLE :: INDXREF ( : )      ! Index of matched xref entries 
        INTEGER, ALLOCATABLE :: MATCHED ( :,: )    ! FIPS/SCC/POL matched source

C...........  character arrays
        CHARACTER(16)                      SEGMENT( MXSEG )
        CHARACTER(256)    , ALLOCATABLE :: METLIST( : )       ! listing of met file names
        CHARACTER(SCCLEN3), ALLOCATABLE :: SCCLIST( : )       ! listing of SCCs
        CHARACTER(256)    , ALLOCATABLE :: CSCCFIP( : )       ! tmp FIPS/SCC x-ref entries

C...........   Parameter array
        CHARACTER(16), PARAMETER :: TPROTYPES( NTPRO ) =
     &      ( / 'MONTHLY         ','DAILY           ','HOURLY          ',
     &          'WEEKLY          ','WEEKEND         ','WEEKDAY         ',
     &          'MONDAY          ','TUESDAY         ','WEDNESDAY       ',
     &          'THURSDAY        ','FRIDAY          ','SATURDAY        ',
     &          'SUNDAY          ','ALLDAY          ' / )
        
C...........   File units and logical names:
        INTEGER      CDEV  ! unit number for co/st/cy file
        INTEGER      DDEV  ! unit number for daily group file
        INTEGER      EDEV  ! unit number for episode group file
        INTEGER      GDEV  ! tmp unit number for individual surrogate file
        INTEGER      IDEV  ! tmp unit number if ENAME is map file
        INTEGER      LDEV  ! unit number for log file
        INTEGER      MDEV  ! unit number for monthly group file
        INTEGER      PDEV  ! unit number for speeds summary file (SPDSUM)
        INTEGER      QDEV  ! unit number for surrogate(s) file
        INTEGER      RDEV  ! unit number for county-specific temp file for RWC
        INTEGER      SDEV  ! unit number for surrogate(s) file
        INTEGER      TDEV  ! unit number for meteorology list file
        INTEGER      WDEV  ! unit number for weekly group file
        INTEGER      ODEV  ! unit number for tmp surrogate file
        INTEGER      XIDEV ! unit number for x-ref  input file
        INTEGER      XODEV ! unit number for x-ref output file
        INTEGER      MODEV ! unit number for met-based monthly temporal profile output file
        INTEGER      DODEV ! unit number for met-based daily   temporal profile output file
        INTEGER      HODEV ! unit number for met-based hourly  temporal profile output file

        CHARACTER(16)    METNAME  ! logical name for meteorology files
        CHARACTER(16) :: HNAME= 'TPRO_HOUR'   ! logical name for hourly output in ncf format

C...........   Other local variables:
        INTEGER    DD, I, IC, J, K, LL, L, L0, L1, L2, L3   ! Counters and pointers
        INTEGER    MM, NN, N, NP, NX, NRH, NR, NS, S, T, NT, PT, TT, T2, V  ! Counters and pointers

        INTEGER    EPI_SDATE   ! episode start date from E.V. (YYYYDDD)
        INTEGER    EPI_STIME   ! episode start time from E.V. (HHMMSS)
        INTEGER    EPI_RUNLEN  ! episode duration   from E.V. (HHMMSS)
        INTEGER    EPI_NSTEPS  ! episode number of time steps
        INTEGER    EPI_EDATE   ! episode ending date based on ERUNLEN
        INTEGER    EPI_ETIME   ! episode ending time based on ERUNLEN

        INTEGER    DAY         ! tmp day of week number
        INTEGER    DST         ! tmp daylight saving time
        INTEGER    EDATE       ! ending input date counter (YYYYDDD) in GMT
        INTEGER    ETIME       ! ending input time counter (HHMMSS)  in GMT
        INTEGER    FIPS        ! tmp inventory county
        INTEGER    FILENUM     ! file number of current meteorology file
        INTEGER    IOS         ! temporary I/O status
        INTEGER    IREC        ! temporary input line number
        INTEGER    IFIP,NTP    ! temporary FIPS code
        INTEGER    ISRGFIP     ! temporary surrogate FIPS code
        INTEGER    HOURIDX     ! current hour of the day
        INTEGER    JDATE       ! input date counter (YYYYDDD) in GMT
        INTEGER    JTIME       ! input time counter (HHMMSS)  in GMT
        INTEGER    LDATE       ! date from previous loop iteration
        INTEGER    SMONTH      ! month of start date
        INTEGER    EMONTH      ! month of end date
        INTEGER    MONTH       ! tmp month
        INTEGER    TYEAR       ! tmp year 
        INTEGER    TDAY        ! tmp day of month
        INTEGER    METNGRID    ! no. grid cells in met data
        INTEGER    NLINES      ! no. lines in met list file
        INTEGER    NDAY        ! no. processing days 
        INTEGER    NMON        ! no. processing month 
        INTEGER    NVAR        ! no. met variables
        INTEGER    NSRG        ! no. surrogates
        INTEGER    NSCC        ! no. processing SCCs
        INTEGER    ISTA        ! current state number
        INTEGER    PSTA        ! previous State number
        INTEGER    NSTA        ! no. of unique State number
        INTEGER    NCNTY       ! no. of counties for RWC method
        INTEGER    NMATCH      ! no. of matched temporal x-ref input file by FIPS/SCC
        INTEGER    MXTREF      ! no. matched temporal x-ref entries 
        INTEGER    MXLINE      ! no. temporal x-ref input file
        INTEGER    NSTEPS      ! number of time steps to process temperature data
        INTEGER    POS         ! position in time step loop
        INTEGER    REFCOUNTY   ! ref. county FIPS code
        INTEGER    INVCOUNTY   ! inv. county FIPS code
        INTEGER    PRCOUNTY    ! previous ref. county
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
        INTEGER    NWARN       ! no of warning messgaes

        INTEGER, SAVE :: MXWARN        ! maximum no of warning messgaes

        REAL       DTEMP               ! RWC default temp (=50.0) 
        REAL       TEMPVAL             ! tmp variable value
        REAL       SLOPE               ! RWC linear euqation slope 
        REAL       CONST               ! RWC linear equation constant 

        LOGICAL :: BASHFLAG = .FALSE.  !  true: processing NH3 option using Bash Equation
        LOGICAL :: EFLAG    = .FALSE.  !  true: error found
        LOGICAL :: COMPLETE = .FALSE.  !  true: program successful complete
        LOGICAL :: GRID_ERR = .FALSE.  !  true: error found in grid settings
        LOGICAL :: CFLAG    = .FALSE.  !  true: Use county-specific min temp setting for RWC eq
        LOGICAL :: NH3FLAG  = .FALSE.  !  true: processing NH3 profile method
        LOGICAL :: MONAVER  = .FALSE.  !  true: monthly averaging
        LOGICAL :: DAYAVER  = .FALSE.  !  true: daily averaging
        LOGICAL :: HOURAVER = .FALSE.  !  true: hourly averaging
        LOGICAL :: OFLAG    = .FALSE.  !  true: ungridding is 0 for some srcs
        LOGICAL :: ZFLAG    = .FALSE.  !  true: use AZ's new equations
        LOGICAL :: FILEOPEN = .FALSE.  !  true: met file is open
        LOGICAL :: FND_DATA = .FALSE.  !  true: found met data for this hour
        LOGICAL :: ALT_DATA = .FALSE.  !  true: using alternate data for this hour

        CHARACTER(SCCLEN3) CSCC        !  SCC code
        CHARACTER(FIPLEN3) CFIPS       !  FIPS code
        CHARACTER(FIPLEN3) TPROID      !  tpro id
        CHARACTER(1000)    CSCCLIST    !  tmp SCC list line buffer 
        CHARACTER(16)      CPOL,TPRO   !  Pollutant code, Profile types 
        CHARACTER(SCCLEN3) TSCC        !  tmp SCC code
        CHARACTER(IOVLEN3) COORUNIT    !  coordinate system projection units
        CHARACTER(IODLEN3) GDESC       !  grid description
        CHARACTER(16)      SRG_CNTRY   !  surrogate country
        CHARACTER(256)     LINE        !  line buffer
        CHARACTER(256)     FULLNAME    !  full file name
        CHARACTER(512)     METFILE     !  tmp physical file name
        CHARACTER(512)     PREVFILE    !  previous physical file name
        CHARACTER(IOVLEN3) PROF_METHOD !  Profile method name
        CHARACTER(IOVLEN3) TPRO_TYPE   !  Temproal profile type name
        CHARACTER(IOVLEN3) TVARNAME    !  temperature variable name
        CHARACTER(IOVLEN3) VAR_NAME    !  NH3 2nd variable name
        CHARACTER(200)     TEMPDIR     !  directory for output files
        CHARACTER(300)     MESG        !  message buffer

        CHARACTER(16), PARAMETER :: PROGNAME = 'GENTPRO'  ! program name
C***********************************************************************
C   begin body of program GENTPRO

C......... Initialize and obtain logical unit for log file.
        LDEV = INIT3()

C.........  Write out copyright, version, web address, HDRer info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Obtain parameter settings from the environment...
C.........  Get file name for country, state, and county file, with time zones
        CDEV = PROMPTFFILE(
     &         'Enter logical name for Country/State/County file', 
     &         .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  Open input surrogates file
        QDEV = PROMPTFFILE(
     &         'Enter logical name for Surrogate Description file',
     &         .TRUE., .TRUE., 'SRGDESC', PROGNAME )

C.........  Open met list file
        TDEV= PROMPTFFILE(
     &         'Enter logical name for Meteorology list file',
     &         .TRUE., .TRUE., 'METLIST', PROGNAME )

C.........  Open temporal cross-reference INPUT file (TREF_IN)
        XIDEV = PROMPTFFILE(
     &          'Enter logical name for temporal x-reference file',
     &          .TRUE., .TRUE., 'TREF_IN', PROGNAME )

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET , ' ', 100, I )
            
C.........  Allocate arrays
        ALLOCATE( SCCLIST( MXVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCLIST', PROGNAME )
        ALLOCATE( SRGIDS( MXVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGIDS', PROGNAME )
        SCCLIST = ' '
        SRGIDS  = 0

C.........  Get a list of SCCs for meteorology processing.
        MESG = 'Specifies a list of SCCs'
        CALL ENVSTR( 'SCC_LIST', MESG, '', CSCCLIST, IOS )
        CALL PARSLINE( CSCCLIST, MXVAR, SCCLIST )

        IF( CSCCLIST == ' ' ) THEN
            MESG = 'ERROR: MUST define a list of SCCs [SCC_LIST]'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NSCC = 0
        DO I = 1, MXVAR
           IF( SCCLIST( I ) /= ' ' ) NSCC = NSCC + 1
        END DO

C.........  Get name of surrogate IDs to use
        MESG = 'Specifies a list of spatial surrogate IDs'
        IF( .NOT. INTLIST( 'SRG_LIST',MESG,MXVAR,NSRG,SRGIDS ) ) THEN
            MESG = 'Could not read list of surrogate IDs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get the type of temporal profiles to produce from the environment
        MESG = 'Specifies the type of temporal profiles to produce'
        CALL ENVSTR( 'TPRO_OUTPUT', MESG, ' ', TPRO_TYPE, IOS )
        CALL UPCASE( TPRO_TYPE )

        MESG = 'ERROR: MUST define TPRO_OUTPUT'
        IF( TPRO_TYPE == ' ' ) CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Obtain profile generation method: 'RWC', 'RC_NH3', 'BASH_NH3', or 'MET'
        MESG ='Specifies profile method for meteorology processing'
        CALL ENVSTR( 'PROFILE_METHOD', MESG, 'MET', PROF_METHOD, IOS )
        CALL UPCASE( PROF_METHOD )

        MESG = 'ERROR: ' // TRIM( PROF_METHOD ) // ' profile method'
     &         // ' you defined is not applicable' 
        IF( .NOT. ( PROF_METHOD == 'RC_NH3' .OR. PROF_METHOD == 'BASH_NH3' .OR.
     &      PROF_METHOD == 'MET' .OR. PROF_METHOD == 'RWC' ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get the name of the temperature variable.
        MESG = 'Specifies meteorological variable name for ' //
     &         TRIM( PROF_METHOD ) // ' profile method'
        CALL ENVSTR( 'TEMP_VAR', MESG, ' ',TVARNAME, IOS )
        CALL UPCASE( TVARNAME )

        IF( TVARNAME == ' ' ) THEN
            MESG = 'ERROR: MUST specifies meteorological variable ' //
     &          'name "TEMP_VAR" for '//TRIM( PROF_METHOD )//
     &          ' profile method'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set default names for additional variables
        IF( PROF_METHOD == 'RC_NH3' ) THEN

            NH3FLAG  = .TRUE.

C.............  Get the name of the wind speed variable.
            MESG = 'Specifies wind speed variable name for ' //
     &             TRIM( PROF_METHOD ) // ' profile method'
            CALL ENVSTR( 'WSPEED_VAR', MESG, ' ',VAR_NAME, IOS )
            CALL UPCASE( VAR_NAME )

           IF( VAR_NAME == ' ' ) THEN
               MESG = 'ERROR: MUST specifies wind speed variable ' //
     &             'name "WSPEED_VAR" for '//TRIM( PROF_METHOD )//
     &             ' profile method'
               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           ELSE
               MESG = 'NOTE: Wind speed variable ('//TRIM( VAR_NAME )//
     &            ') is chosen for "' // TRIM( PROF_METHOD ) //
     &            '" profile method'
               CALL M3MESG( MESG )
           END IF

        ELSE IF( PROF_METHOD == 'BASH_NH3' ) THEN

            NH3FLAG  = .TRUE.
            BASHFLAG = .TRUE.

C.............  Check name of temperature variable
            IF( TVARNAME /= 'TEMPG' ) THEN
                MESG = 'WARNING: Bash Equation is based on surface ' //
     &              'temperature variable.'
                CALL M3MESG( MESG )
            END IF

C.............  Get the name of the aerodynamic resistance  variable.
            MESG = 'Specifies aerodynamic resistance variable name for ' //
     &             TRIM( PROF_METHOD ) // ' profile method'
            CALL ENVSTR( 'AERO_RESISTANCE_VAR', MESG, ' ',VAR_NAME, IOS )
            CALL UPCASE( VAR_NAME )

           IF( VAR_NAME == ' ' ) THEN
               MESG = 'ERROR: MUST specifies aerodynamic resistance ' //
     &             'variable name "AERO_RESISTANCE_VAR" for '//
     &             TRIM( PROF_METHOD )// ' profile method'
               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           ELSE
               MESG = 'NOTE: Aerodynamic resistance variable ('
     &             // TRIM( VAR_NAME ) // ') is chosen for "' 
     &             // TRIM( PROF_METHOD ) // '" profile method'
               CALL M3MESG( MESG )
           END IF

        ELSE
           VAR_NAME = '' 

        END IF

C.........  Defines output types based on profile method
        IF( TPRO_TYPE == 'MONTHLY' ) MONAVER  = .TRUE.
        IF( TPRO_TYPE == 'HOURLY'  ) HOURAVER = .TRUE.
        IF( TPRO_TYPE == 'DAILY'   ) THEN
            MONAVER  = .TRUE.
            DAYAVER  = .TRUE.
        END IF
        IF( TPRO_TYPE == 'ALL' ) THEN
            MONAVER  = .TRUE.
            DAYAVER  = .TRUE.
            HOURAVER = .TRUE.
        END IF
C.........  Determine optional linear equation for RWC profile calculation
        IF( PROF_METHOD == 'RWC' ) THEN

            MESG = 'Enter A for RWC equation: y = Ax + B'
            SLOPE = ENVREAL( 'SLOPE', MESG, 0.79, IOS )
            MESG = 'Enter B for RWC equation: y = Ax + B'
            CONST = ENVREAL( 'CONSTANT', MESG ,42.12, IOS )

            MESG = 'Enter default minimum temp for RWC method'
            DTEMP= ENVREAL( 'DEFAULT_TEMP_RWC', MESG ,50.0, IOS )

            MESG = 'Use county-specific min temperature for RWC'
            CFLAG = ENVYN( 'RWC_COUNTY_TEMP_YN', MESG, .FALSE., IOS )

            MESG = 'Use the alternative RWC equation'
            ZFLAG = ENVYN( 'RWC_ALT_EQ_YN', MESG, .TRUE., IOS )

        END IF

C.........  Error if hourly output setting for RWC profile method.
        IF( PROF_METHOD == 'RWC' ) THEN
            MESG = 'ERROR: Profile Method "'//TRIM( PROF_METHOD )//'" can not' 
     &          // ' process "'//TRIM(TPRO_TYPE)//'" temporal output option'
            IF( TPRO_TYPE == 'HOURLY' .OR. TPRO_TYPE == 'ALL' ) 
     &          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            HOURAVER = .FALSE.   ! override the setting

C.........  Error if daily or monthly output setting for AGNH3 profile method.
        ELSE IF( NH3FLAG ) THEN
            MESG = 'ERROR: Profile Method "'//TRIM( PROF_METHOD )//'" can not' 
     &          // ' process "'//TRIM(TPRO_TYPE)//'" temporal output option'
     &         //CRLF()//BLANK10//':: Set TPRO_OUTPUT to HOURLY'
            IF ( TPRO_TYPE == 'ALL' )   CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
            IF ( DAYAVER .OR. MONAVER ) CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 

        END IF

C.........  Open and obtain logical units for output files
C.........  Open temporal cross-reference OUTPUT file (TREF_OUT)
        MESG = 'Enter logical name for temporal cross-reference ' //
     &         'output file'
        XODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'TREF_OUT', PROGNAME )

C.........  Open output file for monthly temporal profiles
        IF( MONAVER ) THEN
            MESG = 'Enter logical name for TPRO_MON temporal profile '//
     &             'output file'
            MODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'TPRO_MON',
     &                           PROGNAME )
        END IF

C.........  Open output file for daily temporal profiles
        IF( DAYAVER ) THEN
            MESG = 'Enter logical name for TPRO_DAY temporal profile '//
     &             'output file'
            DODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'TPRO_DAY',
     &                           PROGNAME )
        END IF

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Read the surrogates HDRer and initialize the grid description
C.........  Open/read SRGDESC file
        MESG = 'Country of surrogate files'
        CALL ENVSTR( 'SRG_COUNTRY', MESG, 'USA', SRG_CNTRY, IOS )

        CALL RDSRGDESC( QDEV )

C.........  Open output file for temporary combined surrogate file
        FULLNAME = 'TMP_COMBINED_SRG.txt'
        IF( .NOT. SETENVVAR( 'TMP_SRG_FILE', FULLNAME )) THEN
             MESG = 'Could not set logical file name of file '
     &              // TRIM( FULLNAME )
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ODEV = GETEFILE( 'TMP_SRG_FILE', .FALSE., .TRUE., PROGNAME )

C.........  Prompt for and open input surrogate file(s)...
        DO I = 1, NSRG

            L = 0
            DO K = 1,NTSRGDSC
                CALL UPCASE( SRGFREG( K ) )
                IF( SRG_CNTRY == SRGFREG( K ) ) THEN
                    IF( SRGIDS( I )  == SRGFCOD( K ) ) L = K
                END IF
            END DO

            IF( L < 1 ) CYCLE

            CALL GETENV( 'SRGPRO_PATH', TEMPDIR )
            WRITE( FULLNAME, '(3A)' ) TRIM( TEMPDIR ),'/',
     &                                TRIM( SRGFNAM( L ) )

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
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D,
     &       COORUNIT, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,
     &       YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &       YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )
        IF( EFLAG ) THEN
            MESG = 'Problem with gridded input data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for and read the gridding surrogates file,
C           extracting data for a subgrid, if necessary
        CALL M3MSG2( 'Processing gridding surrogate(s)...' )

        CALL RDSRG( .FALSE., SDEV, SRGFMT, SRGNROWS, SRGNCOLS )

        MESG = 'NOTE: A list of surrogates used in the modeling'
        CALL M3MESG( MESG )

        DO I = 1, NSRG
            WRITE( MESG,94010 ) 'Surrogate ID ::', SRGLIST( I )
            CALL M3MESG( MESG )
        END DO

C.........  Read region codes file
        CALL RDSTCY( CDEV, NSRGFIPS, SRGFIPS )

C.........  Allocate arrays for county time zone
        ALLOCATE( TZONES( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONES', PROGNAME )
        TZONES = 0

C.........  Assign time zone to inventory counties
        DO I = 1, NSRGFIPS
            J = FINDC( SRGFIPS( I ), NCOUNTY, CNTYCOD )
            IF( J < 1 ) THEN
                MESG = 'ERROR: Could not find time zone for county '
     &               // SRGFIPS( I ) // ' from COSTCY file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                TZONES( I ) = CNTYTZON( J )
            END IF
        END DO

C.........  count no of states in modeling modain
        NSTA = 0
        PSTA = 0
        DO I = 1, NSRGFIPS
            ISTA = STR2INT( SRGFIPS( I )( 8:9 ) )
            IF( ISTA /= PSTA ) THEN
                NSTA = NSTA + 1
                PSTA = ISTA
            END IF
        END DO

C.........  Allocate array
        ALLOCATE( SRGSTA( NSTA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGSTA', PROGNAME )
        SRGSTA = 0

        NS = 0
        PSTA = 0
        DO I = 1, NSRGFIPS
            ISTA = STR2INT( SRGFIPS( I )( 8:9 ) )
            IF( ISTA /= PSTA ) THEN
                NS = NS + 1
                SRGSTA( NS ) = ISTA
                PSTA = ISTA
            END IF
        END DO

C.........  Define total number of sources (FIPS * SCC)
        NSRC = NSRGFIPS * NSCC

C.........  Allocate arrays for county-specific temperature settings for RWC
        ALLOCATE( RWC_TEMP( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RWC_TEMP', PROGNAME )
        RWC_TEMP = DTEMP

C.........  Open county-specific Temperature setting for RWC
        IF( CFLAG ) THEN

            MESG='Enter logical name for County-specific Temperature '//
     &           'Input file for RWC method'
            RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,'RWC_COUNTY_TEMP',
     &             PROGNAME )

            NLINES= GETFLINE( RDEV,'County-specific Temperatures file' )
             
C..............  Store county-specific temperatures
            IREC = 0
            NCNTY = 0
            DO I = 1, NLINES

                READ( RDEV, 93000, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010)
     &                 'I/O error', IOS, 'reading COUNTY_TEMP_RWC '//
     &                 'file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Skip blank and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Sparse line
                CALL PARSLINE( LINE, 2, SEGMENT )
                
                IFIP  = STR2INT ( SEGMENT( 1 ) )     ! FIP ID
                DTEMP = STR2REAL( SEGMENT( 2 ) )     ! min temp set for RWC eq

                DO S = 1, NSRGFIPS

                    ISTA = STR2INT( SRGFIPS( S )( 8:9 ) )
                    ISRGFIP = STR2INT( SRGFIPS( S ) )

                    IF( IFIP == ISTA .OR. IFIP == ISRGFIP ) THEN
                        RWC_TEMP( S ) = DTEMP
                    END IF

                END DO

            END DO    ! end of loop

             IF( N == 0 ) THEN
                 MESG = 'ERROR: No entries in COUNTY_TEMP_RWC file'
                 CALL M3MSG2( MESG )
             END IF
            
            CLOSE( RDEV )

        END IF

C.........  Get episode starting date and time and ending date
        MESG = 'Episode start date (YYYYDDD)'
        EPI_SDATE = ENVINT( 'STDATE', MESG, 0, IOS )

        MESG = 'Episode start time (HHMMSS)'
        EPI_STIME = ENVINT( 'STTIME', MESG, 0, IOS )

        MESG = 'Episode end date (YYYYDDD)'
        EPI_EDATE = ENVINT( 'ENDATE', MESG, 0, IOS )

        MESG = 'Episode end time (HHMMSS)'
        EPI_ETIME = ENVINT( 'ENDTIME', MESG, 230000, IOS )

C.........  Check start date is Jan 1st for a proper processing.
        TYEAR = INT( EPI_SDATE / 1000 )
        TDATE = TYEAR * 1000 + 1
        IF( EPI_SDATE /= TDATE .AND. .NOT. NH3FLAG ) THEN
            TDATE = TYEAR * 1000 + 1 
            WRITE( MESG,94010 ) 'ERROR: MUST set starting date (STDATE) to ',
     &           TDATE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        TDAY = 1 / YR2DAY( TYEAR )
        TDATE = TYEAR * 1000 + TDAY 
        IF( EPI_EDATE /= TDATE .AND. .NOT. NH3FLAG ) THEN
            TDAY = TYEAR * 1000 + TDAY
            WRITE( MESG,94010 ) 'ERROR: MUST set ending date (ENDATE) to ',
     &             TDATE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Output temporal x-ref output header for TREF_OUT file
        WRITE( XODEV,93000 ) '#TREF'
        WRITE( XODEV,94040 ) '#YEAR=', INT( EPI_SDATE/1000 )
        WRITE( XODEV,93000 ) '#GRDNAME=' // GDNAM3D // COORD
        WRITE( XODEV,93000 ) '#GENTPRO_PROFILE_METHOD = ' // TRIM( PROF_METHOD )
        WRITE( XODEV,94010 ) '#PERIOD=', EPI_SDATE, '-', EPI_EDATE
        WRITE( XODEV,93000 ) '#DESC:Temporal cross-reference file ' //
     &                       '(TREF) generated by GenTPRO program'
        WRITE( XODEV,93000 ) '#SCC,FIPS,FacilityID,UnitID,ReleasePoint,'//
     &                       'ProcessID,PollutantID,ProfileType,ProfileID,Comments...'

C.........  Read data from temporal x-ref input file.
C           Only retain lines that contain both SCC and FIPS information.
C.........  Get number of lines in x-ref input file
        MXLINE=GETFLINE( XIDEV, 'Temporal x-ref input file [TREF_IN]' )

        CALL M3MSG2( 'Reading Temporal x-ref input file' )

C.........  Allocate local arrays
        MXTREF = NSCC * ( NSRGFIPS + NSTA + 1 ) * 2
        ALLOCATE( CSCCFIP( MXTREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCFIP', PROGNAME )
        ALLOCATE( INDXREF( MXTREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXREF', PROGNAME )
        ALLOCATE( MATCHED( NSCC,NSRGFIPS+NSTA+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MATCHED', PROGNAME )
        CSCCFIP = ''
        INDXREF = 0
        MATCHED = 0

C.........  Count a list of x-ref entries in TREF
        NMATCH = 0
        DO I = 1, MXLINE

            READ( XIDEV, '(A)', IOSTAT=IOS ) LINE
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 )
     &             'I/O error', IOS, 'reading temporal cross-' //
     &             'reference file at line', I
                CALL M3MESG( MESG )
            END IF

C............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

            CALL PARSLINE( LINE, MXSEG, SEGMENT )

            CSCC = TRIM   ( SEGMENT( 1 ) )
            FIPS = STR2INT( SEGMENT( 2 ) )
            CFIPS = TRIM  ( SEGMENT( 2 ) )
            CALL PADZERO( CFIPS )
            CPOL = TRIM   ( SEGMENT( 7 ) )
            TPRO = TRIM   ( SEGMENT( 8 ) )

            IF( INDEX1( TPRO, NTPRO, TPROTYPES ) < 1 ) THEN
                MESG = 'ERROR: Older TREF format is no longer supported!'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF            

C.............  temporary limit for supporting older TREF format (v3.5.1 or earlier)
            IF( FIPS < 1 ) FIPS = -9
            IF( CPOL == '' .OR. CPOL == '0' ) CPOL = '-9'

            L0 = INDEX1( CSCC, NSCC, SCCLIST )
            L1 = FINDC ( CFIPS, NSRGFIPS, SRGFIPS )

            LL = 0
            NS = 0
C.............  find matched SCC entries
            IF( L0 > 0 ) THEN

                IF( L1 < 1 ) THEN
                    ISTA = STR2INT( CFIPS( 8:9 ) )
                    NS = FIND1( ISTA, NSTA, SRGSTA )
                    IF( NS   > 0 ) L1 = NSRGFIPS + NS         ! state-specific entry in XREF fle
                    IF( FIPS < 1 ) L1 = NSRGFIPS + NSTA + 1   ! SCC-specific ultimate default (No FIPS)
                END IF

                IF( MONAVER .AND. TPRO == 'MONTHLY' ) THEN
                    LL = 1
                ELSE IF( DAYAVER .AND. TPRO == 'DAILY' ) THEN
                    LL = 1
                ELSE IF( HOURAVER .AND. TPRO == 'HOURLY' ) THEN
                    LL = 1
                END IF
                
                IF( L1 < 0 ) LL = 0  ! L1 can be zero if fips can not be found from XREF input file

                WRITE( MESG, 94010 ) 'ERROR: pollutant-specific '//
     &                'entry for target SCC: ' // CSCC // ' at line', I, 
     &                ' is NOT supported'

                IF( NH3FLAG ) THEN
                    IF( CPOL /= '-9' .AND. CPOL /= 'NH3' )
     &                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    IF( CPOL /= '-9' )
     &                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

C.............  check pollutant-specific x-ref entries while processing AGNH3 mothod
            IF( LL > 0 ) THEN

C.................  Concatenate all segments into a string
                LINE = ''
                DO J = 1, 9
                    LINE = TRIM( LINE ) //  TRIM( SEGMENT( J ) ) // ','
                END DO

C.....................  Store all TREF entries including new ones
                NN = INDEX1( LINE, MXTREF, CSCCFIP )
                IF( NN > 0 ) CYCLE
                NMATCH = NMATCH + 1
                INDXREF( NMATCH ) = I
                MATCHED( L0,L1  ) = NMATCH
                CSCCFIP( NMATCH ) = TRIM( LINE )

            END IF

        END DO
        
        REWIND( XIDEV )

C.........  Concatenate org x-ref entries and add new FIPS if necessary
C.........  Store a list of entries that matches with processing FIPS/SCC
C.........  Find ultimate FIPS/SCC x-ref when no matched in TREF
        CALL M3MSG2('Processing Temporal x-ref input file..........')

        NWARN = 0
        DO I = 1, NSCC
            
            DO J = 1, NSRGFIPS
            
                LL = 0
                L0 = I
                L1 = J

                CSCC = SCCLIST( I )
                CFIPS = SRGFIPS( J )
                CALL PADZERO( CFIPS )

C.................  Check secondly SCC/FIPS specific entries
                LL = MATCHED( L0,L1 )
 
C.................  Check a list of state-specific entry first before using default profiles
                IF( LL < 1 ) THEN
                    ISTA = STR2INT( CFIPS( 8:9 ) )
                    NS = FIND1 ( ISTA, NSTA, SRGSTA )
                    IF( NS > 0  ) L1 = NSRGFIPS + NS         ! state-specific entry in XREF fle
                    LL = MATCHED( L0,L1 )
                END IF

C.................  Check a list of SCC-specific entry first before using default profiles
                IF( LL < 1 ) THEN
                    L1 = NSRGFIPS + NSTA + 1   ! SCC-specific ultimate default (No FIPS)
                    LL = MATCHED( L0,L1 )
                END IF

C.................  Use default profiles when there no matched x-ref entry
                LINE = ''
                IF( LL < 1 ) THEN

                    IF( NH3FLAG ) THEN
                        WRITE( LINE,93000 ) CSCC//','//CFIPS //',,,,,NH3,'
                    ELSE
                        WRITE( LINE,93000 ) CSCC//','//CFIPS //',,,,,-9,'
                    END IF

                    WRITE( MESG,93000 ) 'WARNING: Could not find ' //
     &                  'x-ref entries for SCC: ' // CSCC // 
     &                  ' and FIPS: '// CFIPS // CRLF() // BLANK10
     &                   //'New x-ref entries for these SCC and FIPS '//
     &                   'are added to TREF_OUT output file'
                    NWARN = NWARN + 1
                    IF( NWARN < MXWARN ) CALL M3MSG2( MESG )

C.................  Concatenate all segments into a string
                ELSE

                    SEGMENT = ''
                    CALL PARSLINE ( CSCCFIP( LL ), MXSEG, SEGMENT )

                    SEGMENT( 1 ) = CSCC
                    SEGMENT( 2 ) = CFIPS
                    IF( NH3FLAG ) SEGMENT( 7 ) = 'NH3'

                    DO L = 1, 7
                        LINE = TRIM( LINE )//TRIM( SEGMENT( L ) ) // ','
                    END DO

                END IF

                TPROID = CFIPS
                
C.................  Append new MONTHLY temporal profile ID to new/existing x-ref entry
                IF( MONAVER ) THEN
                    WRITE( XODEV,'( A )' ) TRIM(LINE)//'MONTHLY,'//TPROID//',""'
                END IF

C.................  Append new DAILY temporal profile ID to new/existing x-ref entry
                IF( DAYAVER ) THEN
                    WRITE( XODEV,'( A )' ) TRIM(LINE)//'DAILY,'//TPROID//',""'
                END IF
                    
C.................  Append new HOURLY temporal profile ID to new/existing x-ref entry
                IF( HOURAVER ) THEN
                    WRITE( XODEV,'( A )' ) TRIM(LINE)//'HOURLY,'//TPROID//',""' 
                END IF

            ENDDO

        ENDDO

C.........  Store the rest of unmatched original x-ref entries        
        DO I = 1, MXLINE

            LL = 0
            L2 = 1
            IF( NH3FLAG ) L2 = 2

            READ( XIDEV, '(A)', IOSTAT=IOS ) LINE

C............  Parse line into substrings (segments)
            IF( BLKORCMT( LINE ) ) CYCLE
               
            CALL PARSLINE ( LINE, MXSEG, SEGMENT )

C.............  Skip all matched org x-ref entries
            DO J = 1, MXTREF
                IF( I == INDXREF( J ) ) THEN
                    LL = J
                END IF
            END DO

            IF( LL > 0 ) CYCLE

C.............  Output all unmatched Xref entries
            WRITE( XODEV,93000 ) TRIM( LINE )

        END DO

C.........  Deallocate unnecessary local arrays
        DEALLOCATE( CSCCFIP )

C.........  Find the total number of time steps
        EPI_NSTEPS = 1 +
     &      SECSDIFF( EPI_SDATE,EPI_STIME,EPI_EDATE,EPI_ETIME ) / 3600
 
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

C.........  Calculate required ending date and time based on episode settings
        EDATE = EPI_EDATE + 1
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

            IF( NH3FLAG ) THEN
                J = INDEX1( VAR_NAME, NVARS3D, VNAME3D )
                IF( J <= 0 ) THEN
                    MESG = 'ERROR: Could not find ' //TRIM( VAR_NAME )//
     &                     ' in file ' // TRIM( METFILE )
                    CALL M3MESG( MESG )

                    IF( BASHFLAG ) THEN
                        MESG = 'NOTE: Both ground temperature and '
     &                      // 'aerodynamic resistance variables need '
     &                      // 'to be in the smae MET file'
                    ELSE
                        MESG = 'NOTE: Both temperature and wind speed ' //
     &                     'variables need to be in the same MET file ' //
     &                     CRLF() // BLANK10 //
     &                     'NOTE: User may need to use Metcombine '//
     &                     'SMOKE utility program to combine met '//
     &                     'variables from MET_CRO_2D and MET_CRO_3D'
                    END IF

                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF
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

C.............  Find starting position in METDAYS array
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

C.........  Setting up all output headers
C.........  Write HDRer to month-of-year temporal profile output file
        IF( MONAVER ) THEN
            WRITE( MODEV,93000 ) '#FORMAT=TPRO_MON'
            WRITE( MODEV,94040 ) '#NUM_PROFILES = ', NSRGFIPS
            WRITE( MODEV,94040 ) '#YEAR=', INT( EPI_SDATE/1000 )
            WRITE( MODEV,93000 ) '#GRDNAME=' // GDNAM3D // COORD
            WRITE( MODEV,94010 ) '#PERIOD=', EPI_SDATE ,'-', EPI_EDATE
            WRITE( MODEV,93000 ) '#DESC: Month of year profiles from '//
     &          'GenTPRO using profile method :: '//TRIM(PROF_METHOD)
            WRITE( MODEV,93000 ) '#PROFID,JAN,FEB,MAR,APR,MAY,,,,DEC'
        END IF

C.........  Write HDRer to day-of-month temporal profile output file
        IF( DAYAVER ) THEN
            WRITE( DODEV,93000 ) '#FORMAT=TPRO_DAY'
            WRITE( DODEV,94040 ) '#NUM_PROFILES = ', NSRGFIPS
            WRITE( DODEV,94040 ) '#YEAR=', INT( EPI_SDATE/1000 )
            WRITE( DODEV,93000 ) '#GRDNAME=' // GDNAM3D // COORD
            WRITE( DODEV,94010 ) '#PERIOD=', EPI_SDATE ,'-', EPI_EDATE
            WRITE( DODEV,93000 ) '#DESC: Day of month profiles from'//
     &          'GenTPRO using profile method :: '//TRIM(PROF_METHOD)
            WRITE( DODEV,93000 ) '#PROFID,MON,DAY1,DAY2,DAY3,...,DAY31'
        END IF

C.........  Allocate met variable arrays
        ALLOCATE( TA( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( VARTMP( METNGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VARTMP', PROGNAME )

C.........  Source met variable arrays
        ALLOCATE( TASRC( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
        ALLOCATE( VARSRC( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VARSRC', PROGNAME )
        ALLOCATE( MINTEMP( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MINTEMP', PROGNAME )
        TASRC = 0.0
        VARSRC = 0.0
        MINTEMP = -1 * BADVAL3

C.........  Allocate memory for storing hourly/annual meteorology profiles
        ALLOCATE( HRLSRC( NSRGFIPS,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HRLSRC', PROGNAME )
        ALLOCATE( TOTSRC( NSRGFIPS,NSTEPS,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOTSRC', PROGNAME )
        ALLOCATE( DAYSRC( NSRGFIPS,MXDAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYSRC', PROGNAME )
        ALLOCATE( MONSRC( NSRGFIPS,NMONTH ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MONSRC', PROGNAME )
        ALLOCATE( ANNSRC( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ANNSRC', PROGNAME )
        ALLOCATE( TMPDSRC( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPDSRC', PROGNAME )
        ALLOCATE( TMPMSRC( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPMSRC', PROGNAME )
        TOTSRC = 0.
        HRLSRC = 0.
        DAYSRC = 0.
        MONSRC = 0.
        ANNSRC = 0.
        TMPDSRC = 0.
        TMPMSRC = 0.

C.........  dates/daylight saving arrays
        ALLOCATE( DAYBEGT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        ALLOCATE( DAYENDT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        ALLOCATE( PRVBEGT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRVBEGT', PROGNAME )
        ALLOCATE( PRVENDT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRVENDT', PROGNAME )
        ALLOCATE( PROCDAYS( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROCDAYS', PROGNAME )
        ALLOCATE( DT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DT', PROGNAME )
        ALLOCATE( MT( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MT', PROGNAME )
        PRVBEGT  = 0
        PRVENDT  = 0
        DAYBEGT  = 0
        DAYENDT  = 0
        PROCDAYS = 0
        DT       = 0
        MT       = 0

C.........  Process meteorology data...
        MESG = 'Processing meteorology data using variables ' //
     &         TRIM( TVARNAME ) // ', ' // TRIM( VAR_NAME ) // '...'
        CALL M3MSG2( MESG )
        
        JDATE = SDATE
        JTIME = STIME
        LDATE = -9

C.....  Loop through days/hours (timesteps) of meteorology files
        DO T = 1, NSTEPS

C.........  When new day...
           IF( JDATE /= LDATE ) THEN

C.................  Write message for day of week and date
                DAY = WKDAY( JDATE )
                MESG = 'Processing ' // DAYS( DAY ) // ' ' //
     &                 MMDDYY( JDATE )
                CALL M3MSG2( MESG )

C.................  Set start and end hours of previous day for all sources
                CALL SETSRCDY( NSRGFIPS, LDATE, TZONES, USEDAYLT, .TRUE.,
     &                         PRVBEGT, PRVENDT )

C.................  Set start and end hours of current day for all sources
                CALL SETSRCDY( NSRGFIPS, JDATE, TZONES, USEDAYLT, .TRUE.,
     &                         DAYBEGT, DAYENDT )

                DO S = 1, NSRGFIPS
                    IF( USEDAYLT( S ) ) THEN
                        DST = DAYBEGT( S ) - PRVBEGT( S )   ! -1 roll-back one hour for a proper
                        EXIT
                    END IF
                END DO
            END IF

C.............  Determine input file for this hour
            POS = T

C.............  Get file number for current timestep
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
                    END IF   ! ( FILEOPEN )

                    PREVFILE = METFILE

                END IF   ! ( METFILE .NE. PREVFILE )

C.................  Set logical file name
                IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                    MESG = 'Could not set logical file name for ' //
     &                     'file ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Open the meteorology data file
                IF( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                    MESG = 'Could not open meteorology file ' //
     &                     TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                ELSE
                    FILEOPEN = .TRUE.
                END IF

C.................  Read current meteorology file
                IF( .NOT. READ3( METNAME, TVARNAME, 1,
     &                            JDATE, JTIME, TA ) ) THEN
                    MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                     ' from ' // TRIM( METFILE )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Average the temperatures of the grid cells that 
C                   intersect with the selected
C                   surrogate. Results are stored in TASRC.
C                   Equation (1) in design document.
                CALL GRDFIPS( NSRGFIPS, SRGFIPS, TA, TASRC, .FALSE. )

C.................  Read current 2nd variable for AGNH3 profile method processing
                IF( NH3FLAG ) THEN
                    IF( .NOT. READ3( METNAME, VAR_NAME, 1,
     &                                JDATE, JTIME, VARTMP ) ) THEN
                        MESG = 'Could not read ' // TRIM( VAR_NAME ) //
     &                     ' from ' // TRIM( METFILE )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

C.....................  Apply ungridding matrix to average the 2nd variable
C                   of the grid cells that intersect with the selected
C                   surrogate. Results are stored in VARSRC.
                    CALL GRDFIPS( NSRGFIPS, SRGFIPS, VARTMP, VARSRC, .FALSE. )

                END IF   ! ( NH3FLAG )

            ELSE
                MESG = 'ERROR: Missing meteorology file on '//
     &                  MMDDYY( JDATE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF  ! check for using alternate data or day averaging

C.............  If last day of month, process monthly averages
            CALL DAYMON( JDATE, TMPMNTH, TDAY )   ! next day of processing date
            CALL DAYMON( JDATE-1, MONTH, DAY )       ! processing month and date

C.............  Processing profile methods (RWC, AGNH3 or MET)
C.............  Loop over source (S) to estimate hourly values and temporal sums
            DO S = 1, NSRGFIPS

C.................  Convert GMT to local time
                HOURIDX = ( JTIME-DAYBEGT( S ) ) / 10000
                IF( HOURIDX < 0 ) HOURIDX = HOURIDX + 24
                HOURIDX = HOURIDX + 1

C.................  Skip if begining/ending hours are out of range
                IF( T <= TSPREAD .AND. HOURIDX > TSPREAD ) CYCLE
                IF( ( T > EPI_NSTEPS .AND. HOURIDX <= TSPREAD ) .OR.
     &                T > EPI_NSTEPS + TSPREAD ) THEN
                    NT = T - EPI_NSTEPS
                    HRLSRC( S,T )   = HRLSRC( S,NT )
                    TOTSRC( S,T,: ) = TOTSRC( S,NT,: )
                    CYCLE
                END IF

C.................  Skip if data is missing
                IF( TASRC( S ) == 0.0 ) THEN  ! temp in Kevin 
                     MESG = 'ERROR: Incorrect temperature '//
     &                  'value on'//MMDDYY( JDATE )//' of county '//
     &                  SRGFIPS( S )
                   CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Choose profile method calcuation
                IF( PROF_METHOD == 'RWC' ) THEN

                    TEMPVAL = 1.8 * TASRC( S ) - 459.67   ! K --> F

C.................  When new day...
                    IF( HOURIDX == 24 ) THEN
C.....................  RWC: Compute percentage of emissions (PE) from
C                       the UNC RWC regression equation (Adelman, 2010)
C                       to calculate temporal profiles for residential
C                       wood combustion sources:
C                          PE = (42.12 - 0.79*T) / Sum(42.12 - 0.79*T) (default equation)
C                       where T is the minimum ambient temperature
C                       in a county (not allowed to exceed 50 F).
C                       Numerator in Equation (2) of the design docment.
C                       Note: Even though this is labeled as an hourly
C                             array, it is only calculated once per day.
C                             Values for the other hours are 0.

C 1. original equation with cap T:
C    if Td <= 50: use equation 42.12 - 0.79*Td
C    if Td > Tt: 0
C    if 50 < Td < Tt: use equation 42.12 - 0.79*50

C    - adv:  preserves shape of profile below 50 degrees
C 	   still based on original regression equation

C    - disadv:  puts only a very small % of profile in Td > 50

C 2. new equation
C    if Td >= Tt: 0
C    if Td < Tt: use equation 0.79*(Tt -Td)

C    - adv:  works for any Tt w/o an additional if clause
C 	   places significant % of profile in Td > 50

C    - disadv: departs from the origial regression equation
C 	     maybe overly flattens the profile for Tt > 50

                        IF( ZFLAG ) THEN
                            HRLSRC( S,T ) = SLOPE * ( RWC_TEMP(S) - MINTEMP(S) )
                        ELSE
                            IF( MINTEMP(S) > 50.0 .AND. MINTEMP(S) <= RWC_TEMP(S) ) THEN
                                MINTEMP( S ) = 50.0
                            END IF
                            HRLSRC( S,T ) = CONST - ( SLOPE *  MINTEMP( S ) )
                        END IF
                        
                        IF( MINTEMP( S ) > RWC_TEMP( S ) ) THEN
                            HRLSRC( S,T ) = 0.0   ! set it to zero when mintemp > 50F
                        END IF

                        MINTEMP( S ) = -1 * BADVAL3   ! Reset MIN(temp) back to flag value

                    ELSE   ! Another hour, but not a new day yet.
C.....................  Continue finding the lowest hourly temperature
C                       in the day for the RWC regression equation.

                        HRLSRC( S,T ) = 0.
                        IF( TEMPVAL == 0. ) CYCLE
                        MINTEMP( S ) = MIN( MINTEMP( S ), TEMPVAL )

                    END IF   ! IF  (new day)

                ELSE IF( NH3FLAG ) THEN

                    IF( .NOT. BASHFLAG ) THEN        ! Russell and Cass Equation for Agricultural NH3

C.....................  AGNH3: Use Russell and Cass (1986) algorithm for
C                              predicting diurnal NH3 emission variations
C                              to calculate temporal profiles for
C                              agricultural ammonia sources using Temp(C) and WS(m/s).
                        TEMPVAL = TASRC( S ) - 273.15
                        HRLSRC( S,T ) = 2.36**(TEMPVAL/10.) * VARSRC( S )

                    ELSE           ! Bash Equation for Livestock NH3

                        TEMPVAL = 161500./TASRC(S) * EXP( -1*10380./TASRC(S) ) ! Henry's equilibrium  
                        HRLSRC( S,T ) = TEMPVAL * VARSRC( S )   ! H(t) / Ra (aerodynamic resistance)
                    
                    END IF


                ELSE
C.................  MET: Use the time series of the selected
C                        meteorological variable to calculate
C                        temporal profiles.
                    HRLSRC( S,T ) = TASRC( S )    ! kevin temp

                END IF

C.................  Calculate tmp daily total from hourly values
                TMPDSRC( S ) = TMPDSRC( S ) + HRLSRC( S,T )

C..................  add first hour=1 value to previous daily total as 
C                    DST adjustment for a proper Temporal processing. 
                 IF( HOURIDX == 1 .AND. DST == -10000 ) THEN
                     DAYSRC( S,NDAY ) = DAYSRC( S,NDAY ) + HRLSRC( S,T )
                     PT = DT( S ) - 23
                     TOTSRC( S,PT:T-1,1 ) = DAYSRC( S,NDAY )
                END IF

                IF( HOURIDX == 24 .AND. DST ==  10000 ) THEN
                    TMPDSRC( S ) = TMPDSRC( S ) - HRLSRC( S,T )  ! skip one hour value
                END IF

                IF( HOURIDX == 24 ) THEN

C.....................  County-specific procesing days
                    PROCDAYS( S ) = PROCDAYS( S ) + 1
                    NDAY = PROCDAYS( S )

C.....................  Sum hourly to daily total in local time
                    DAYSRC( S,NDAY ) = TMPDSRC( S )
                    
C.....................  Store daily total to TOTSRC array for ncf output file
                    IF( DT( S ) == 0 ) THEN
                        DT( S ) = T - 23
                        MT( S ) = T - 23
                    END IF

                    PT = DT( S )
                    TOTSRC( S,PT:T,1 ) = DAYSRC( S,NDAY )
                    DT( S ) = T + 1                             ! store starting time step for next day total

                    TMPDSRC( S ) = 0.0   ! reset tmp daily total array

C.....................  Sum daily to monthly total in local time
                    TMPMSRC( S ) = TMPMSRC( S ) + DAYSRC( S,NDAY )

                    IF( MONTH /= TMPMNTH ) THEN

C......,,,,...............  Sum daily to monthly total in local time
                        MONSRC( S,MONTH ) = TMPMSRC( S )

C.........................  Store daily total to TOTSRC array for ncf output file
                        PT = MT( S )
                        TOTSRC( S,PT:T,2 ) = MONSRC( S,MONTH )
                        MT( S ) = T + 1                           ! store starting time step for next day total

C.........................  Sum Monthly to Annual total in local time
                        ANNSRC( S ) = ANNSRC( S ) + MONSRC( S,MONTH )

                        TMPMSRC( S ) = 0.0    ! reset tmp monthly total array

                    END IF

                END IF

            END DO   ! source loop: S = 1, NSRGFIPS

C.............  Increment loop time
            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, 10000 )

        END DO   ! End T-loop on hours of met files

C.........  Reset annual total to 1.0 if it is zero
        DO S = 1, NSRGFIPS 
            IF( ANNSRC( S ) == 0.0 ) THEN
                ANNSRC( S ) = 1.0
                DO I = 1, NSCC
		     MESG = 'CRITICAL WARNING: All temporal'//
     &               ' profiles are ZERO for county: '//SRGFIPS( S )// 
     &               ' and SCC: ' // SCCLIST( I )
                    CALL M3MSG2( MESG ) 
                END DO
                MESG = ':: It will result in zeroing out inventory '//
     &                 'if you use these temporal profiles'
                CALL M3MSG2( MESG )
            END IF
        END DO

C.........  Deallocate arrays no-longer needed.
        DEALLOCATE( TMPDSRC, TMPMSRC )

C.........  Output Temporal profiles and new/updated x-ref entries
C.........  Compute daily/monthly temporal profile fractions and generate output

C.........  Allocate memory for temporal profile data
        ALLOCATE( PROF_MON( NSRGFIPS,NMONTH ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROF_MON', PROGNAME )
        ALLOCATE( PROF_DAY( NSRGFIPS,31 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROF_DAY', PROGNAME )
        PROF_MON = 0.
        PROF_DAY = 0.

C......... Output monthly temporal profiles
C.........  Compute month of year temporal profiles
        IF( MONAVER ) THEN

            DO MM = 1, NMONTH
                PROF_MON( :,MM ) = MONSRC( :,MM ) / ANNSRC( : )
            END DO

C.............  Output monthly profiles by county
            DO S = 1, NSRGFIPS
                WRITE( MODEV, "(A,12(',',E10.3))" ) SRGFIPS( S ),
     &              ((PROF_MON( S,NP )), NP = 1,12 )
            END DO

        END IF

C.........  Deallocate local arrays        
        DEALLOCATE( PROF_MON )

C......... Output daily temporal profiles
C.........  Compute day of month temporal profiles
        IF( DAYAVER ) THEN

            DO DD = 1, MAXVAL( PROCDAYS )

                JDATE = INT( SDATE/1000 ) * 1000 + DD
                CALL DAYMON( JDATE+1, TMPMNTH, TDAY )
                CALL DAYMON( JDATE  , MONTH  , DAY  )

                DO S = 1, NSRGFIPS
                    IF( MONSRC( S,MONTH ) /= 0.0 ) THEN
                        PROF_DAY( S,DAY ) = DAYSRC( S,DD ) / MONSRC( S,MONTH )
                    ELSE
                        PROF_DAY( S,DAY ) = 0.0 
                    END IF
                END DO

C.................  Output daily profiles by county
                IF( MONTH /= TMPMNTH ) THEN

                    DO S = 1, NSRGFIPS
                        WRITE( DODEV, "(A,A,I2.2,31(',',E10.3))" )
     &                      SRGFIPS( S ), ',', MONTH,
     &                      ( (PROF_DAY( S,NP ) ), NP = 1,31 )
                    END DO

                    PROF_DAY = 0.0    ! re-initializing 

                END IF

            END DO

        END IF

C.........  Deallocate local arrays        
        DEALLOCATE( PROF_DAY )

C.........  Output daily/monthly/annual total values for Temporal to compute
C           hour-of-day, hour-of-month, hour-of-year hourly profiles as an option

C.........  Open output file for hourly temporal profiles
        IF( HOURAVER ) THEN

C.............  Initialize I/O API output file headers
            CALL HDRMISS3

            FDESC3D( 1 ) = 'GENTPRO-based meteorology profiles file'
            FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
            FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
            WRITE( FDESC3D( 4 ), 93000 ) '/PROFILE_METHOD/ ' // PROF_METHOD
            WRITE( FDESC3D( 5 ), 93000 ) '/TPRO_OUTPUT/ '// ' HOURLY'
            WRITE( FDESC3D( 6 ), 93000 ) '/T_UNITS/ "deg K"'
            FDESC3D( 7 ) = '/T_VNAME/ ' // TVARNAME
            FDESC3D( 8 ) = '/NOTE/ Time 000000 in file represents ' //
     &                     '0 hour in output time zone'
            WRITE( FDESC3D( 9 ), 94010 ) '/END DATE/ ', EDATE

            FDESC3D( 21 ) = '/INVEN FROM/ ' // 'N/A'
            FDESC3D( 22 ) = '/INVEN VERSION/ ' // 'N/A'

C.............  Set header values that cannot be default
            JDATE = SDATE
            JTIME = 0000

            SDATE3D = JDATE
            STIME3D = JTIME
            TSTEP3D = 10000
            NROWS3D = NSRGFIPS
            NLAYS3D = 1
            NVARS3D = 5

            J = 1
            VNAME3D( J ) = 'COUNTIES'
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'County FIPS code'
            VTYPE3D( J ) = M3INT

            J = 2
            VNAME3D( J ) = 'HRLSRC'
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Hourly total by source'
            VTYPE3D( J ) = M3REAL

            J = 3
            VNAME3D( J ) = 'DAYTOT'
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Daily total by source'
            VTYPE3D( J ) = M3REAL

            J = 4
            VNAME3D( J ) = 'MONTOT'
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Monthly total by source'
            VTYPE3D( J ) = M3REAL

            J = 5
            VNAME3D( J ) = 'ANNTOT'
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 'Annual total by source'
            VTYPE3D( J ) = M3REAL

C.............  Open new file
            IF( .NOT. OPEN3( HNAME, FSUNKN3, PROGNAME ) ) THEN
                 MESG = 'Could not create new output file ' //
     &                 TRIM( HNAME )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            JDATE = SDATE
            JTIME = STIME

            DO T = 1, NSTEPS

C.................  Write county codes to file
                IF( .NOT. WRITE3( HNAME, 'COUNTIES', JDATE, JTIME,
     &                    SRGFIPS ) ) THEN
                     MESG = 'Could not write county codes to "' //
     &                      TRIM( HNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Write daily total by source to file
                IF( .NOT. WRITE3( HNAME, 'HRLSRC', JDATE, JTIME,
     &                   HRLSRC( :,T ) ) ) THEN
                     MESG = 'Could not write hourly total by sources'//
     &                      ' data to "' // TRIM( HNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Write daily total by source to file
                IF( .NOT. WRITE3( HNAME, 'DAYTOT', JDATE, JTIME,
     &                   TOTSRC( :,T,1 ) ) ) THEN
                     MESG = 'Could not write daily total by sources'//
     &                      ' data to "' // TRIM( HNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Write monthly total by source to file
                IF( .NOT. WRITE3( HNAME, 'MONTOT', JDATE, JTIME,
     &                   TOTSRC( :,T,2 ) ) ) THEN
                     MESG = 'Could not write daily total by sources'//
     &                      ' data to "' // TRIM( HNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Write monthly total by source to file
                IF( .NOT. WRITE3( HNAME, 'ANNTOT', JDATE, JTIME,
     &                   ANNSRC ) ) THEN
                     MESG = 'Could not write daily total by sources'//
     &                      ' data to "' // TRIM( HNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Increment loop time
                CALL NEXTIME( JDATE, JTIME, 10000 )

            END DO    ! loop over time

        END IF


C......... End program successfully

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I8, 1X, F13.5, 1X, F13.5, 1X, I8, 1X, A )

C...........   Internal buffering fosrmats............ 94xxx
94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, :, 2( I7.7, :, A ) )

94030   FORMAT( A, I5 )

94040   FORMAT( A, I4 )

94070   FORMAT( A, F5.1, A )


        END PROGRAM GENTPRO
