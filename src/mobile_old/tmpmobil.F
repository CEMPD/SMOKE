
C Version "@(#)$Id$ $Source$ $Date$ 


        PROGRAM TMPMOBIL

C***********************************************************************
C  program body starts at line 265
C
C  DESCRIPTION:
C       Produces source-level temporally allocated VMT using the 
C       mobile-source inventory vector produced by RAWMOBIL
C       and the temporal profiles and cross-reference tables for mobile
C       sources (extracted as the mobile-only portions of the AIRS tables
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level mobile source emissions in canonical order
C       produced by RAWMOBIL.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND1, FIND2, FIND3 GETEFILE, GETNUM, GETYN, DAYMON, 
C       TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C***********************************************************************

      IMPLICIT NONE


C...........   INCLUDES:

        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         CHKXLIN
        INTEGER         CSVNAME
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETEFILE
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDMPLIST
        INTEGER         RDTPROF
        INTEGER         STR2INT
        INTEGER         SEC2TIME
        INTEGER         TIME2SEC
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL    CHKXLIN, CSVNAME,
     &              FIND1, FIND2, FIND3, GETEFILE, GETNUM, GETYN, 
     &              INDEX1, ISDSTIME, PROMPTFFILE, PROMPTMFILE,
     &              RDMPLIST, RDTPROF, STR2INT, SEC2TIME, TIME2SEC, 
     &              TRIMLEN, WKDAY


C...........   PARAMETERS and their descriptions

        CHARACTER*16    BLANK16
        REAL*8          FTOC
        PARAMETER     ( BLANK16 = ' ' ,
     &                  FTOC    = 5.0D0 / 9.0D0 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobile Sources input and output arrays

        REAL        TOUT  ( NMSRC, NVTYPE, MXEMIS ) ! timestepped output emssions
        REAL        VMT   ( NMSRC )        ! vehicle miles traveled

        INTEGER     IFIP  ( NMSRC )        ! FIP codes for sources
        INTEGER     ILINK ( NMSRC )        ! link number or 0 for non-link
        INTEGER     IRCLAS( NMSRC )        ! EPA road classification code
        INTEGER     TFLAG ( NMSRC )        ! governs profile type usage
        INTEGER     TZONES( NMSRC )        ! time zones 

        INTEGER*2   DDEX  ( NMSRC, NVTYPE, MXEMIS ) ! wkday diurnal prof subscrpt
        INTEGER*2   EDEX  ( NMSRC, NVTYPE, MXEMIS ) ! wkend-diurnal prof subscrpt
        INTEGER*2   MDEX  ( NMSRC, NVTYPE, MXEMIS ) ! monthly       prof subscrpt
        INTEGER*2   WDEX  ( NMSRC, NVTYPE, MXEMIS ) ! weeky         prof subscrpt

        INTEGER     DREF  ( NVTYPE, MXEMIS )        ! diurnal temporary prof sub
        INTEGER     MREF  ( NVTYPE, MXEMIS )        ! monthly temporary prof sub
        INTEGER     WREF  ( NVTYPE, MXEMIS )        ! weekly  temporary prof sub

C...........   Vehicle Mix transform matrix

        REAL        MXMAT( NMSRC, NVTYPE )

C...........   MPLIST emissions factor cross reference tables:

        INTEGER   NMPLIST
        INTEGER   MPFIP( MXPLST )
        INTEGER   MPRDT( MXPLST )
        INTEGER   MPLNK( MXPLST )
        INTEGER   MPREF( 0:23, MXPLST )   !  hour&source-specific xrefs
        INTEGER   REFCNT                  !  count of unique list of PSIs
        INTEGER   REFLIST( MXPLST )       !  unique list of PSIs

        INTEGER   EFDEX( 0:23, NMSRC )    ! subscript into emis fac tables

C...........   Emissions factor tables:

        REAL    TEMP0		!  starting em factor table temperature
        REAL    TEMP1		!  ending em factor table temperature
        REAL    DTEMP		!  em factor table temperature increment (K)
        REAL    DDTEMP		!  inverse non-diurnal temperature increment (K)
        REAL    DNLTI           !  inverse diurnal temperature increment (K)

        INTEGER NTEMP		!  number of active elements in temp grid

        REAL	EXH_CO ( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	EXH_NOX( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	EXH_VOC( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	EVP_VOC( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	RNL_VOC( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	RST_VOC( MAXTMP , NVTYPE, MXPSI )  !  non-drl emissions factor
        REAL	WDL_VOC( NMINMAX, NVTYPE, MXPSI )  !  diurnal emissions factor
        REAL	DNL_VOC( NMINMAX, NVTYPE, MXPSI )  !  diurnal emissions factor
        REAL	HOT_VOC( NMINMAX, NVTYPE, MXPSI )  !  diurnal emissions factor
        REAL	CRC_VOC( NMINMAX, NVTYPE, MXPSI )  !  diurnal emissions factor

C.......   Temporal profiles:

        INTEGER         NMON, NWEK, NDIU, NEND	!  numbers of profiles

        INTEGER         MONREF( NMCOD )     ! Monthly         XREF table
        INTEGER         WEKREF( NWCOD )     ! Weekly          XREF table
        INTEGER         DIUREF( NDCOD )     ! Weekday-diurnal XREF table
        INTEGER         ENDREF( NDCOD )     ! Weekend-diurnal XREF table

        REAL            MONFAC( 12, NMCOD ) ! Monthly         profile coeffs
        REAL            WEKFAC(  7, NWCOD ) ! Weekly          profile coeffs
        REAL            DIUFAC( 24, NDCOD ) ! Weekday-diurnal profile coeffs
        REAL            ENDFAC( 24, NDCOD ) ! Weekend-diurnal profile coeffs

C...............   Temporal x-ref arrays

        INTEGER NPRFA               ! Actual number of XREF entries
        INTEGER INDXA( MXREF )      ! Index for sorted XREF list
        INTEGER FIPFA( MXREF )      ! FIPS in valid XREF list
        INTEGER RDTFA( MXREF )      ! Road class code in XREF list
        INTEGER LNKFA( MXREF )      ! Link number (user defined) in XREF list
        INTEGER VHTFA( MXREF )      ! Vehicle class in XREF list
        INTEGER CODEA( MXREF )      ! Emission type code in XREF list 
        INTEGER MPRFA( MXREF )      ! Monthly temporal profile for XREF
        INTEGER WPRFA( MXREF )      ! Weekly temporal profile for XREF
        INTEGER DPRFA( MXREF )      ! Daily temporal profile for XREF

        COMMON / TXREF / INDXA, FIPFA, RDTFA, LNKFA, VHTFA, 
     &                   CODEA, MPRFA, WPRFA, DPRFA

C...............  After sorting...
C.......  Default no-FIP, no-RDT
        INTEGER MPR00( NVTYPE, MXEMIS )
        INTEGER WPR00( NVTYPE, MXEMIS )
        INTEGER DPR00( NVTYPE, MXEMIS )

C.......  no-FIP, just RDT
        INTEGER NPRNF                           ! Number of actual entries

        INTEGER MPRNF( NRCLAS, NVTYPE, MXEMIS )
        INTEGER WPRNF( NRCLAS, NVTYPE, MXEMIS )
        INTEGER DPRNF( NRCLAS, NVTYPE, MXEMIS )
        INTEGER RDTNF( NRCLAS )                 ! RDT assoc w/each profile set

C.......  State, RDT
        INTEGER NPRNC                           ! Number of actual entries
        INTEGER MPRNC( MXSTRC, NVTYPE, MXEMIS )
        INTEGER WPRNC( MXSTRC, NVTYPE, MXEMIS )
        INTEGER DPRNC( MXSTRC, NVTYPE, MXEMIS )
        INTEGER FIPNC( MXSTRC )                 ! State assoc w/each profile set
        INTEGER RDTNC( MXSTRC )                 ! RDT assoc w/each profile set

C.......  FIP, RDT
        INTEGER NPRFC                           ! Number of actual entries
        INTEGER MPRFC( MXFRC, NVTYPE, MXEMIS )
        INTEGER WPRFC( MXFRC, NVTYPE, MXEMIS )
        INTEGER DPRFC( MXFRC, NVTYPE, MXEMIS )
        INTEGER FIPFC( MXFRC )                  ! FIP assoc w/each profile set
        INTEGER RDTFC( MXFRC )                  ! RDT assoc w/each profile set

C.......  FIP, LNK, (RDT)
        INTEGER NPRFL                           ! Number of actual entries
        INTEGER MPRFL( MXFLK, NVTYPE, MXEMIS )
        INTEGER WPRFL( MXFLK, NVTYPE, MXEMIS )
        INTEGER DPRFL( MXFLK, NVTYPE, MXEMIS )
        INTEGER FIPFL( MXFLK )                  ! FIP assoc w/each profile set
        INTEGER RDTFL( MXFLK )                  ! RDT assoc w/each profile set
        INTEGER LNKFL( MXFLK )                  ! LNK assoc w/each profile set

c...........   Time-stepping Matrix.  
C.......   TMAT( S,V,E,H ) is profile coeffs for 
C.......                    source S, pollutant V, emis type E, hour H

        REAL            TMAT( NMSRC, NVTYPE, MXEMIS, 24 )    ! coeff profiles

C...........   Use equivalence statement to reduce the required memory

        EQUIVALENCE ( TMAT, INDXA )

C...........   Ungridding matrix:

        INTEGER     NU( NMSRC )
        INTEGER     IU( NMATX )
        REAL        CU( NMATX )

        COMMON  / UGRIDMAT / NU, IU, CU

C...........   Ungridded minimum/maximum temperature arrays
        INTEGER     TMMI ( NMSRC )      ! Min/max indice
        REAL        TKMIN( NMSRC )      ! daily minimum temperature [deg K]
        REAL        TKMAX( NMSRC )      ! daily maximum temperature [deg K]

        COMMON  / MINMAX / TMMI, TKMIN, TKMAX

C...........   Other local variables

        LOGICAL  DAYLIT  !  true when TZONES have been changed to daylight time
        DATA     DAYLIT / .FALSE. /

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         OFLAG   !  sources outside grid
        LOGICAL         MFLAG   !  perform monthly adjustments flag
        LOGICAL         WFLAG   !  perform weekly adjustments flag

        INTEGER         COD     !  Temporary pollutant code
        INTEGER         EMT     !  Temporary MV emissions type
        INTEGER         ICOD    !  Index to pollutant code/emission type in list
        INTEGER         IDIU    !  Temporary diurnal profile number
        INTEGER         IMON    !  Temporary monthly profile number
        INTEGER         IOS     !  input status placeholder
        INTEGER         IWEK    !  Temporary weekly profile number
        INTEGER         IZONE   !  Zone of time stamps in temperature file
        INTEGER         FIP     !  Temporary state/county code
        INTEGER         FDATE
        INTEGER         JDATE	!  for output file
        INTEGER         JTIME	!  for output file
        INTEGER         JMDATE	!  Met date using output zone
        INTEGER         JMTIME	!  Met time using output zone
        INTEGER         JWDATE	!  Min/max temp date using output zone
        INTEGER         JWTIME	!  Min/max temp time using output zone
        INTEGER         JYEAR
        INTEGER         LNK     !  Temporary user defined link code
        INTEGER         MAXTL   !  Intermediate temperature indice variable
        INTEGER         MDATE	!  for met input file
        INTEGER         MTIME	!  for met input file
        INTEGER         NEMODE  !  Number of scenario-specific emission modes
        INTEGER         NETYPE  !  Number of scenario-specific emission types
        INTEGER         NIPOL   !  Number of scenario-specific inventory spcs
        INTEGER         NRDT    !  Index to road class in list
        INTEGER         NSTEPS
        INTEGER         OSRC    !  Number of sources outside grid
        INTEGER         RDT     !  Temporary EPA road class number
        INTEGER         TSTEP
        INTEGER         TZONE
        INTEGER         VHT     !  Temporary vehicle type number
        INTEGER         VID     !  Index for volatile inventory pollutant
        INTEGER         WDATE	!  Date of min/max temperature file
        INTEGER         WTIME	!  Time of min/max temperature file

        INTEGER         LDEV
        INTEGER         FDEV    ! unit number for emissions factors xref file
        INTEGER         VDEV    ! unit number for VMT mix file
        INTEGER         XDEV    ! ...for cross-reference input file
        INTEGER         PDEV    ! ...temporal profile input file
        INTEGER         TDEV    ! ...speciation list file

        INTEGER         H, HCORR, I, II, IREC  
        INTEGER         J, K, KU, L, M, N1, N2, F, S, T, V
        INTEGER         J1, J2, J3, J4   
        INTEGER         MONTH, DAY, LDATE
        REAL            TMMIDV

        REAL            SCR, PP, QQ
        REAL            R1, R2, S1, S2     
        REAL            TMIN0, TMAX0

        REAL            TA( NGRID )	!  one-layer slice of temperature
        REAL            TASRC           !  temperature, interpolated to srcs

        CHARACTER*3     EINAM ( MXIPOL ) ! Scenario-specific inventory pollutants
        CHARACTER*3     MEMODE( MXEPRC ) ! Scenario-specific MV processes
        CHARACTER*7     EMTNAM( MXEMIS ) ! Scenario-specific emission type names
  
        CHARACTER*3     INV     !  temporary string for inventory species name
        CHARACTER*3     PROC    !  temporary string for process name
        CHARACTER*16    ANAME   !  logical name for mobile-source     input file
        CHARACTER*16    DNAME   !  logical name for diurnal EFs       input file
        CHARACTER*16    FNAME   !  logical name for non-diurnal EFs input file
        CHARACTER*16    MNAME   !  logical name for met_cro_3d        input file
        CHARACTER*16    TNAME   !  logical name for timestepped      output file
        CHARACTER*16    UNAME   !  logical name for ungridding-matrix input file
        CHARACTER*16    WNAME   !  logical name for ungridded min/max tmpr input
        CHARACTER*256   LINE    !  buffer for reading XREF file
        CHARACTER*256   MESG    !  buffer for M3EXIT()


C***********************************************************************
C   begin body of program TMPMOBIL

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program TMPMOBIL to take a sorted source level mobile source',
     &  'emissions file, the ASCII temporal profile file, the ',
     &  'ASCII temporal cross-reference file, a met data file',
     &  'with temperature, an ungridding matrix, an ASCII vehicle mix',
     &  'file, a master emission factors file, and an emission factors',
     &  'cross-reference file, and produce temporally allocated',
     &  '(time-stepped) VMT for each vehicle class for the user-',
     &  'requested time period.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>"). The XREF file must',
     &  'have been sorted as indicated, prior to program execution.',
     &  'You will be prompted for which species to allocate -- these',
     &  'should only be the chemical species, and not the FIP-related',
     &  'or ASC-related variables.)', 
     &  ' '
        WRITE( *,92000 ) 
     &  'You will also be prompted for the TIME PERIOD to be covered',
     &  'by this run of the program, as specified by STARTING DATE',
     &  'in format YYYYDDD (= 1000*year + Julian day n0 (1...365,6) )',
     &  'STARTING TIME in format HHMMSS (= 10000*hr + 100*min + sec) )',
     &  'and DURATION in format HHMMSS.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 'Ending program', 2 )
        END IF

C.......   Get file name; open input mobile sources, temporal cross-reference,
C.......   and temporal profiles files

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for VMT INVENTORY input file',
     &          FSREAD3, 'MOBL', 'TMPMOBIL' )

        FNAME = PROMPTMFILE( 
     &          'Enter logical name for NON-DIURNAL EMISSIONS ' //
     &          'FACTOR input file', FSREAD3, 'MEFSND', 'TMPMOBIL' )

        DNAME = PROMPTMFILE( 
     &          'Enter logical name for DIURNAL EMISSIONS ' //
     &          'FACTOR input file', FSREAD3, 'MEFSD', 'TMPMOBIL' )

C.........  Prompt for time zone of temperature data.
        IZONE  = 0
        IZONE  = GETNUM( -12, 12, IZONE,
     &             'Enter time zone of temperature files (0 for GMT)' )

        MNAME = PROMPTMFILE( 
     &          'Enter logical name for SURFACE TEMPERATURE input file',
     &          FSREAD3, 'TKDAT01', 'TMPMOBIL' )

        UNAME = PROMPTMFILE( 
     &          'Enter logical name for UNGRIDDING MATRIX input file',
     &          FSREAD3, 'MUGMAT', 'TMPMOBIL' )

        WNAME = PROMPTMFILE( 
     &          'Enter logical name for UNGRIDDED MIN/MAX TEMPERATURE '
     &          // 'input file', FSREAD3, 'MINMAXT', 'TMPMOBIL' )

        VDEV = PROMPTFFILE( 
     &           'Enter logical name for VMT MIX file',
     &           .TRUE., .TRUE., 'MVMTM', 'TMPMOBIL' )

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS XREF file',
     &           .TRUE., .TRUE., 'MPLIST', 'TMPMOBIL' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., 'MTREF', 'TMPMOBIL' )

        PDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., 'MTPRO', 'TMPMOBIL' )

        TDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION LIST file',
     &           .TRUE., .TRUE., 'MSLST', 'TMPMOBIL' )


C.......   Get time period for output file:

        IF ( ISDSTIME( SDATE3D ) ) THEN
            TZONE = 4
        ELSE 
            TZONE = 5
        END IF
        TZONE  = GETNUM( -12, 12, TZONE, 
     &           'Enter time zone for output emissions (0 for GMT)' )
        JDATE  = GETNUM( 0, 9999999, 1988200, 
     &                   'Enter simulation starting date (YYYYDDD)' )
        JTIME  = GETNUM( 0, 999999, 0, 
     &                   'Enter simulation starting time  (HHMMSS)' )
        NSTEPS = GETNUM( 1, 999999, 24,
     &                   'Enter simulation duration (hours)' )

        TSTEP  = 10000
        JYEAR  = JDATE / 1000

C.........  Get starting date and time from met data file
        IF ( .NOT. DESC3( MNAME ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // MNAME( 1:TRIMLEN( MNAME ) ) // '"',
     &                  2 )
        END IF

C.........  Set   date/time for temperature based on simulation
C........   start date/time
        MDATE = JDATE
        MTIME = JTIME
        CALL NEXTIME( MDATE, MTIME, (TZONE - IZONE )*10000 )

C.........  Calculate start meterology date/time using TZONE
        JMDATE = SDATE3D
        JMTIME = STIME3D
        CALL NEXTIME( JMDATE, JMTIME, (IZONE - TZONE )*10000 )

C.........  Abort if (start time) < (first time for temperatures)
        IF( ( MDATE .EQ. SDATE3D .AND. MTIME .LT. STIME3D ) .OR.
     &        MDATE .LT. SDATE3D ) THEN
 
            WRITE( *,92000 ) ' '
            WRITE( *,92000 ) 'For file ' // MNAME // ' ...'
            WRITE( *,92000 ) ' '
            WRITE( *,92010 ) 'Given time zone       :', IZONE
            WRITE( *,92010 ) 'Requested date & time :', MDATE  , MTIME
            WRITE( *,92010 ) 'Actual    date & time :', SDATE3D, STIME3D

            CALL M3EXIT( 'TMPMOBIL', 0, 0,
     &                   'Requested starting date & time too early', 2 )
 
C.........  Give warning if modeling start date/time is not the same as 
C.........  the meterology start date/time

        ELSEIF( JDATE .NE. JMDATE .OR. JTIME .NE. JMTIME ) THEN

            WRITE( *,92000 ) ' '
            WRITE( *,92000 ) 'For file ' // MNAME // ' ...'
            WRITE( *,92000 ) ' '
            WRITE( *,92010 ) 'Given time zone       :', TZONE
            WRITE( *,92010 ) 'Requested date & time :', JDATE , JTIME
            WRITE( *,92010 ) 'Actual    date & time :', JMDATE, JMTIME

            CALL M3WARN( 'TMPMOBIL', 0, 0, 
     &                   'Inconsistent starting date and/or time' )

        ENDIF

C.........  Get starting date and time min/max temperature file (in GMT)
        IF ( .NOT. DESC3( WNAME ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // WNAME( 1:TRIMLEN( WNAME ) ) // '"',
     &                  2 )
        END IF

        WDATE = SDATE3D   ! From latest DESC3
        WTIME = STIME3D   ! From latest DESC3

C.........  Calculate start min/max temperature date/time using TZONE
C.........  the min/max temperature file is guaranteed to be in GMT
        JWDATE = WDATE
        JWTIME = WTIME
        CALL NEXTIME( JWDATE, JWTIME, (0 - TZONE)*10000 )

C.........  Abort if simulation date is less than date in MINMAXT file
        IF( JDATE .LT. JWDATE ) THEN   

            WRITE( *,92000 ) ' '
            WRITE( *,92000 ) 'For file ' // WNAME // ' ...'
            WRITE( *,92000 ) ' '
            WRITE( *,92010 ) 'Given time zone :', TZONE
            WRITE( *,92010 ) 'Requested date  :', JDATE
            WRITE( *,92010 ) 'Actual    date  :', JWDATE

            CALL M3EXIT( 'TMPMOBIL', 0, 0,
     &                   'Requested starting date too early', 2 )

        ELSEIF( WTIME .NE. 0 ) THEN

            CALL M3EXIT( 'TMPMOBIL', 0, 0, 'Error: time in ' // 
     &                   WNAME( 1:TRIMLEN( WNAME ) ) // 
     &                   ' is not 00000!', 2 )

        ELSEIF( JDATE .NE. JWDATE .OR. JTIME .NE. JWTIME ) THEN

            WRITE( *,92000 ) ' '
            WRITE( *,92000 ) 'For file ' // WNAME // ' ...'
            WRITE( *,92000 ) ' '
            WRITE( *,92010 ) 'Given time zone       :', TZONE
            WRITE( *,92010 ) 'Requested date & time :', JDATE , JTIME
            WRITE( *,92010 ) 'Actual    date & time :', JWDATE, JWTIME

            CALL M3WARN( 'TMPMOBIL', 0, 0, 
     &                   'Inconsistent starting date and/or time' )

        ENDIF

C.........  Update first min/max read date based on JDATE in GMT
C.........  Don't update time because time should always be zero
        WDATE = JWDATE

C.........  Get list of scenario-specific emission types in compiled order

        NETYPE  = CSVNAME( TDEV, .TRUE., .FALSE., .TRUE., EMTNAM )

        IF( NETYPE .LE. 0 ) THEN
            WRITE( MESG,94010 )
     &             'Error reading SPECIATION LIST file'
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        ENDIF

C.........  Get description of MRAW file for initialization of MTMP
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ANAME( 1:TRIMLEN( ANAME ) ) // '"',
     &                  2 )
        END IF

C.........  Get file name; open output time-stepped area sources file
        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = TSTEP
        NROWS3D = NMSRC
        NCOLS3D = 1
        NTHIK3D = IMISS3
        NLAYS3D = NVTYPE
        VGTOP3D = BADVAL3

C.........  Set/process the MTMP variable names
        NVARS3D = NETYPE
        DO 16 I = 1, NVARS3D

            VNAME3D( I ) = EMTNAM( I )

            L = TRIMLEN( VNAME3D( I ) )
            J = INDEX  ( VNAME3D( I )( 1:L ), '_' )

            PROC = VNAME3D( I )(   1 : J-1 )
            INV  = VNAME3D( I )( J+1 : L   ) 

C.............  Set other variable properties
            UNITS3D( I ) = 'gm/hour'
            VDESC3D( I ) = 'Hourly emissions for inventory species ' //
     &                     INV ( 1:TRIMLEN( INV ) ) // ' from process '
     &                     // PROC ( 1:TRIMLEN( PROC ) )
            VTYPE3D( I ) = M3REAL

16      CONTINUE

        FDESC3D( 1 ) = 
     &  'Temporally allocated source level mobile source VMT values.'
        DO  20 V = 1, NVTYPE
            WRITE( FDESC3D( V+1 ), 94010 )
     &      'Layer', V, 'represents vehicle type ' // VTYPE3( V )
20      CONTINUE
        DO  21 V = NVTYPE+2, MXDESC3
            FDESC3D( V ) = ' '
21      CONTINUE

        TNAME = PROMPTMFILE( 
     &  'Enter logical name for output TIME-STEPPED MOBILE EMIS file',
     &          FSUNKN3, 'MTMP', 'TMPMOBIL' )

C...........   Read in emissions indexing data IFIP, IRCLAS, ILINK.

        WRITE( *,92000 ) ' ', 'Reading VMT INVENTORY file.'

        IF ( .NOT. READ3( ANAME, 'IFIP', ALLAYS3, 0,0, IFIP ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading IFIP from file' // ANAME, 2)
        END IF

        IF ( .NOT. READ3( ANAME, 'IRCLAS', ALLAYS3, 0,0, IRCLAS ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading IRCLAS from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ILINK', ALLAYS3, 0,0, ILINK ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading ILINK from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'TZONES', ALLAYS3, 0,0, TZONES ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading TZONES from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'TPFLAG', ALLAYS3,0,0, TFLAG ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading TFLAG from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'VMT', ALLAYS3, 0, 0, VMT ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading "VMT" from file "' // 
     &                   ANAME( 1:TRIMLEN( ANAME ) ) // '".', 2 )
        END IF

C.........  Read vehicle mix fractions tables, and multiply with
C.........  source VMT totals
        CALL RDVMIX( VDEV, JYEAR, IFIP, IRCLAS, ILINK, VMT, MXMAT )

C.........  Read temporal-profiles file:  
C.........  4 parts (monthly, weekly, weekday diurnal, weekend diurnal)

        WRITE( *,92000 ) ' ', 'Reading TEMPORAL PROFILES file...', ' '

        NMON = RDTPROF( PDEV, 'MONTHLY', 12, NMCOD, MONREF, MONFAC )
        NWEK = RDTPROF( PDEV, 'WEEKLY',   7, NWCOD, WEKREF, WEKFAC )
        NDIU = RDTPROF( PDEV, 'WEEKDAY', 24, NDCOD, DIUREF, DIUFAC )
        NEND = RDTPROF( PDEV, 'WEEKEND', 24, NDCOD, ENDREF, ENDFAC )

C.........  Get list of scenario's emissions modes
        NEMODE = CSVNAME( TDEV, .FALSE., .FALSE., .TRUE. , MEMODE )

C.........  Get list of scenario's inventory pollutants
        NIPOL  = CSVNAME( TDEV, .TRUE. , .FALSE., .FALSE., EINAM  )

C.........   Decide which inventory pollutnat is volatile
        VID  = 0
        I    = 0
   33   CONTINUE  ! Loop to find volatile inventory pollutant

            I = I + 1
            VID = INDEX1( EINAM( I ), NUMVOL, VOLTYP )

            IF( VID .GT. 0 ) THEN 
                GO TO 35

            ELSEIF( I .EQ. NIPOL ) THEN

                I = I + 1
                CALL M3WARN( 'SPCMMAT', 0, 0,
     &                       'No valid volatile inventory pollutant' )
                GO TO 35

            ENDIF

        GO TO 33  ! To head of finding volatile pollutant type

   35   CONTINUE


C.........  Read and process temporal XREF 

        WRITE( *,92000 ) ' ', 'Reading TEMPORAL XREF file...', ' '

        IREC  = 0
        K     = 0
        EFLAG = .FALSE.
 
44      CONTINUE
 
            IREC = IREC + 1

            READ( XDEV, 93000, END=50, IOSTAT=IOS ) LINE

            IF( CHKXLIN( LINE, IOS, IREC, 'MOBILE' ) .NE. 0 ) THEN
                GO TO 44
            END IF

C.............  Ensure valid pollutant code
            COD = STR2INT( LINE( 32:36 ) )
            IF( COD .LE. 0 ) THEN
                ICOD = 0    ! Profile for all relevant emission types
            ELSE
                ICOD = FIND1( COD, MXIPOL, PCODE ) 
                IF ( ICOD .LT. 0 ) THEN    ! Skip if pollutant not on list
                    GO TO  44

                END IF
            END IF

C.............   Perform mobile-specific checks on XREF entry
            LNK  = MAX ( STR2INT( LINE( 13:15 ) ), 0 ) ! Use stack portion
            VHT  = STR2INT( LINE( 23:24 ) )
            RDT  = STR2INT( LINE( 25:27 ) )
            PROC =          LINE( 28:30 )
C.............  Inspect emission mode index

C.............  If process is invalid (numeric, like 000)
            IF( PROC .LT. 'A' ) THEN   ! Error
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Invalid emission process "' // LINE( 28:30 ) //
     &              '" at line', IREC, 'of TEMPORAL XREF file'
                CALL M3MESG( MESG )
                GO TO 44   !  to head of the XDEV-read loop

C.............  Set EMT and counter so will set specific emission type
            ELSEIF( ICOD .NE. 0 ) THEN

                EMT = INDEX1 ( PROC, NEMODE, MEMODE )

C.................  Skip if process is not used
                IF( EMT .LE. 0 ) THEN
                    GO TO 44

                ELSE
                    N1 = INDEX1( INVSPEC( ICOD ), MXIPOL, EINAM )
                    N2 = N1
                ENDIF

C.............  Unless ICOD is set to apply to all relevant processes, then
C.............  set parameters to expand xref records to several emission type
            ELSE

                EMT = INDEX1 ( PROC, NEMODE, MEMODE )

C.................  Skip if process is not used
                IF( EMT .LE. 0 ) THEN
                    GO TO 44

C.................  For exhaust, expand to all valid inv pollutants
                ELSEIF( PROC .EQ. 'EXH' ) THEN
                    N1 = 1
                    N2 = NIPOL

C.................  For non-exhaust, use volatile pollutant, if valid
                ELSEIF( VID .GT. 0 ) THEN
                    N1 = VID
                    N2 = VID

C.................  Otherwise, record is useless
                ELSE            
                    GO TO 44

                ENDIF

            ENDIF

C.............  Inspect road class index 
            IF( RDT .NE. 0 ) THEN  ! Make sure we keep the defaults
                NRDT = FIND1( RDT, NRCLAS, MROADS3 )
            ELSE
                NRDT = 1
            ENDIF

C.............  Warn and skip profile if MV-ASCT is blank
            IF( VHT .LT. 0 .AND. RDT .LT. 0 .AND. EMT .EQ. 0 ) THEN
 
                WRITE( MESG,94010 ) 'Temporal x-ref record', IREC,
     &                              'has blank motor-vehicles code.'
                CALL M3MESG( MESG )
 
                GO TO 44   !  to head of the XDEV-read loop

C.............  Skip profile if road class not used
            ELSE IF( NRDT .LT. 0 ) THEN
                GO TO 44   !  to head of the XDEV-read loop
 
C.............  Store XREF entry
            ELSE


                DO 47 I = N1, N2

                    ICOD = INDEX1( MEMODE( EMT ) // '_' // 
     &                             EINAM( I ), NETYPE, EMTNAM )

C.....................  Check table size for overflow
                    IF ( K .GE. MXREF ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                      'Max XREF table size ', MXREF,
     &                      'exceeded at line', IREC
                        CALL M3MESG( MESG )
 
                    ELSE

                        K = K + 1

                        INDXA( K ) = K
                        FIPFA( K ) = STR2INT( LINE( 1 : 5 ) )
                        RDTFA( K ) = RDT
                        LNKFA( K ) = LNK
                        VHTFA( K ) = VHT
                        CODEA( K ) = ICOD
                        MPRFA( K ) = STR2INT( LINE( 38:40 ) )
                        WPRFA( K ) = STR2INT( LINE( 42:44 ) )
                        DPRFA( K ) = STR2INT( LINE( 45:48 ) )

                    ENDIF   ! x-ref overflow or not

   47          CONTINUE

            ENDIF

        GO TO 44  ! to head of XDEV read loop

   50   CONTINUE

        IF ( EFLAG ) THEN
            MESG = 'Error reading TEMPORAL XREF file'
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF

        NPRFA = K

C.......   Process cross-reference file
        WRITE( *,92000 ) ' ', 'Processing TEMPORAL XREF file... ', ' '

C.......   Sort and Process temporal cross-references, according to category
        CALL SORTI4( NPRFA, INDXA, FIPFA, RDTFA, LNKFA, VHTFA )

C.......   Initialize category-specific xref profile number tables
        DO 56 J = 1, NETYPE
        DO  55  I = 1, NVTYPE
            MPR00( I,J ) = IMISS3
            WPR00( I,J ) = IMISS3
            DPR00( I,J ) = IMISS3
            DO  51 S = 1, NRCLAS
                MPRNF( S,I,J ) = IMISS3
                WPRNF( S,I,J ) = IMISS3
                DPRNF( S,I,J ) = IMISS3
51          CONTINUE
            DO  52 S = 1, MXSTRC
                MPRNC( S,I,J ) = IMISS3
                WPRNC( S,I,J ) = IMISS3
                DPRNC( S,I,J ) = IMISS3
52          CONTINUE
            DO  53 S = 1, MXFRC
                MPRFC( S,I,J ) = IMISS3
                WPRFC( S,I,J ) = IMISS3
                DPRFC( S,I,J ) = IMISS3
53          CONTINUE
            DO  54  S = 1, MXFLK
                MPRFL( S,I,J ) = IMISS3
                WPRFL( S,I,J ) = IMISS3
                DPRFL( S,I,J ) = IMISS3
54          CONTINUE
55      CONTINUE
56      CONTINUE

C.........  Group xref profiles by how specifically they are referenced.

        EFLAG = .FALSE.         !  for detecting XREF errors
        DO 66  K = 1, NPRFA

            J    = INDXA( K )
            FIP  = FIPFA( J )
            RDT  = RDTFA( J )
            VHT  = VHTFA( J )
            LNK  = LNKFA( J )
            ICOD = CODEA( J )
            IMON = MPRFA( J )
            IWEK = WPRFA( J )
            IDIU = DPRFA( J )

            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles
 
                IF ( RDT .EQ. 0 ) THEN  !  ultimate fallback profiles
 
                    CALL GRPXREF( VHT, ICOD, IMON, IWEK, IDIU, 
     &                      MPR00, WPR00, DPR00, 1, NVTYPE, MXEMIS,
     &                      1, NVTYPE, NETYPE )

                ELSE                    !  RDT?-dependent-only profiles
 
                    IF( NPRNF .EQ. 0              .OR.        ! First instance
     &                  RDT   .NE. RDTNF( NPRNF )      )      ! Not repeat
     &                  NPRNF = NPRNF + 1

                    IF ( NPRNF .LE. NRCLAS ) THEN

                        RDTNF( NPRNF ) = RDT

                        CALL GRPXREF( VHT, ICOD, IMON, IWEK, IDIU, 
     &                              MPRNF, WPRNF, DPRNF, NRCLAS, NVTYPE, 
     &                              MXEMIS, NPRNF, NVTYPE, NETYPE )
                    END IF
                END IF

            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     ! State dependent
 
                FIP = FIP/1000
 
                IF( NPRNC .EQ. 0              .OR.         ! First instance
     &              FIP   .NE. FIPNC( NPRNC ) .OR.         ! Not previous
     &              RDT   .NE. RDTNC( NPRNC )      )       ! Not previous
     &              NPRNC = NPRNC + 1

                IF ( NPRNC .LE. MXSTRC ) THEN

                    FIPNC( NPRNC ) = FIP
                    RDTNC( NPRNC ) = RDT

                    CALL GRPXREF( VHT, ICOD, IMON, IWEK, IDIU, 
     &                          MPRNC, WPRNC, DPRNC, MXSTRC, NVTYPE,
     &                          MXEMIS, NPRNC, NVTYPE, NETYPE )

                END IF
 
            ELSE IF ( LNK .EQ. 0 ) THEN    !  FIP-RDT dependent profiles
 
                IF ( NPRFC .EQ. 0              .OR.        ! First instance
     &               FIP   .NE. FIPFC( NPRFC ) .OR.        ! Not previous
     &               RDT   .NE. RDTFC( NPRFC )      )      ! Not previous
     &               NPRFC = NPRFC + 1

                IF ( NPRFC .LE. MXFRC ) THEN

                    FIPFC( NPRFC ) = FIP
                    RDTFC( NPRFC ) = RDT
 
                    CALL GRPXREF( VHT, ICOD, IMON, IWEK, IDIU, 
     &                          MPRFC, WPRFC, DPRFC, MXFRC, NVTYPE,
     &                          MXEMIS, NPRFC, NVTYPE, NETYPE )

                END IF

            ELSE                 !  FIP-LNK dependent profiles
 
                IF ( NPRFL .EQ. 0              .OR.        ! First instance
     &               FIP   .NE. FIPFL( NPRFL ) .OR.        ! Not previous
     &               LNK   .NE. LNKFL( NPRFL )      )      ! Not previous
     &               NPRFL = NPRFL + 1

                IF ( NPRFL .LE. MXFLK ) THEN

                    FIPFL( NPRFL ) = FIP
                    RDTFL( NPRFL ) = RDT
                    LNKFL( NPRFL ) = LNK
 
                    CALL GRPXREF( VHT, ICOD, IMON, IWEK, IDIU, 
     &                          MPRFL, WPRFL, DPRFL, MXFLK, NVTYPE,
     &                          MXEMIS, NPRFL, NVTYPE, NETYPE )

                END IF

            ENDIF

   66   CONTINUE   ! End of 1st XREF processing loop

C...........   Check for overflows
        IF( NPRNF .GT. NRCLAS ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max XREF NRCLAS=', NRCLAS,
     &             'exceeded in XREF file: count=', NPRNF
            CALL M3MESG( MESG )
        ENDIF

        IF( NPRNC .GT. MXSTRC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max XREF MXSTRC=', MXSTRC,
     &             'exceeded in XREF file: count=', NPRNC
            CALL M3MESG( MESG )
        ENDIF

        IF( NPRFC .GT. MXFRC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max XREF MXFRC=', MXFRC,
     &             'exceeded in XREF file: count=', NPRFC
            CALL M3MESG( MESG )
        ENDIF

        IF( NPRFL .GT. MXFLK ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max XREF MXFLK=', MXFLK,
     &             'exceeded in XREF file: count=', NPRFL
            CALL M3MESG( MESG )
        ENDIF
 
C...........   Exit if overflows found
        IF ( EFLAG ) THEN
            MESG = 'Error processing TEMPORAL XREF file'
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF

C...........   Construct mapping of source number to effective
C...........   composite temporal profile number and to vehicle mixes.

        DO 299 S = 1, NMSRC

            MFLAG = ( MOD( TFLAG( S ), MTPRFAC ) .EQ. 0 )
            WFLAG = ( MOD( TFLAG( S ), WTPRFAC ) .EQ. 0 )
           
            IF ( .NOT. MFLAG ) THEN
                DO 88 M = 1, NETYPE
                DO 87 V = 1, NVTYPE
                    MDEX( S,V,M ) = 0
   87           CONTINUE
   88           CONTINUE
            END IF
 
            IF ( .NOT. WFLAG ) THEN
                DO 99 M = 1, NETYPE
                DO 98 V = 1, NVTYPE
                    WDEX( S,V,M ) = 0
   98           CONTINUE
   99           CONTINUE
            END IF

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )
 
C...........   Search for FIP-LNK match; else for FIP-RDT match;
C...........   else for STATE-RDT match; else for     RDT match;
C...........   else use ultimate fall back default.

C............. Set using FIP-RDT-LNK specific profiles
            F = FIND3( FIP, RDT, LNK, NPRFL, FIPFL, RDTFL, LNKFL )
            IF ( F .GT. 0 ) THEN
                DO 101 J = 1, NETYPE
                    DO 100 I = 1, NVTYPE
                        MREF( I,J ) = MPRFL( F,I,J )
                        WREF( I,J ) = WPRFL( F,I,J )
                        DREF( I,J ) = DPRFL( F,I,J )
  100               CONTINUE
  101           CONTINUE
                GO TO 266

            ENDIF

C............. Set using FIP-RDT specific profiles
            F = FIND2( FIP, RDT, NPRFC, FIPFC, RDTFC )
            IF ( F .GT. 0 ) THEN
                DO 103 J = 1, NETYPE
                    DO 102 I = 1, NVTYPE
                        MREF( I,J ) = MPRFC( F,I,J )
                        WREF( I,J ) = WPRFC( F,I,J )
                        DREF( I,J ) = DPRFC( F,I,J )
  102               CONTINUE
  103           CONTINUE
                GO TO 266

            END IF
 
C............. Set using STATE-RDT specific profiles
            F = FIND2( FIP/1000, RDT, NPRNC, FIPNC, RDTNC )
            IF ( F .GT. 0 ) THEN
                DO 105 J = 1, NETYPE
                    DO 104 I = 1, NVTYPE
                        MREF( I,J ) = MPRNC( F,I,J )
                        WREF( I,J ) = WPRNC( F,I,J )
                        DREF( I,J ) = DPRNC( F,I,J )
  104               CONTINUE
  105           CONTINUE
                GO TO 266

            END IF
 
C............. Set using STATE specific profiles
            F = FIND2( FIP/1000, 0, NPRNC, FIPNC, RDTNC )
            IF ( F .GT. 0 ) THEN
                DO 107 J = 1, NETYPE
                    DO 106 I = 1, NVTYPE
                        MREF( I,J ) = MPRNC( F,I,J )
                        WREF( I,J ) = WPRNC( F,I,J )
                        DREF( I,J ) = DPRNC( F,I,J )
  106               CONTINUE
  107           CONTINUE
                GO TO 266

            END IF
 
C............. Set using RDT specific profiles
            F = FIND1( RDT, NPRNF, RDTNF )    ! RDT search
            IF ( F .GT. 0 ) THEN
                DO 109 J = 1, NETYPE
                    DO 108 I = 1, NVTYPE
                        MREF( I,J ) = MPRNF( F,I,J )
                        WREF( I,J ) = WPRNF( F,I,J )
                        DREF( I,J ) = DPRNF( F,I,J )
  108               CONTINUE
  109           CONTINUE

C............. Set using fallback profiles
            ELSE 
                DO 111 J = 1, NETYPE
                    DO 110 I = 1, NVTYPE
                        MREF( I,J ) = MPR00( I,J )
                        WREF( I,J ) = WPR00( I,J )
                        DREF( I,J ) = DPR00( I,J )
  110               CONTINUE
  111           CONTINUE

            END IF

266         CONTINUE  ! End of XREF selection

C............. Now look up these refs in profiles table
            DO  277  M = 1, NETYPE
            DO  276  V = 1, NVTYPE
 
                IF ( MFLAG ) THEN
                    MDEX( S,V,M ) = FIND1( MREF( V,M ), NMON, MONREF )
                    IF( MDEX( S,V,M ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'XREF index not found for MREF=', MREF( V,M )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF
 
                IF ( WFLAG ) THEN
                    WDEX( S,V,M ) = FIND1( WREF( V,M ), NWEK, WEKREF )
                    IF( WDEX( S,V,M ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'XREF index not found for WREF=', WREF( V,M )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF
 
                DDEX( S,V,M ) = FIND1( DREF( V,M ), NDIU, DIUREF )
                IF( DDEX( S,V,M ) .LT. 0 .AND. NDIU .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'XREF index not found for DREF=', DREF( V,M )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
 
                EDEX( S,V,M ) = FIND1( DREF( V,M ), NEND, ENDREF )
                IF( EDEX( S,V,M ) .LT. 0 .AND. NEND .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'XREF index not found for EREF=', DREF( V,M )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
 
  276       CONTINUE!  end index lookup loop on species V
  277       CONTINUE!  end index lookup loop on emissions type M

  299   CONTINUE!  end loop on sources

C...........   Read and process emissions factors cross reference file MPLIST

        WRITE( *,92000 ) 
     &         ' ', 'Reading EMISSION FACTORS XREF file...', ' '

        NMPLIST = RDMPLIST( FDEV, MPFIP, MPRDT, MPLNK,
     &                      MPREF, REFCNT, REFLIST )
        IF ( NMPLIST .LE. 0 ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0,0, 'Error reading MPLIST', 2 )
        END IF

C...........   Map emissions factor indexes onto sources:

        EFLAG = .FALSE.
        DO  355  S = 1, NMSRC

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )

C............   Link-specific matches:

            F = FIND3( FIP, RDT, LNK, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  301  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
301             CONTINUE
                GO TO  355
            END IF

C............   FIP-RDT matches:

            F = FIND3( FIP, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  311  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
311             CONTINUE
                GO TO  355
            END IF

C............   state--RDT matches:

            FIP = 1000 * ( FIP / 1000 )
            F = FIND3( FIP, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  322  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
322             CONTINUE
                GO TO  355
            END IF

C............   State-only matches:
            F = FIND3( FIP, 0, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  327  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
327             CONTINUE
                GO TO  355
            END IF

C............   RDT-only matches:

            F = FIND3( 0, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  333  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
333             CONTINUE
                GO TO  355
            END IF

C............   fallback

            F = FIND3( 0, 0, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  344  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
344             CONTINUE
                GO TO  355
            END IF

C...........   If you get to here:  error

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &            'No emission factor xref for FIP', FIP,
     &            'RCLASS', RDT, 'Link', LNK
            CALL M3MESG( MESG )

355     CONTINUE

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Bad emission factor matchup', 2 )
        END IF                                               


C...........   Read non-diurnal emissions factors params from MEFSND

        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description for ' // FNAME
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF

C.......   Get starting point, increment, number of items for
C.......   non-diurnal emissions factor temperature tables:

        NTEMP  = STR2INT( FDESC3D( 5 ) )
        TEMP0  = CTOK  +  SNGL( FTOC * ( XORIG3D - 32.0D0 ) )
        DTEMP  = SNGL( FTOC * XCELL3D )
        TEMP1  = TEMP0 +  DTEMP * FLOAT( NTEMP - 1 )
        DDTEMP = 1.0 / DTEMP

C.......   Get constants and increment for
C.......   non-diurnal emissions factor temperature tables:

        MAXTL  = INT( ( MAXT_MAX - MAXT_MIN ) / TINTV ) + 1 
        DNLTI  = 1.0 / ( FTOC * TINTV )

C.......   Set date for reading emission factors
        FDATE = JYEAR * 1000 + 1

C.......   Read emissions factor tables:

        DO  366  T = 1, REFCNT

            IF ( .NOT. READ3( FNAME, 'EXH_CO', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_CO( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_CO from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EXH_NOX', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_NOX( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_NOX from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EXH_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EVP_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EVP_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EVP_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'RNL_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        RNL_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable RNL_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'RST_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        RST_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable RST_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( DNAME, 'WDL_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        WDL_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable WDL_VOC from ' // DNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( DNAME, 'DNL_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        DNL_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable DNL_VOC from ' // DNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( DNAME, 'HOT_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        HOT_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable HOT_VOC from ' // DNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( DNAME, 'CRC_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        CRC_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable CRC_VOC from ' // DNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

366     CONTINUE

C.......   Read ungridding matrix:

        IF ( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 0, 0, NU ) ) THEN
            MESG = 'Could not read ungridding matrix from ' // UNAME
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF

C......................................................................
C...........   Construct TMAT -- array of effective composite 
C...........   profile coefficients
C......................................................................

C...........   For each time step, write out the area source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out emissions values...', ' '

        LDATE = -1   !  initialization forces TMAT calculation on first iter.

        DO  899  T = 1, NSTEPS

            IF ( .NOT. READ3( MNAME, 'TA', 1, MDATE, MTIME, TA ) ) THEN
                MESG = 'Could not read TA from ' // MNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

C...............   If this is a new month, or new day, calculate
C...............   new TMAT

            IF ( LDATE .NE. JDATE ) THEN

                LDATE = JDATE

C...............   Adjust sources' time zones to account for daylight time
C..................   Add 1 if date is daylight and TZONES is not
C..................   Subtract 1 if date is not daylight and TZONES is daylight
                IF( ISDSTIME( JDATE ) .AND. .NOT. DAYLIT ) THEN
 
                    DAYLIT = .TRUE.
 
                    DO 777 S = 1, NMSRC
                        TZONES( S ) = TZONES( S ) + 1
  777               CONTINUE
 
                ELSEIF( .NOT. ISDSTIME( JDATE ) .AND. DAYLIT ) THEN
 
                    DAYLIT = .FALSE.
 
                    DO 778 S = 1, NMSRC
                        TZONES( S ) = TZONES( S ) - 1
  778               CONTINUE
 
                ENDIF

C.................  Create temporal adjustment matrix for current day

                CALL DAYMON( JDATE, MONTH, DAY )  !  get month
                DAY = WKDAY( JDATE )              !  get day-of-week

                IF ( NEND .EQ. 0 ) THEN  !  no weekend packet

                    WRITE( *,92011 ) 'Processing day', DAY
                    CALL MKVMAT( NMSRC,  NETYPE, JDATE, MONTH,  DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, NMON, MDEX,
     &                           WEKFAC, NWEK, WDEX,
     &                           DIUFAC, NDIU, DDEX, MXMAT, TMAT )

                ELSE IF ( DAY .GE. 6 ) THEN     !  day is Saturday or Sunday
                                                !  in SAI's brain-damaged
                                                !  day-of-week conventions.

                    WRITE( *,92011 ) 'Processing weekend day', DAY
                    CALL MKVMAT( NMSRC,  NETYPE, JDATE, MONTH, DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, NMON, MDEX,
     &                           WEKFAC, NWEK, WDEX,
     &                           ENDFAC, NEND, EDEX, MXMAT, TMAT )

                ELSE

                    WRITE( *,92011 ) 'Processing midweek day', DAY
                    CALL MKVMAT( NMSRC,  NETYPE, JDATE, MONTH, DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, NMON, MDEX,
     &                           WEKFAC, NWEK, WDEX,
     &                           DIUFAC, NDIU, DDEX, MXMAT, TMAT )

                END IF


C.................  Read ungridded MIN/MAX temperature indices for the day

                IF( .NOT. READ3( 
     &              WNAME, ALLVAR3, 1, WDATE, WTIME, TMMI ) ) THEN

                     MESG = 'Could not read TMMI from ' // 
     &                      WNAME( 1:TRIMLEN( WNAME ) )
                     CALL M3EXIT( 'TMPMOBIL', MDATE, 0, MESG, 2 )
                ENDIF
             

            END IF      ! if jdate not ldate

            H  = 1 + MOD( JTIME / 10000 , 24 )
            KU = 0
            OFLAG = .FALSE.
            OSRC = 0

            DO  888  S = 1, NMSRC

C...............   Ungrid temperature:

                TASRC = 0.0
                DO  811  J = 1, NU( S )
                    KU = KU + 1
                    TASRC = TASRC + TA( IU( KU ) ) * CU( KU )
811             CONTINUE

C...............   Get factors for interpolating non-diurnal emission factors

                IF( NU( S ) .EQ. 0 ) THEN          !  source outside grid
 
                    OFLAG = .TRUE.
                    OSRC = OSRC + 1
 
                    PP = 0.0
                    QQ = 0.0
                    II = 1
 
                ELSEIF ( TASRC .LT. TEMP0 ) THEN   !  trap TASRC against TEMP0

                    WRITE( MESG,94030 )
     &                  'Temperature', TASRC, 'at source', S,
     &                  'is below MEFS table min:', TEMP0
                    CALL M3WARN( 'TMPMOBIL', JDATE, JTIME, MESG )
                    PP = 1.0
                    QQ = 0.0
                    II = 1

                ELSE IF ( TASRC .GT. TEMP1 ) THEN   !  and against TEMP1

                    WRITE( MESG,94030 )
     &                  'Temperature', TASRC, 'at source', S,
     &                  'is above MEFS table max:', TEMP1
                    CALL M3WARN( 'TMPMOBIL', JDATE, JTIME, MESG )
                    PP = 0.0
                    QQ = 1.0
                    II = NTEMP - 1

                ELSE			!  get table entry, interp coeffs

                    SCR = DDTEMP * ( TASRC - TEMP0 )
                    II  = 1 + INT( SCR )
                    QQ  = AMOD( SCR, 1.0 )
                    PP  = 1.0 - QQ

                END IF		!  tasrc trapped; coeffs & index found.

C...............   Get factors for interpolating diurnal emission factors
C...................   Note that min/max has already been trapped by PREDIUR

                J1     = TMMI( S )

C.................  Set up if diurnal emissions are non-zero
                IF( J1 .NE. 0 ) THEN

                    J2     = J1 + MAXTL
                    J3     = J1 + 1
                    J4     = J2 + 1

                    TMMIDV = INT( J1 / MAXTL )
                    TMIN0  = MINT_MIN + TINTV * TMMIDV
                    TMIN0  = CTOK  +  SNGL( FTOC * ( TMIN0 - 32.0D0 ) )
                    TMAX0  = MAXT_MIN + TINTV * ( J1 - MAXTL*TMMIDV - 1)
                    TMAX0  = CTOK  +  SNGL( FTOC * ( TMAX0 - 32.0D0 ) )

                    R1     = AMOD ( DNLTI * ( TKMIN(S) - TMIN0 ), 1.0 )
                    R2     = 1.0 - R1
                  
                    S1     = AMOD ( DNLTI * ( TKMAX(S) - TMAX0 ), 1.0 )
                    S2     = 1.0 - S1

C.................  Set up if diurnal emissions are zero
                ELSE
                    J1 = 1
                    J2 = 1
                    J3 = 1
                    J4 = 1
                    R1 = 0.
                    R2 = 0.
                    S1 = 0.
                    S2 = 0.

                ENDIF
                  
C...............   Adjust which hourly EF to use based on time zone

                HCORR = MOD( H + TZONE + 25 - TZONES( S ), 24 )
                M     = EFDEX( HCORR, S )

                I = 1 

                IF( EMTNAM( I ) .EQ. 'EXH_CO' ) THEN

                    DO  818  V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * EXH_CO( II  ,V,M ) +
     &                                QQ * EXH_CO( II+1,V,M ) )
  818               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'EXH_NOX' ) THEN

                    DO 819 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * EXH_NOX( II  ,V,M ) +
     &                                QQ * EXH_NOX( II+1,V,M ) )
  819               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'EXH_VOC' ) THEN

                    DO 820 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * EXH_VOC( II  ,V,M ) +
     &                                QQ * EXH_VOC( II+1,V,M ) )
  820               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'EVP_VOC' ) THEN

                    DO 821 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * EVP_VOC( II  ,V,M ) +
     &                                QQ * EVP_VOC( II+1,V,M ) )
  821               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'RNL_VOC' ) THEN

                    DO 822 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * RNL_VOC( II  ,V,M ) +
     &                                QQ * RNL_VOC( II+1,V,M ) )
  822               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'RST_VOC' ) THEN

                    DO 823 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                              ( PP * RST_VOC( II  ,V,M ) +
     &                                QQ * RST_VOC( II+1,V,M ) )
  823               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'WDL_VOC' ) THEN

                    DO 824 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                             ( S2 * ( R2 * WDL_VOC( J1,V,M ) +
     &                                      R1 * WDL_VOC( J2,V,M )   ) +
     &                               S1 * ( R2 * WDL_VOC( J3,V,M ) +
     &                                      R1 * WDL_VOC( J4,V,M )   ) )
  824               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'DNL_VOC' ) THEN

                    DO 825 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                             ( S2 * ( R2 * DNL_VOC( J1,V,M ) +
     &                                      R1 * DNL_VOC( J2,V,M )   ) +
     &                               S1 * ( R2 * DNL_VOC( J3,V,M ) +
     &                                      R1 * DNL_VOC( J4,V,M )   ) )
  825               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'HOT_VOC' ) THEN

                    DO 826 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                             ( S2 * ( R2 * HOT_VOC( J1,V,M ) +
     &                                      R1 * HOT_VOC( J2,V,M )   ) +
     &                               S1 * ( R2 * HOT_VOC( J3,V,M ) +
     &                                      R1 * HOT_VOC( J4,V,M )   ) )
  826               CONTINUE

                    I = I + 1
                ENDIF

                IF( EMTNAM( I ) .EQ. 'CRC_VOC' ) THEN

                    DO 827 V = 1, NVTYPE

                        TOUT( S,V,I ) = TMAT( S,V,I,H ) *
     &                             ( S2 * ( R2 * CRC_VOC( J1,V,M ) +
     &                                      R1 * CRC_VOC( J2,V,M )   ) +
     &                               S1 * ( R2 * CRC_VOC( J3,V,M ) +
     &                                      R1 * CRC_VOC( J4,V,M )   ) )
  827               CONTINUE

                    I = I + 1
                ENDIF

  888       CONTINUE		!  end loop on sources S

            IF ( .NOT.WRITE3( TNAME, 'ALL', JDATE, JTIME, TOUT ) ) THEN

                CALL M3EXIT( 'TMPMOBIL', JDATE, JTIME,
     &                       'Could not write timestep to "' // 
     &                       TNAME( 1 : TRIMLEN( TNAME ) ) // '".', 2 )

            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )
            CALL NEXTIME( MDATE, MTIME, TSTEP )

899     CONTINUE

        IF( OFLAG ) THEN

            WRITE( MESG, 94010 ) 
     &         'During ungridding, ', OSRC,' sources excluded from grid'

            CALL M3WARN( 'TMPMOBIL', 0, 0, MESG )

        ENDIF

999   CONTINUE          !  exit program

      CALL M3EXIT( 'TMPMOBIL',0,0, 'Normal completion of TMPMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, 1X, I7.7, 1X, I6.6, 1X, I2 )

92011   FORMAT( /5X, A, I10 )

92020   FORMAT( 5X, A, I9 )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I5, 1X, A, I7, 1X, 10( A, I7, :, 1X ) )

94030   FORMAT( 10 ( A, :, 1PG12.5, :, 2X , A, :, I5, :, 2X ) )

        END

