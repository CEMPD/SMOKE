
        PROGRAM LAYPOINT

C***********************************************************************
C  program body starts at line 262
C
C  DESCRIPTION:
C     This program computes the layer fractions for point sources.  It uses
C     a modified Briggs algorithm to compute plume rise.  The plume is
C     allocated to multiple layers when necessary.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKMETEM
        CHARACTER*2     CRLF
        LOGICAL         DSCM3GRD
        LOGICAL         DSCM3LAY
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FIND1
        CHARACTER*50    GETCFDSC
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         WKDAY

        EXTERNAL   CHKMETEM, CRLF, DSCM3GRD, DSCM3LAY, ENVINT, ENVYN, 
     &             FIND1, GETCFDSC, HHMMSS, INDEX1, MMDDYY, 
     &             PROMPTFFILE, PROMPTMFILE, VERCHAR, WKDAY

C...........  LOCAL PARAMETERS and their descriptions:

        REAL        , PARAMETER :: USTARMIN  = 0.1  ! Min valid value for USTAR

        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C.........  Indicator for which public inventory arrays need to be read
        INTEGER               , PARAMETER :: NINVARR = 8
        CHARACTER(LEN=IOVLEN3), PARAMETER :: IVARNAMS( NINVARR ) = 
     &                                 ( / 'IFIP           '
     &                                   , 'XLOCA          '
     &                                   , 'YLOCA          '
     &                                   , 'STKHT          '
     &                                   , 'STKDM          '
     &                                   , 'STKTK          '
     &                                   , 'STKVE          '
     &                                   , 'CSOURC         ' / )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point source stack parameters:

C.........  Allocatable, per-source meteorology variables
        REAL   , ALLOCATABLE :: HFX  ( : )    !  sensible heat flux (M K / S )
        REAL   , ALLOCATABLE :: HMIX ( : )    !  mixing height (m)
        REAL   , ALLOCATABLE :: TSFC ( : )    !  surface temperature (deg K)
        REAL   , ALLOCATABLE :: USTAR( : )    !  friction velocity (m/s)

C.........  Allocatable, per-source and per layer meteorology variables. 
C           Dimensioned by layers, then sources
        REAL   , ALLOCATABLE :: DDZH ( :,: )  !  1/( zh(l) - zh(l-1) )
        REAL   , ALLOCATABLE :: DDZF ( :,: )  !  1/( zh(l) - zh(l-1) )
        REAL   , ALLOCATABLE :: PRES ( :,: )  !  pressure (Pa)
        REAL   , ALLOCATABLE :: QV   ( :,: )  !  mixing ratio (kg/kg)
        REAL   , ALLOCATABLE :: TA   ( :,: )  !  temperature (K)
        REAL   , ALLOCATABLE :: UWIND( :,: )  !  wind speed (m/s)
        REAL   , ALLOCATABLE :: VWIND( :,: )  !  wind speed (m/s)
        REAL   , ALLOCATABLE :: ZF   ( :,: )  !  layer surface height (m)
        REAL   , ALLOCATABLE :: ZH   ( :,: )  !  layer center  height (m)
        REAL   , ALLOCATABLE :: ZSTK ( :,: )  !  zf( l,s ) - stkht(s) (m)

C.........  Allocatable, temporary per-layer variables from 1:EMLAYS
        REAL   , ALLOCATABLE :: WSPD ( : )    !  wind speed (m/s)
        REAL   , ALLOCATABLE :: DTHDZ( : )    !  gradient of THETV
        REAL   , ALLOCATABLE :: TFRAC( : )    !  Temporary LFRAC

C.........  Allocatable, temporary per-layer variables from 0:EMLAYS
        REAL   , ALLOCATABLE :: PRESF( : )    !  pressure at full-levels
        REAL   , ALLOCATABLE :: ZZF  ( : )    !  elevation at full-levels

C.........  Allocatable cross- OR dot-point meteorology input buffers
        REAL   , ALLOCATABLE :: XBUF( :,: )   ! cross-point
        REAL   , ALLOCATABLE :: DBUF( :,: )   ! dor-point

C.........  Allocatable un-gridding matrices (uses bilinear interpolation)
C           Dimensioned 4 by NSRC
        INTEGER, ALLOCATABLE :: ND( :,: )     !  dot-point, cell indexes
        INTEGER, ALLOCATABLE :: NX( :,: )     !  cross-point, cell indexes
   
        REAL   , ALLOCATABLE :: CD( :,: )     !  dot-point, coefficients
        REAL   , ALLOCATABLE :: CX( :,: )     !  cross-point, coefficients

C.........  Output layer fractions, dimensioned NSRC, emlays
        REAL   , ALLOCATABLE :: LFRAC( :, : )

C.........  Input/output hour-specific data index, dimensioned by NSRC and
C           by EMLAYS, so that index can be written to PLAY_EX file
        INTEGER, ALLOCATABLE :: LOCINDXH( :,: )

C.........  Fixed-dimension arrays
        REAL         LFULLHT( 0:MXLAYS3 )     !  full-level heights [m]
        REAL         LHALFHT( 1:MXLAYS3 )     !  half-level heights [m]
        REAL         SIGH   ( 0:MXLAYS3-1 )   !  half-level sigma values
        REAL         VGLVSXG( 0:MXLAYS3 )     !  vertical coord values
        REAL         WEIGHTS( 1:MXLAYS3 )     !  tmp weights for vertical aloc

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  log file
        INTEGER      :: PDEV = 0!  elevated/PinG source file
        INTEGER      :: RDEV = 0!  optional report iff REP_LAYER_MAX is set
        INTEGER         SDEV    !  ASCII part of inventory file

        CHARACTER*16    ANAME   !  ASCII point-source inventory file
        CHARACTER*16    DNAME   !  dot-point layered met file name
        CHARACTER*16    ENAME   !  point-source inventory input file
        CHARACTER*16    GNAME   !  cross-point layered grid file name
        CHARACTER*16    HNAME   !  hourly input file name
        CHARACTER*16    LNAME   !  layer fractions matrix output file
        CHARACTER*16    SNAME   !  cross-point surface met file name
        CHARACTER*16    XNAME   !  cross-point layered met file name

C...........   Other local variables

        INTEGER          I, J, K, L, L1, L2, S, T  ! counters and indices

        INTEGER          EMLAYS    ! number of emissions layers
        INTEGER          IOS       ! tmp i/o status
        INTEGER          JDATE     ! Julian date (YYYYDDD)
        INTEGER          JTIME     ! time (HHMMSS)
        INTEGER          LBOT      ! plume bottom layer
        INTEGER          LDATE     ! previous date
        INTEGER          LPBL      ! first L: ZF(L) above mixing layer
        INTEGER          LSTK      ! first L: ZF(L) > STKHT
        INTEGER          LTOP      ! plume top    layer
        INTEGER          METNCOLS  ! met grid number of columns
        INTEGER          METNGRID  ! met grid number of cells
        INTEGER          METNROWS  ! met grid number of rows
        INTEGER          NDOTS     ! dot grid number of cells
        INTEGER          NHR       ! no. hour-specific sources for current hour
        INTEGER          NMAJOR    ! no. major sources
        INTEGER          NPING     ! no. plume-in-grid sources
        INTEGER       :: NSTEPS= 1 ! mumber of time steps
        INTEGER          REP_LAYR  ! layer for reporting srcs w/ high plumes
        INTEGER       :: SDATE = 0 ! Julian start date (YYYYDDD)
        INTEGER       :: STIME = 0 ! start time (HHMMSS)
        INTEGER          TSTEP     ! output time step

        REAL             X, Y, P, Q, PP, QQ
        REAL             DM, HT, TK, VE, FL  ! tmp stack parameters
        REAL             XBEG, XEND, XL  ! tmp x-coords
        REAL             YBEG, YEND, YL  ! tmp y-coords
        REAL             FAC       !  tmp factor for renormalizing
        REAL             PSFC      !  surface pressure (Pa)
        REAL             SURFACE   !  tmp weight at surface
        REAL             TDIFF     !  tmp layer frac diff for renormalizing
        REAL             TSTK      !  temperature at top of stack (K)
        REAL             TSUM      !  tmp layer frac sum for renormalizing
        REAL             USTMP     !  tmp Ustar
        REAL             WSTK      !  wind speed  at top of stack (m/s)
        REAL             ZZ0, ZZ1, ZF0, ZF1
        REAL             ZBOT      !  plume bottom elevation (m)
        REAL             ZTOP      !  plume top    elevation (m)

        REAL*8           METXORIG  ! cross grid X-coord origin of met grid 
        REAL*8           METYORIG  ! cross grid Y-coord origin of met grid
        REAL*8           XCELLDG   ! dot grid X-coordinate cell dimension
        REAL*8           YCELLDG   ! dot grid Y-coordinate cell dimension
        REAL*8           XORIGDG   ! dot grid X-coordinate origin of grid 
        REAL*8           YORIGDG   ! dot grid Y-coordinate origin of grid

        LOGICAL       :: BFLAG = .FALSE.  ! true: use plume bottom and top
        LOGICAL       :: CFLAG = .FALSE.  ! true: recalc vel w/ flow & diam
        LOGICAL     :: COMPUTE = .FALSE.  ! true: compute plume rise 
        LOGICAL       :: EFLAG = .FALSE.  ! error flag
        LOGICAL       :: FFLAG = .FALSE.  ! true: use hourly flow rate
        LOGICAL       :: HFLAG = .FALSE.  ! true: hourly input used
        LOGICAL       :: IFLAG = .FALSE.  ! true: hr data okay for timestep
        LOGICAL       :: LFLAG = .FALSE.  ! true: use hourly layer 1 fraction
        LOGICAL       :: PFLAG = .FALSE.  ! true: compute plm ris for iteration
        LOGICAL       :: TFLAG = .FALSE.  ! true: use hourly temperatures
        LOGICAL       :: VFLAG = .FALSE.  ! true: use elevated file (PELV)
        LOGICAL       :: XFLAG = .FALSE.  ! true: process ONLY explicit sources
        LOGICAL       :: YFLAG = .FALSE.  ! true: use hourly velocities
        LOGICAL       :: ZSTATIC = .TRUE. ! true: Get heights from GRID_CRO file
        LOGICAL          LFG( 9 )          ! true: source characteristic is valid

        CHARACTER*50     CHARS( 9 )!  tmp source characeristics 
        CHARACTER*50  :: METSCEN   !  temporary string for met scenario name
        CHARACTER*50  :: CLOUDSHM  !  temporary string for cloud scheme name
        CHARACTER*80  :: GDESC     !  grid description
        CHARACTER*256    OUTFMT    !  output format for RDEV report
        CHARACTER*256    BUFFER    !  source characteristics buffer
        CHARACTER*256    MESG      !  buffer for M3EXIT() messages

        CHARACTER(LEN=IOVLEN3) VNAME      ! variable name buffer 
        CHARACTER(LEN=IOVLEN3) COORD3D    ! coordinate system name
        CHARACTER(LEN=IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER*16  :: PROGNAME = 'LAYPOINT'   !  program name

C***********************************************************************
C   begin body of program LAYPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........   Get setting from environment variables
        EMLAYS = ENVINT( 'SMK_EMLAYS', 'Number of emission layers',
     &                   -1, IOS )

        MESG = 'Use Elevpoint output to determine elevated sources'
        VFLAG = ENVYN( 'SMK_SPECELEV_YN', MESG, .FALSE., IOS )

        MESG = 'Indicator for defining hourly plume rise data'
        HFLAG = ENVYN( 'HOUR_PLUMEDATA_YN', MESG, .FALSE., IOS )

        MESG = 'Indicator for processing ONLY explicit plume ' //
     &         'rise sources'
        XFLAG = ENVYN( 'EXPLICIT_PLUMES_YN', MESG, .FALSE., IOS )

C.........  Must have HOUR_PLUMEDATA_YN = Y to have EXPLICIT_PLUMES_YN = Y
        IF ( XFLAG .AND. .NOT. HFLAG ) THEN
            HFLAG = .TRUE.
            MESG = 'NOTE: Setting HOUR_PLUMEDATA_YN to Y because '//
     &             'EXPLICIT_PLUMES_YN is Y'
            CALL M3MSG2( MESG )
        END IF

        CFLAG = ENVYN( 'VELOC_RECALC', 
     &                 'Flag for recalculating velocity', .FALSE., IOS )

C.........  Cannot use default and cannot set to less than 4 because of
C           limits of plume rise algorithm
        IF( EMLAYS .LT. 4 ) THEN
            MESG = 'Environment variable SMK_EMLAYS must be set to ' //
     &             'a number from 4 to the ' // CRLF() // BLANK10 //
     &             'number of layers in the meteorology inputs.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF
            
        REP_LAYR = ENVINT( 'REP_LAYER_MAX', 
     &                     'Layer number for reporting high plume rise',
     &                     -1, IOS )

        IF( IOS .EQ. 0 ) THEN

            IF( REP_LAYR .LT. 1 ) THEN

                MESG = 'NOTE: Environment variable REP_LAYR_MAX is ' //
     &                 'less than 1.  Turning off reporting...'
 
            ELSE IF ( REP_LAYR .GT. EMLAYS ) THEN
                WRITE( MESG,94010 ) 
     &                 'NOTE: Environment variable REP_LAYR_MAX is '//
     &                 'greater than the number of emissions ' // 
     &                 CRLF() //BLANK10 // 'layers (', EMLAYS, '). '//
     &                 'Resetting to equal number of emissions layers.'
            END IF

            CALL M3MSG2( MESG )

        END IF

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Make sure only run for point sources
        IF( CATEGORY .NE. 'POINT' ) THEN
            MESG = 'ERROR: ' // PROGNAME( 1:LEN_TRIM( PROGNAME ) ) //
     &             ' is not valid for ' // CATEGORY( 1:CATLEN ) // 
     &             ' sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.......   Get file name; open input point sources, temporal cross-reference,
C.......   and temporal profiles files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINT I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for the ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

C.........  Get header description of inventory file 
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:LEN_TRIM( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO
            IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
            IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

        END IF

C.........  Get file name and open daily input inventory file
        IF( HFLAG ) THEN
            HNAME = PROMPTMFILE( 
     &               'Enter logical name for HOUR-SPECIFIC file',
     &               FSREAD3, CRL // 'HOUR', PROGNAME )

C.............  Check to see if appropriate variable list exists
            CALL RETRIEVE_IOAPI_HEADER( HNAME )

            NHRSRC = NROWS3D

C.............  Check input variables and allocate memory...
C.............  Check for layer-1 fraction
            I = INDEX1( SPDATNAM( 1 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                LFLAG = .TRUE.
                WRITE( MESG,94010 ) 'NOTE: Layer-1 fraction ' //
     &                 'hourly input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( LAY1F( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LAY1F', PROGNAME )

            END IF
            
C.............  Check for plume top and plume bottom
            J = INDEX1( SPDATNAM( 2 ), NVARS3D, VNAME3D )
            K = INDEX1( SPDATNAM( 3 ), NVARS3D, VNAME3D )
            IF ( J .GT. 0 .AND. K .LE. 0 ) THEN
                MESG = 'WARNING: Plume bottom in hourly input file '//
     &                 'will not be used '// CRLF()// BLANK10//
     &                 'because plume top is not also present.'
                CALL M3MSG2( MESG )

            ELSE IF ( J .LE. 0 .AND. K .GT. 0 ) THEN
                MESG = 'WARNING: Plume top in hourly input file '//
     &                 'will not be used '// CRLF()// BLANK10//
     &                 'because plume bottom is not also present.'
                CALL M3MSG2( MESG )

            ELSE IF ( J .GT. 0 .AND. K .GT. 0 ) THEN
                BFLAG = .TRUE.
                WRITE( MESG,94010 ) 'NOTE: Plume top and bottom in ' //
     &                 'hourly input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( PLMBOT( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PLMBOT', PROGNAME )
                ALLOCATE( PLMTOP( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PTOP', PROGNAME )

            END IF

C.............  Check for temperatures
            I = INDEX1( SPDATNAM( 4 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                TFLAG = .TRUE.
                WRITE( MESG,94010 ) 'NOTE: Temperatures ' //
     &                 'hourly input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( HRSTKTK( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'HRSTKTK', PROGNAME )

            END IF
         
C.............  Check for velocity
            I = INDEX1( SPDATNAM( 5 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                YFLAG = .TRUE.
                WRITE( MESG,94010 ) 'NOTE: Velocities ' //
     &                 'hourly input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( HRSTKVE( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'HRSTKVE', PROGNAME )

            END IF
         
C.............  Check for flow rate
            I = INDEX1( SPDATNAM( 6 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                FFLAG = .TRUE.
                WRITE( MESG,94010 ) 'NOTE: Flow rate ' //
     &                 'hourly input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( HRSTKFL( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'HRSTKFL', PROGNAME )

            END IF

C.............  If no correct variables, then ignore file
            HFLAG= ( LFLAG .OR. BFLAG .OR. TFLAG .OR. YFLAG .OR. FFLAG )

C.............  Give warning if no valid data
            IF( .NOT. HFLAG ) THEN
                MESG = 'WARNING: No hourly data used because ' //
     &                 'no correct variables names ' // CRLF() // 
     &                 BLANK10 // '(defined in EMCNST3.EXT) were found.'
                CALL M3MSG2( MESG )
            END IF
         
        END IF      ! End if hourly data use was requested by E.V. settings

        IF( VFLAG ) THEN
            PDEV = PROMPTFFILE( 
     &          'Enter logical name for the ELEVATED POINT SOURCE file',
     &          .TRUE., .TRUE., CRL // 'ELV', PROGNAME )
        END IF

C.........  If not explicit plume rise only, open and process other met files
        IF ( .NOT. XFLAG ) THEN

            SNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT SURFACE MET file',
     &          FSREAD3, 'MET_CRO_2D', PROGNAME )

            GNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED GRID file',
     &          FSREAD3, 'GRID_CRO_3D', PROGNAME )

            XNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED MET file',
     &          FSREAD3, 'MET_CRO_3D', PROGNAME )

            DNAME = PROMPTMFILE( 
     &          'Enter name for DOT-POINT LAYERED MET file',
     &          FSREAD3, 'MET_DOT_3D', PROGNAME )

C.............  Check multiple met files for consistency
            EFLAG = ( .NOT. CHKMETEM( 'NONE',SNAME,GNAME,XNAME,DNAME ) )

            IF ( EFLAG ) THEN

                MESG = 'Input met files have inconsistent grids or ' //
     &                 'layers.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Get grid parameters from 3-d cross-point met file and store 
C               needed header information.  Use time parameters for time 
C               defaults.
            CALL RETRIEVE_IOAPI_HEADER( XNAME )

C.............  Initialize reference grid with met file
            CALL CHKGRID( XNAME, 'GRID', 0, EFLAG )

            SDATE  = SDATE3D
            STIME  = STIME3D
            NSTEPS = MXREC3D
            NLAYS  = NLAYS3D
            VGTYP  = VGTYP3D
            VGTOP  = VGTOP3D
            METNCOLS = NCOLS
            METNROWS = NROWS
            METNGRID = NGRID
            METXORIG = XORIG
            METYORIG = YORIG

            NDOTS = ( NCOLS + 1 ) * ( NROWS + 1 )

            METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. ) 
            CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. ) 

C.........  Determine whether height information is time dependent or time
C           independent. Non-hydrostatic is time-independent and hydrostatic
C           is time-dependent.
            SELECT CASE( VGTYP )
            CASE ( VGSGPH3, VGHVAL3 ) 
                ZSTATIC = .FALSE.

            CASE ( VGSGPN3 )
                ZSTATIC = .TRUE.

            CASE DEFAULT
                WRITE( MESG,94010 ) 'Cannot process vertical ' //
     &                 'coordinate type', VGTYP
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SELECT

            CALL RETRIEVE_IOAPI_HEADER( DNAME )
            XCELLDG = XCELL3D
            YCELLDG = YCELL3D 
            XORIGDG = XORIG3D
            YORIGDG = YORIG3D

C.........  If not using met data (for explicit plume rise only...)
        ELSE
 
C.............  Get vertical layer structure from the G_GRIDPATH file
            IF ( .NOT. DSCM3LAY( NLAYS, VGTYP, VGTOP, VGLVS3D ) )
     &           THEN
                MESG = 'Could not get vertical layer structure from '//
     &                 'Models-3 grid description file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check to make sure input vertical structure has been provided
C               that is "meters above ground."
            IF ( VGTYP .NE. VGHVAL3 ) THEN
                WRITE( MESG,94010 ) 'Explicit plume rise requires ' //
     &                 'vertical type ', VGHVAL3, 'in grid ' //
     &                 'description' // CRLF() // BLANK10 //
     &                 'file, but type', VGTYP, 'was found.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C...........  The following useless do loop is so that the SGI compiler will not
C             give a core dump.  Yes, very strange.
            DO I = 1, NLAYS3D
                J = J + 1
                SIGH   ( I-1 ) = 0.5 * ( VGLVSXG( J ) + VGLVSXG( J-1 ) )
            END DO

        END IF      ! If using met data or not (not only for explicit plumes)

C.........  Get horizontal grid structure from the G_GRIDPATH file
        IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D, COORUN3D,
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,
     &                     YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &                     YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set subgrid if using met files, define grid if not using met files
        CALL CHKGRID( 'GRIDDESC', 'GRIDDESC' , 1, EFLAG ) 

C.........  Store local layer information
        J = LBOUND( VGLVS3D, 1 )
        VGLVSXG( 0 ) = VGLVS3D( J )
        DO I = 1, NLAYS
            J = J + 1
            VGLVSXG( I ) = VGLVS3D( J )
            SIGH   ( I-1 ) = 0.5 * ( VGLVS3D( J ) + VGLVS3D( J-1 ) )
        END DO

C.........  Compare number of meteorology layers to number of emissions layers
        IF( EMLAYS .LE. NLAYS ) THEN
            WRITE( MESG,94010 ) 'NOTE: The number of emission layers '//
     &             'is', EMLAYS, ', and the maximum '// CRLF()// 
     &             BLANK10//'possible layers is', NLAYS
            CALL M3MSG2( MESG )

        ELSE
            WRITE( MESG,94010 ) 'Resetting number of emission layers '//
     &             'from', EMLAYS, 'to number of '// CRLF()// BLANK10 //
     &             'layers in the meteorology file,', NLAYS
            CALL M3WARN( PROGNAME, 0, 0, MESG )

            EMLAYS = NLAYS

        END IF

C.........  Abort if error found analyzing inputs
        IF ( EFLAG ) THEN
            MESG = 'Problem with inputs.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Update start date/time and duration from the environment
        CALL GETM3EPI( -1, SDATE, STIME, NSTEPS )
        TSTEP = 10000

C.........  Set up and open output file, which will primarily using I/O API 
C           settings from the cross-point met file (XNAME), which are 
C           already retrieved
        CALL OPENLAYOUT( SDATE, STIME, TSTEP, EMLAYS, REP_LAYR, XFLAG, 
     &                   IFDESC2, IFDESC3, METSCEN, CLOUDSHM, VGLVSXG, 
     &                   LNAME, RDEV )

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( 'POINT', ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Call subroutine to convert grid coordinates from lat-lon to
C           coordinate system of the destination grid
        CALL CONVRTXY( NSRC, GDTYP, P_ALP, P_BET, P_GAM, 
     &                 XCENT, YCENT, XLOCA, YLOCA )

C.........  Call elevated sources indicator file, even thought it might not
C           be opened - routine will initialize LMAJOR and LPING regardless
C           of whether the file is available.
        CALL RDPELV( PDEV, NSRC, .FALSE., NMAJOR, NPING )

C.........  If explicit plume rise, only explicit plume sources will be
C           output, but LMAJOR needs to be true for error checking.  So, set it
        IF( XFLAG ) LMAJOR = .TRUE.

C.........  Allocate memory for all remaining variables using dimensions 
C           obtained previously...

C.........  Allocate per-source arrays
        ALLOCATE( HFX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HFX', PROGNAME )
        ALLOCATE( HMIX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HMIX', PROGNAME )
        ALLOCATE( TSFC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSFC', PROGNAME )
        ALLOCATE( USTAR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'USTAR', PROGNAME )

C.........  Allocate per-source and per-layer arrays
        ALLOCATE( DDZH( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DDZH', PROGNAME )
        ALLOCATE( DDZF( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DDZF', PROGNAME )
        ALLOCATE( PRES( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )
        ALLOCATE( QV( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'QV', PROGNAME )
        ALLOCATE( TA( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TA', PROGNAME )
        ALLOCATE( UWIND( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UWIND', PROGNAME )
        ALLOCATE( VWIND( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VWIND', PROGNAME )
        ALLOCATE( ZF( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ZF', PROGNAME )
        ALLOCATE( ZH( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ZH', PROGNAME )
        ALLOCATE( ZSTK( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ZSTK', PROGNAME )

C.........  If hourly data input, allocate index array
        IF( HFLAG ) THEN
            ALLOCATE( LOCINDXH( NHRSRC,EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LOCINDXH', PROGNAME )
            LOCINDXH = 0   ! array
        END IF

C.........  Allocate layer fractions array: by source if not explicit, by
C           hour-specific source if it is explicit
        IF( XFLAG ) THEN
            ALLOCATE( LFRAC( NHRSRC,EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )
            ALLOCATE( TFRAC( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TFRAC', PROGNAME )

C.........  If computing plume rise...
        ELSE

C.............  Layer fractions for all sources
            ALLOCATE( LFRAC( NSRC,EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )

C.............  Allocate ungridding arrays
            ALLOCATE( ND( 4,NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ND', PROGNAME )
            ALLOCATE( NX( 4,NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NX', PROGNAME )
            ALLOCATE( CD( 4,NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CD', PROGNAME )
            ALLOCATE( CX( 4,NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CX', PROGNAME )

C.............  Allocate per-layer arrays from 1:EMLAYS
            ALLOCATE( WSPD( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WSPD', PROGNAME )
            ALLOCATE( DTHDZ( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTHDZ', PROGNAME )
            ALLOCATE( TFRAC( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TFRAC', PROGNAME )

C.............  Allocate per-layer arrays from 0:EMLAYS
            ALLOCATE( PRESF( 0:EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PRESF', PROGNAME )
            ALLOCATE( ZZF( 0:EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ZZF', PROGNAME )

C.............  Allocate array for tmp gridded, layered cross-point met data
            ALLOCATE( XBUF( METNGRID,NLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XBUF', PROGNAME )
            ALLOCATE( DBUF( NDOTS,NLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DBUF', PROGNAME )

C.............  Compute un-gridding matrices for dot and cross point met data
            CALL UNGRIDB( METNCOLS+1, METNROWS+1, 
     &                    XORIGDG, YORIGDG, XCELLDG, YCELLDG,
     &                    NSRC, XLOCA, YLOCA, ND, CD )

            CALL UNGRIDB( METNCOLS, METNROWS, METXORIG, METYORIG, 
     &                    XCELL, YCELL, NSRC, XLOCA, YLOCA, NX, CX )

C.............  Read time-independent ZF and ZH for non-hydrostatic Met data
C.............  Compute per-source heights
            IF( ZSTATIC ) THEN

                CALL RETRIEVE_IOAPI_HEADER( GNAME )
                CALL GET_VARIABLE_NAME( 'ZH', VNAME )
                CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D, 
     &                           STIME3D, XBUF )
                CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, ZH )

                CALL GET_VARIABLE_NAME( 'ZF', VNAME )
                CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D,
     &                           STIME3D, XBUF )
                CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, ZF )

C.................  Pre-process ZF and ZH to compute DDZH and DDZF
                CALL COMPUTE_DELTA_ZS

            END IF

        END IF     ! if explicit plume rise or not

C.........  Write out header to report, if any. This includes generating
C           format statement for the 
        IF( REP_LAYR .GT. 0 ) THEN

            MESG = 'Cy/St/Co, Plant'
            DO I = 1, NCHARS - 2
                L = LEN_TRIM( MESG )
                WRITE( MESG,'(A,I1,A)' ) MESG( 1:L ) // ', Char', I
            END DO
            L = LEN_TRIM( MESG )

            WRITE( RDEV,93040 ) REP_LAYR, MESG( 1:L )

        END IF

C.........  Set logical array for setting valid source characeristics columns
        LFG( 1:NCHARS ) = .TRUE.   ! array
        IF( NCHARS .LE. 8 ) LFG( NCHARS+1:9 ) = .FALSE.  ! array

C.........  Get variable names from surface meteorology file
        IF ( .NOT. XFLAG ) CALL RETRIEVE_IOAPI_HEADER( SNAME )

C.........  For each time step, compute the layer fractions...

        MESG = 'Calculating hourly layer fractions...'
        CALL M3MSG2( MESG )

        XBEG  = XORIG
        YBEG  = YORIG
        XEND  = XORIG + NCOLS * XCELL
        YEND  = YORIG + NROWS * YCELL
        LDATE = 0
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

            IF ( LDATE .NE. JDATE ) THEN
 
C.................  Write day and date message to stdout and log file
                CALL WRDAYMSG( JDATE, MESG )

C.................  Write day and date message to report file
                IF( RDEV .GT. 0 ) THEN
                    WRITE( RDEV,93000 ) MESG( 1:LEN_TRIM( MESG ) )
                END IF

                LDATE = JDATE
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

C.............  Write to report file if report feature is on
            IF( RDEV .GT. 0 ) THEN
                WRITE( RDEV,93020 ) HHMMSS( JTIME )
            END IF

C.............  Initialize layer fraction array
            LFRAC = 0.    ! 2-d array

C.............  If needed, read hourly plume rise and/or stack parameters...
C.............  Read source index
            IF ( HFLAG ) THEN

C.................  Do not give an error if could not read data, because it 
C                   might not be there
                IFLAG = .TRUE.
                IF( .NOT. READ3( HNAME, 'INDXH', ALLAYS3, 
     &                           JDATE, JTIME, LOCINDXH( 1,1 ) ) ) THEN
                    L1 = LEN_TRIM( HNAME )
                    WRITE( MESG,94010 ) 'WARNING: Could not read '//
     &                     '"IDXH" from file "' //
     &                     HNAME( 1:L1 ) // '", at', JDATE, ':', JTIME
                    CALL M3MESG( MESG )

                    LOCINDXH = 0       ! 2-d array
                    IFLAG = .FALSE. 

                END IF

C.................  Determine the number of valid hour-specific sources for
C                   the current hour
                DO I = NHRSRC, 1, -1
                    IF ( LOCINDXH( I,1 ) .NE. 0 ) EXIT
                END DO
                NHR = I

            END IF

C.............  Layer-1 fraction
            IF ( LFLAG .AND. IFLAG ) 
     &           CALL SAFE_READ3( HNAME, SPDATNAM(1), ALLAYS3, 
     &                            JDATE, JTIME, LAY1F )

C.............  Plume bottom and top
            IF ( BFLAG .AND. IFLAG ) THEN
                CALL SAFE_READ3( HNAME, SPDATNAM(2), ALLAYS3, 
     &                           JDATE, JTIME, PLMBOT )
                CALL SAFE_READ3( HNAME, SPDATNAM(3), ALLAYS3, 
     &                           JDATE, JTIME, PLMTOP )
            END IF

C.............  Temperatures
            IF ( TFLAG .AND. IFLAG ) 
     &           CALL SAFE_READ3( HNAME, SPDATNAM(4), ALLAYS3, 
     &                            JDATE, JTIME, HRSTKTK )

C.............  Velocity
            IF ( CFLAG .AND. IFLAG ) 
     &           CALL SAFE_READ3( HNAME, SPDATNAM(5), ALLAYS3, 
     &                            JDATE, JTIME, HRSTKVE )

C.............  Flow rate
            IF ( FFLAG .AND. IFLAG ) 
     &           CALL SAFE_READ3( HNAME, SPDATNAM(6), ALLAYS3, 
     &                            JDATE, JTIME, HRSTKFL )

C.............  Read time-dependent ZF and ZH for hydrostatic Met data
C.............  Compute per-source heights
            IF( .NOT. XFLAG .AND. .NOT. ZSTATIC ) THEN

                CALL SAFE_READ3( XNAME,'ZH',ALLAYS3,JDATE,JTIME,XBUF )
                CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, ZH )

                CALL SAFE_READ3( XNAME,'ZF',ALLAYS3,JDATE,JTIME,XBUF )
                CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, ZF )

C.................  Pre-process ZF and ZH to compute DDZH and DDZF
                CALL COMPUTE_DELTA_ZS

            END IF

C.............  Read and transform meteorology:
            IF ( .NOT. XFLAG ) THEN
            CALL SAFE_READ3( SNAME, 'HFX', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( METNGRID, NSRC, 1, NX, CX, XBUF, HFX )

            CALL SAFE_READ3( SNAME, 'PBL', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( METNGRID, NSRC, 1, NX, CX, XBUF, HMIX )

            CALL GET_VARIABLE_NAME( 'TGD', VNAME )
            CALL SAFE_READ3( SNAME, VNAME, ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( METNGRID, NSRC, 1, NX, CX, XBUF, TSFC )

            CALL SAFE_READ3( SNAME, 'USTAR', ALLAYS3, JDATE,JTIME,XBUF )
            CALL BMATVEC( METNGRID, NSRC, 1, NX, CX, XBUF, USTAR )

            CALL SAFE_READ3( XNAME, 'TA', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, TA )

            CALL SAFE_READ3( XNAME, 'QV', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, QV )

            CALL SAFE_READ3( XNAME, 'PRES', ALLAYS3, JDATE, JTIME,XBUF )
            CALL BMATVEC( METNGRID, NSRC, EMLAYS, NX, CX, XBUF, PRES )

            CALL SAFE_READ3( DNAME, 'UWIND', ALLAYS3, JDATE,JTIME,DBUF )
            CALL BMATVEC( NDOTS, NSRC, EMLAYS, ND, CD, DBUF, UWIND )

            CALL SAFE_READ3( DNAME, 'VWIND', ALLAYS3, JDATE,JTIME,DBUF )
            CALL BMATVEC( NDOTS, NSRC, EMLAYS, ND, CD, DBUF, VWIND )
            END IF

C.............  Precompute constants before starting source loop
            P  = ( SIGH(0) - VGLVSXG(0) ) / ( SIGH( 1 ) - SIGH( 0 ) )
            PP = 1.0 + P 
            QQ =     - P 

C.............  Loop through sources and compute plume rise
            K = 0
            DO S = 1, NSRC

                DM = STKDM( S )
                HT = STKHT( S )
                TK = STKTK( S )
                VE = STKVE( S )
                FL = 0.          ! initialize flow
	        XL = XLOCA( S )
	        YL = YLOCA( S )

C.................  Find source in index list if hourly data or used
                IF ( HFLAG ) THEN
                    K = FIND1( S, NHR, LOCINDXH( 1,1 ) )
                END IF

C.................  Skip source if explicit processing and source not on list
                IF ( XFLAG .AND. K .LE. 0 ) THEN
                    CYCLE

C.................  Skip source if it is outside output grid
                ELSE IF( XL .LT. XBEG .OR. XL .GT. XEND .OR.
     &              YL .LT. YBEG .OR. YL .GT. YEND     ) THEN
                    CYCLE

C.................  Skip source if it is minor source and assign layer fractions
C                   that put all emissions in layer 1
                ELSE IF( .NOT. LMAJOR( S ) ) THEN
                    IF( XFLAG ) THEN
                       WRITE( MESG,94010 ) 
     &                      'INTERNAL ERROR: LMAJOR(S) = FALSE for'//
     &                      'explicit plume source number', S
                       CALL M3MSG2( MESG )
                       CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                    ELSE
                        LFRAC( S,1 ) = 1.
                    END IF
                    CYCLE

                END IF

C.................  If hourly data available, check if source has hourly data
C                   for the current hour, then read hourly stack parameters
                IF ( HFLAG ) THEN

C.....................  If source has hourly data...
                    IF( K .GT. 0 ) THEN

C.........................  If hourly temperatures are available, reset
                        IF( TFLAG ) THEN
                            IF( HRSTKTK( K ) .GT. 0 ) TK = HRSTKTK( K )
                        END IF

C.........................  If hourly velocities are available, reset
                        IF( YFLAG ) THEN
                            IF( HRSTKVE( K ) .GT. 0 ) VE = HRSTKVE( K )
                        END IF

C.........................  If hourly flow rates are available, reset and
C                           recompute velocity if requested
                        IF( FFLAG ) THEN
                            IF( HRSTKFL( K ) .GT. 0 ) THEN
                                FL = HRSTKFL( K )
                                IF ( CFLAG ) THEN
                                    VE = FL / ( 0.25 * PI * DM * DM )
                                END IF  ! if velocity recalculation requested.
                            END IF      ! if flow valid for source
                        END IF          ! if flow in hourly file
                    END IF              ! if source is hourly
                END IF                  ! if hourly used

C.................  For explicit plume rise, assign plume top and
C                   plume bottom from hourly data, and setup weights for 
C                   allocation to the layers using static layer heights. Weight
C                   by penetration of plume into layer and the layer thickness.
C.................  This is the approach for UAM-style processing
                IF ( XFLAG ) THEN

C.....................  If plume bottom, and plume top are available, set 
C                       these and set to skip plume rise computation
                    IF( PLMBOT( K ) .GE. 0. .AND.
     &                  PLMTOP( K ) .GT. 0.       ) THEN
                        ZBOT = PLMBOT( K )
                        ZTOP = PLMTOP( K )

C.....................  Otherwise, set top and bottom of plume to be in layer 1.
                    ELSE
                        ZBOT = VGLVSXG( 1 ) * 0.5
                        ZTOP = ZBOT

                    END IF

                    SURFACE = 100.             ! percent to surface
                    LFULLHT( 0 ) = 0.
                    DO L = EMLAYS, 1, -1
                        LFULLHT( L ) = VGLVSXG( L )
                        LHALFHT( L ) = VGLVSXG( L-1 ) +
     &                                 0.5 * ( VGLVSXG(L)-VGLVSXG(L-1) )
                    
                        WEIGHTS( L ) = 100. * ( LFULLHT( EMLAYS ) - 
     &                                          LHALFHT( L )        ) /
     &                                          LFULLHT( EMLAYS )
                    END DO

C.................  For non-explicit plume rise, preprocess met data...
                ELSE

C.....................  Compute surface pressure (and convert to mb from Pa)
                    PSFC = 1.0E-2 * ( PP * PRES( 1,S ) + 
     &                                QQ * PRES( 2,S )   )

C.....................  Compute derived met vars needed before layer assignments
                    CALL PREPLM( EMLAYS, HMIX( S ), STKHT( S ), PSFC, 
     &                           TSFC( S ), DDZH( 1,S ), QV( 1,S ), 
     &                           TA( 1,S ), UWIND( 1,S ), VWIND( 1,S ), 
     &                           ZH( 1,S ), ZF( 1,S ), ZSTK( 1,S ), 
     &                           PRES( 1,S ), LSTK, LPBL, TSTK, WSTK,
     &                           DTHDZ, WSPD, ZZF )

C.....................  Trap USTAR at a minimum realistic value
                    USTMP = MAX( USTAR( S ), USTARMIN )

                    COMPUTE = .TRUE.

C.....................  If available, assign hourly plume top and plume bottom
                    IF ( BFLAG .AND. K .GT. 0 ) THEN

C.........................  If plume bottom, and plume top are available, set 
C                           these and set to skip plume rise computation
                        IF( PLMBOT( K ) .GE. 0. .AND.
     &                      PLMTOP( K ) .GT. 0.       ) THEN
                            ZBOT = PLMBOT( K )
                            ZTOP = PLMTOP( K )
                            COMPUTE = .FALSE.
                        END IF

                    END IF

C.....................  Compute plume rise for this source, if needed
                    IF ( COMPUTE ) THEN

                        CALL PLMRIS( EMLAYS, LPBL, LSTK, HFX(S), 
     &                           HMIX(S), DM, HT, TK, VE, TSTK, USTMP, 
     &                           DTHDZ, TA(1,S), WSPD, ZZF(0), ZH(1,S), 
     &                           ZSTK(1,S), WSTK, ZTOP, ZBOT )
                    END IF

C.....................  Setup for computing plume fractions, assuming uniform
C                       distribution in pressure (~mass concentration -- minor 
C                       hydrostatic assumption) from bottom to top.
                    SURFACE = PSFC
                    LFULLHT( 0:EMLAYS ) = ZZF ( 0:EMLAYS   )
                    LHALFHT( 1:EMLAYS ) = ZH  ( 1:EMLAYS,S )
                    WEIGHTS( 1:EMLAYS ) = PRES( 1:EMLAYS,S )

                END IF  ! if computing plume rise

C.................  Check plume rise for nonsense values
                IF( ZTOP .LT. STKHT( S ) .AND. K .LE. 0 ) THEN

                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 ) 
     &                     'WARNING: Top of plume found to be ' //
     &                     'less than top of stack for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )

                    CALL M3MESG( MESG )

                END IF

C.................  Allocate plume to layers
                CALL POSTPLM( EMLAYS, S, SURFACE, ZBOT, ZTOP, WEIGHTS, 
     &                        LFULLHT, LHALFHT, LTOP, TFRAC )

C.................  If hourly layer-1 fraction is present, reset this and re-
C                   normalize
C.................  Must account for the case where LAY1F value is 
C                   missing
                IF( LFLAG .AND. K .GT. 0 ) THEN
                    IF( LAY1F( K ) .GT. 0. .AND.
     &                  TFRAC( 1 ) .LT. 1.       ) THEN
                        TSUM = SUM( TFRAC( 2:EMLAYS ) )
                        TDIFF = TSUM + TFRAC( 1 ) - LAY1F( K )
                        FAC = TDIFF / TSUM

                        TFRAC( 1 ) = LAY1F( K )
                        TFRAC( 2:EMLAYS ) = TFRAC( 2:EMLAYS ) * FAC
                    END IF
                END IF

C.................  Check if layer fractions are negative and reset
C                   to output in the first layer if they are.
                X = MINVAL( TFRAC( 1:EMLAYS ) )
                IF( X .LT. 0 ) THEN

                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                    WRITE( MESG,94010 ) 
     &                     'WARNING: One or more negative plume ' //
     &                     'fractions found for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )//'.'//
     &                     CRLF() // BLANK10 // 'Plume reset to '//
     &                     'have all emissions in surface layer.'
                    CALL M3MESG( MESG )

                    TFRAC( 1 ) = 1.0
                    TFRAC( 2:EMLAYS ) = 0.0

                END IF

C.................  Store layer fractions
                IF( XFLAG ) THEN
                    LFRAC( K,1:EMLAYS ) = TFRAC( 1:EMLAYS )  ! array
                ELSE 
                    LFRAC( S,1:EMLAYS ) = TFRAC( 1:EMLAYS )  ! array

                END IF

C.................  Check if LTOP out of range, and report (will only work
C.................    if REP_LAYR env var has been set b/c default is -1
                IF( REP_LAYR .GT. 0 .AND. LTOP .GT. REP_LAYR ) THEN

                    CALL PARSCSRC( CSOURC( S ), NCHARS, SC_BEGP, 
     &                             SC_ENDP, LFG, I, CHARS )

                    WRITE( OUTFMT, 93042 ) PLTLEN3, NCHARS-2, CHRLEN3
                    WRITE( RDEV,OUTFMT ) S, IFIP( S ),
     &                   ( CHARS( I ), I = 2,NCHARS ), STKHT( S ), 
     &                     STKVE ( S ), STKTK( S ), TSTK, 
     &                     WSTK, LPBL, LTOP 

                END IF

            END DO    !  end loop on sources S

C.............  Write out layer fractions
            IF ( .NOT. WRITE3( LNAME, 'LFRAC', 
     &                         JDATE, JTIME, LFRAC ) ) THEN

                MESG = 'Problem writing "LFRAC" to file "' // 
     &                 LNAME( 1:LEN_TRIM( LNAME ) ) // '."'
     
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

C.............  For explicit plume rise, also write out source numbers
            IF ( XFLAG ) THEN
                IF ( .NOT. WRITE3( LNAME, 'INDXH', JDATE, 
     &                             JTIME, LOCINDXH( 1,1 ) ) ) THEN

                    MESG = 'Problem writing "LFRAC" to file "' // 
     &                     LNAME( 1:LEN_TRIM( LNAME ) ) // '."'

                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF
            END IF


            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO     !  end loop on time steps T


C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( 8X, 'at time ', A8 )

93040   FORMAT( 'Sources with top of plume greater than layer', I3, //,
     &          'Src ID, ', A, ', H[m], ', 'V[m/s], ', 'Ts[K], ', 
     &          'Ta[K], ', 'U[m/s], ', 'LPBL, ', 'LTOP' )

93042   FORMAT( '( I6, ",", I6.6, ",", A', I2.2, ', ","', I2.2, '(A', 
     &          I2.2, ',", ") , F6.1, ", ", F6.2, ", ", F6.1, ", ",', 
     &          'F5.1, ", ", F6.2, ", ", I3, ", ", I3 )' )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 12( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS

C.............  This internal subprogram tries to retrieve the I/O API header
C               and aborts if it was not successful
            SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C----------------------------------------------------------------------

            IF ( .NOT. DESC3( FILNAM ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            END SUBROUTINE RETRIEVE_IOAPI_HEADER

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram resolves the differences in 
C               variable names for different version of the Met files
            SUBROUTINE GET_VARIABLE_NAME( INNAME, OUTNAME )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: INNAME    ! variable name to check
            CHARACTER(*), INTENT(OUT) :: OUTNAME   ! variable name to read

C.............  Local variables
            INTEGER J

C----------------------------------------------------------------------

C.............  Search for variable name in the list of names
            J = INDEX1( INNAME, NVARS3D, VNAME3D )

C.............  If the input name is there, then set output name and return
            IF( J .GT. 0 ) THEN
                OUTNAME = INNAME
                RETURN
            END IF

C.............  Set output name
C.............  Currently there is only one alternative for each
            SELECT CASE( INNAME )
            CASE( 'ZH' )
                OUTNAME = 'X3HT0M'
            CASE( 'ZF' )
                OUTNAME = 'X3HT0F'
            CASE( 'TGD' ) 
                OUTNAME = 'TEMPG'
            CASE DEFAULT
                MESG = 'INTERNAL ERROR: Do not have an alternative ' //
     &                 'name for met variable ' // INNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END SELECT

            RETURN

            END SUBROUTINE GET_VARIABLE_NAME

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram tries to read a variable from an
C               I/O API file, and aborts if not successful.
            SUBROUTINE SAFE_READ3( FILNAM, VARNAM, LAYER, 
     &                             JDATE, JTIME, XBUF     )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM    ! logical file name
            CHARACTER(*) VARNAM    ! variable name
            INTEGER      LAYER     ! layer number (or ALLAYS3)
            INTEGER      JDATE     ! Julian date
            INTEGER      JTIME     ! time
            REAL         XBUF( * ) ! read buffer

C----------------------------------------------------------------------

            IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,
     &                        JDATE, JTIME, XBUF ) ) THEN

                L1 = LEN_TRIM( VARNAM )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not read "' // VARNAM( 1:L1 ) // 
     &                 '" from file "' // FILNAM( 1:L2 ) // '."'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

            END SUBROUTINE SAFE_READ3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram computes DDZH and DDZF
            SUBROUTINE COMPUTE_DELTA_ZS

C----------------------------------------------------------------------

            DO S = 1, NSRC

                ZZ0 = ZF( 1,S )
                ZSTK ( 1,S ) = ZZ0 - STKHT( S )
                ZF0 = ZF( 1,S )
                DDZF( 1,S ) = 1.0 / ZF0

                DO L = 2, EMLAYS

                    ZZ1 = ZF( L,S )
                    ZSTK( L  ,S ) = ZZ1 - STKHT( S )
                    DDZH( L-1,S ) = 1.0 / ( ZZ1 - ZZ0 )
                    ZZ0 = ZZ1
                    ZF1 = ZF( L,S )
                    DDZF( L,S ) = 1.0 / ( ZF1 - ZF0 )
                    ZF0 = ZF1

                END DO

            END DO  ! End processing for intermediate layer height calcs

            END SUBROUTINE COMPUTE_DELTA_ZS

        END PROGRAM LAYPOINT

