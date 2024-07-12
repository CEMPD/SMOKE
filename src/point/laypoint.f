
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
C     Updated Feb. 2005 with changes from J. Godowitch & G. Pouliot
C     Updated to read either ACRESBURNED(EPA) or AREA(Bluesky) as acres burned variables
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: XLOCA, YLOCA, STKDM, STKHT, STKTK, STKVE,
     &                      CSOURC, CIFIP, CSCC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NHRSRC, HRSTKTK, HRSTKVE, HRSTKFL, LMAJOR,
     &                     LAY1F, PLMBOT, PLMTOP

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, CATLEN, SC_BEGP, SC_ENDP,
     &                     NSRC, NCHARS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NROWS, NCOLS, NLAYS, VGTYP, VGTOP,
     &                     XORIG, YORIG, COORD, GDTYP, P_ALP, P_BET,
     &                     P_GAM, XCENT, YCENT, XCELL, YCELL, GRDNM

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        USE MODGRDLIB

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKMETEM
        CHARACTER(2)    CRLF
        LOGICAL         DSCM3GRD
        LOGICAL         DSCM3LAY
        INTEGER         ENVINT
        REAL            ENVREAL
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         MFIND1
        CHARACTER(50)   GETCFDSC
        CHARACTER(10)   HHMMSS
        INTEGER         INDEX1
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         WKDAY
        REAL            STR2REAL

        EXTERNAL   CHKMETEM, CRLF, DSCM3GRD, DSCM3LAY, ENVINT, ENVREAL,
     &             ENVYN, FIND1, GETCFDSC, HHMMSS, INDEX1, MMDDYY,
     &             PROMPTFFILE, PROMPTMFILE, VERCHAR, WKDAY, STR2REAL

C...........  LOCAL PARAMETERS and their descriptions:

        REAL, PARAMETER :: USTARMIN  = 0.1     ! Min valid value for USTAR
        REAL, PARAMETER :: CONVPA    = 1.0E-2  ! conversion factor for Pa to mb
        REAL, PARAMETER :: ZERO      = 0.0     ! dummy zero value
        REAL, PARAMETER :: BTU2M4PS3 = 0.00000258  ! conv. factor for bouyancy flux

        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

C.........  Indicator for which public inventory arrays need to be read
        INTEGER,            PARAMETER :: NINVARR = 9
        CHARACTER(IOVLEN3), PARAMETER :: IVARNAMS( NINVARR ) =
     &                                 ( / 'CIFIP          '
     &                                   , 'XLOCA          '
     &                                   , 'YLOCA          '
     &                                   , 'STKHT          '
     &                                   , 'STKDM          '
     &                                   , 'STKTK          '
     &                                   , 'STKVE          '
     &                                   , 'CSCC           '
     &                                   , 'CSOURC         ' / )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point source stack parameters:

C.........  Allocatable, per-source meteorology variables
        REAL   , ALLOCATABLE :: HFX  ( : )    !  sensible heat flux (watts/m^2)
        REAL   , ALLOCATABLE :: HMIX ( : )    !  mixing height (m)
        REAL   , ALLOCATABLE :: TSFC ( : )    !  surface temperature (deg K)
        REAL   , ALLOCATABLE :: USTAR( : )    !  friction velocity (m/s)
        REAL   , ALLOCATABLE :: PRSFC( : )    !  surface pressure (Pascals)

C.........  Allocatable, per-source fire data variables
        REAL   , ALLOCATABLE :: BFLX ( : )    !  Briggs bouyancy flux (m^4/s^3)
        REAL   , ALLOCATABLE :: ACRES( : )    !  area burned (acres)

C.........  Allocatable, per-source and per layer meteorology variables.
C           Dimensioned by layers, then sources
        REAL   , ALLOCATABLE :: DDZH ( :,: )  !  1/( zh(l) - zh(l-1) )
        REAL   , ALLOCATABLE :: DDZF ( :,: )  !  1/( zf(l) - zf(l-1) )
        REAL   , ALLOCATABLE :: PRES ( :,: )  !  pressure (Pa)
        REAL   , ALLOCATABLE :: QV   ( :,: )  !  mixing ratio (kg/kg)
        REAL   , ALLOCATABLE :: TA   ( :,: )  !  temperature (K)
        REAL   , ALLOCATABLE :: UWIND( :,: )  !  wind speed (m/s)
        REAL   , ALLOCATABLE :: VWIND( :,: )  !  wind speed (m/s)
        REAL   , ALLOCATABLE :: ZF   ( :,: )  !  layer surface height (m)
        REAL   , ALLOCATABLE :: ZH   ( :,: )  !  layer center  height (m)
        REAL   , ALLOCATABLE :: ZSTK ( :,: )  !  zf( l,s ) - stkht(s) (m)
        REAL   , ALLOCATABLE :: DENS ( :,: )  !  air density (kg/m^3)

C........... added variables by George Pouliot 2/8/07 (for bug fix wildfires)
        INTEGER :: DAY_NSRC  !number of sources in daily file (not the same as in inventory)
        REAL, ALLOCATABLE    ::  DAY_ACRES(:)   ! acres in daily file by source
        INTEGER, ALLOCATABLE ::  DAY_INDEX(:)   ! index of sources in daily file
        INTEGER :: MY_LOOP, MY_INDEX   ! for linear search


C.........  Allocatable, temporary per-layer variables from 1:EMLAYS
        REAL   , ALLOCATABLE :: WSPD ( : )    !  wind speed (m/s)
        REAL   , ALLOCATABLE :: DTHDZ( : )    !  gradient of THETV
        REAL   , ALLOCATABLE :: TFRAC( : )    !  Temporary LFRAC
        REAL   , ALLOCATABLE :: PFRAC( : )    !  Previous LFRAC

C.........  Allocatable, temporary per-layer variables from 0:EMLAYS
        REAL   , ALLOCATABLE :: PRESF( : )    !  pressure at full-levels
        REAL   , ALLOCATABLE :: ZZF  ( : )    !  elevation at full-levels

C.........  Allocatable cross- OR dot-point meteorology input buffers
        REAL   , ALLOCATABLE :: XBUF( :,: )   ! cross-point
        REAL   , ALLOCATABLE :: DBUF( :,: )   ! dot-point

C.........  Allocatable cross-point surface grid coordinates
        REAL   , ALLOCATABLE :: XXVALS( :,: )   ! x values
        REAL   , ALLOCATABLE :: XYVALS( :,: )   ! y values

C.........  Allocatable dot-point surface grid coordinates
        REAL   , ALLOCATABLE :: DXVALS( :,: )   ! x values
        REAL   , ALLOCATABLE :: DYVALS( :,: )   ! y values

C.........  Allocatable un-gridding matrices
        INTEGER, ALLOCATABLE :: ND( : )     !  records per cell (unused) 
        INTEGER, ALLOCATABLE :: NX( : )     !  records per cell (unused)
        INTEGER, ALLOCATABLE :: GD( : )     !  dot-point, cell index 
        INTEGER, ALLOCATABLE :: GDS( : )    !  dot-point, cell index shift
        INTEGER, ALLOCATABLE :: GX( : )     !  cross-point, cell index
        INTEGER, ALLOCATABLE :: SN( : )        ! record index (unused)
        INTEGER, ALLOCATABLE :: INDX( : )      ! sorting index (unused)

C.........  Output layer fractions, dimensioned NSRC, emlays
        REAL   , ALLOCATABLE :: LFRAC( :, : )

C.........  Input/output hour-specific data index, dimensioned by NSRC and
C           by EMLAYS, so that index can be written to PLAY_EX file
        INTEGER, ALLOCATABLE :: LOCINDXH( :,: )

C.........  Fixed-dimension arrays
        REAL         LFULLHT( 0:MXLAYS3 )     !  full-level heights [m]
        REAL         LHALFHT( 1:MXLAYS3 )     !  half-level heights [m]
        REAL         TEMPS  ( 1:MXLAYS3 )     !  half-level temps (K)
        REAL         SIGH   ( 0:MXLAYS3-1 )   !  half-level sigma values
        REAL         VGLVSXG( 0:MXLAYS3 )     !  vertical coord values
        REAL         WEIGHTS( 1:MXLAYS3 )     !  tmp weights for vertical aloc

C...........   Logical names and unit numbers

        INTEGER         IDEV    !  tmp unit number if ENAME is map file
        INTEGER         LDEV    !  log file
        INTEGER      :: PDEV = 0!  elevated/PinG source file
        INTEGER      :: RDEV = 0!  optional report iff REP_LAYER_MAX is set
        INTEGER         SDEV    !  ASCII part of inventory file

        CHARACTER(16)   ANAME   !  ASCII point-source inventory file
        CHARACTER(16)   CNAME   !  cross-point surface grid file
        CHARACTER(16)   DNAME   !  dot-point layered met file name
        CHARACTER(16)   ENAME   !  point-source inventory input file
        CHARACTER(16)   GNAME   !  cross-point layered grid file name
        CHARACTER(16)   HNAME   !  hourly input file name
        CHARACTER(16)   INAME   !  tmp name for inven file of unknown fmt
        CHARACTER(16)   LNAME   !  layer fractions matrix output file
        CHARACTER(16)   SNAME   !  cross-point surface met file name
        CHARACTER(16)   TGNAM   !  Ground temperature variable name
        CHARACTER(16)   TNAME   !  dot-point surface grid file
        CHARACTER(16)   XNAME   !  cross-point layered met file name
        CHARACTER(16)   MNAME   !  temporalized data file name
        CHARACTER(16) DAYNAME   !  daily inventory file name

C...........   Other local variables
        INTEGER          I, J, K, L, L1, L2, S, T   ! counters and indices

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
        INTEGER          IPVERT    ! plume vertical spread method
        INTEGER          NEXCLD    ! number of sources excluded on grid
        INTEGER          COL       ! column number
        INTEGER          ROW       ! row number
        INTEGER          NCEL      ! cell index number

        REAL             X, Y, P, Q
        REAL             DM, HT, TK, VE, FL  ! tmp stack parameters
        REAL*8           XBEG, XEND, XL  ! tmp x-coords
        REAL*8           YBEG, YEND, YL  ! tmp y-coords
        REAL             FAC       !  tmp factor for renormalizing
        REAL             PSFC      !  surface pressure (Pa)
        REAL             PBL       !  PBL height (m)
        REAL             SURFACE   !  tmp weight at surface
        REAL             TDIFF     !  tmp layer frac diff for renormalizing
        REAL             TSTK      !  temperature at top of stack (K)
        REAL             TSUM      !  tmp layer frac sum for renormalizing
        REAL             WSTK      !  wind speed  at top of stack (m/s)
        REAL             ZZ0, ZZ1, ZF0, ZF1
        REAL             ZBOT,TBOT !  plume bottom elevation (m)
        REAL             ZTOP,TTOP !  plume top    elevation (m)
        REAL             ZPLM      !  plume centerline height above stack (m)
        REAL             USTMP     !  tmp storage for ustar (m/s)
        REAL             TMPBFLX   !  tmp Briggs bouyancy (m^4/s^4)
        REAL             TMPACRE   !  tmp area burned (acres)
        REAL             ALTITUDE  !  aircraft altitude (m)

        REAL(8)          METXORIG  ! cross grid X-coord origin of met grid
        REAL(8)          METYORIG  ! cross grid Y-coord origin of met grid
        REAL(8)          XCELLDG   ! dot grid X-coordinate cell dimension
        REAL(8)          YCELLDG   ! dot grid Y-coordinate cell dimension
        REAL(8)          XORIGDG   ! dot grid X-coordinate origin of grid
        REAL(8)          YORIGDG   ! dot grid Y-coordinate origin of grid

        LOGICAL       :: BFLAG = .FALSE.  ! true: use plume bottom and top
        LOGICAL       :: CFLAG = .FALSE.  ! true: recalc vel w/ flow & diam
        LOGICAL     :: COMPUTE = .FALSE.  ! true: compute plume rise
        LOGICAL       :: EFLAG = .FALSE.  ! error flag
        LOGICAL       :: FFLAG = .FALSE.  ! true: use hourly flow rate
        LOGICAL       :: GFLAG = .FALSE.  ! true: use variable grid
        LOGICAL       :: HFLAG = .FALSE.  ! true: hourly input used
        LOGICAL       :: IFLAG = .FALSE.  ! true: hr data okay for timestep
        LOGICAL       :: LFLAG = .FALSE.  ! true: use hourly layer 1 fraction
        LOGICAL       :: PFLAG = .FALSE.  ! true: compute plm ris for iteration
        LOGICAL       :: TFLAG = .FALSE.  ! true: use hourly temperatures
        LOGICAL       :: VFLAG = .FALSE.  ! true: use elevated file (PELV)
        LOGICAL       :: XFLAG = .FALSE.  ! true: process ONLY explicit sources
        LOGICAL       :: YFLAG = .FALSE.  ! true: use hourly velocities
        LOGICAL       :: ZSTATIC = .TRUE. ! true: Get heights from GRID_CRO file
        LOGICAL          LFG( 9 )         ! true: source characteristic is valid
        LOGICAL       :: AIRFLAG = .FALSE.! true: calculate plumes for aircraft/EDMS data
        LOGICAL       :: FIREFLAG =.FALSE.! true: calculate plumes for fire data
        LOGICAL       :: FPB1FLAG =.FALSE.! true: set fire plume bottom to layer 1
        LOGICAL       :: HOURFIRE =.FALSE.! true: use hourly fire data
        LOGICAL, SAVE :: WILDFLAG =.TRUE. ! true: calculate plumes for fire data (NO BLUESKY USE)

        CHARACTER(50)    CHARS( 9 )!  tmp source characeristics
        CHARACTER(50) :: METSCEN   !  temporary string for met scenario name
        CHARACTER(50) :: CLOUDSHM  !  temporary string for cloud scheme name
        CHARACTER(80) :: GDESC     !  grid description
        CHARACTER(256)   OUTFMT    !  output format for RDEV report
        CHARACTER(256)   BUFFER    !  source characteristics buffer
        CHARACTER(512)   MESG      !  buffer for M3EXIT() messages

        CHARACTER(IOVLEN3) VNAME      ! variable name buffer
        CHARACTER(IOVLEN3) COORD3D    ! coordinate system name
        CHARACTER(IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER(16) :: PROGNAME = 'LAYPOINT'   !  program name

C***********************************************************************
C   begin body of program LAYPOINT

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
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

        MESG = 'Vertical spread method'
        IPVERT = ENVINT( 'VERTICAL_SPREAD', MESG, 0, IOS )

        MESG = 'Ground temperature'
        CALL ENVSTR( 'PLUME_GTEMP_NAME', MESG, 'TEMP2', TGNAM, IOS )

        MESG = 'Set the fire plume bottom to layer 1'
        FPB1FLAG = ENVYN( 'FIRE_BOTTOM_LAYER_1_YN', MESG, .TRUE., IOS )

        MESG = 'Use aircraft inventory data'
        AIRFLAG = ENVYN( 'USE_EDMS_DATA_YN', MESG, .FALSE., IOS )

C.........  Reset flags for aircraft emissions
        IF( AIRFLAG ) THEN
            VFLAG = .FALSE.    ! no SMK_SPECELEV_YN
            HFLAG = .FALSE.    ! no HOUR_PLUMDATA_YN
            XFLAG = .FALSE.    ! no EXPLICIT_PLUMES_YN
            MESG = 'NOTE: To process EDMS aircraft inventory, ' //
     &         'SMK_SPECELEV_YN, HOUR_PLUMEDATA_YN, and  ' //
     &          CRLF() //BLANK10 // 'EXPLICIT_PLUMES_YN were set to N'
            CALL M3MSG2( MESG )
        END IF

        MESG = 'Use fire plume rise calculations'
        FIREFLAG = ENVYN( 'FIRE_PLUME_YN', MESG, .FALSE., IOS )

        IF( FIREFLAG ) THEN
            MESG = 'Use hourly fire data'
            HOURFIRE = ENVYN( 'HOURLY_FIRE_YN', MESG, .FALSE., IOS )
        END IF

C.........  If processing fire data without hourly data, need area
C           and heat flux values
        IF( FIREFLAG .AND. .NOT. HOURFIRE ) THEN
            MESG = 'Hourly heat flux value (BTU/hr)'
            TMPBFLX = ENVREAL( 'FIRE_HFLUX', MESG, 0., IOS )

C.............  Convert BTU/hr to Briggs bouyancy
            TMPBFLX = TMPBFLX * BTU2M4PS3  ! array

            MESG = 'Daily area burned (acres/day)'
            TMPACRE = ENVREAL( 'FIRE_AREA', MESG, 0., IOS )
        END IF

        MESG = 'Use variable grid definition'
        GFLAG = ENVYN( 'USE_VARIABLE_GRID', MESG, .FALSE., IOS )

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

C.........  Prompt for and open inventory file
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

C.........  Store source-category-specific header information,
C           including the inventory pollutants in the file (if any).
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

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
                MESG = 'NOTE: Layer-1 fraction hourly input ' //
     &                 'will be used '//CRLF()// BLANK10//
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
                MESG = 'NOTE: Plume top and bottom in hourly ' //
     &                 'input will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( PLMBOT( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PLMBOT', PROGNAME )
                ALLOCATE( PLMTOP( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PLMTOP', PROGNAME )

            END IF

C.............  Check for temperatures
            I = INDEX1( SPDATNAM( 4 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                TFLAG = .TRUE.
                MESG = 'NOTE: Temperatures hourly input ' //
     &                 'will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( HRSTKTK( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'HRSTKTK', PROGNAME )

            END IF

C.............  Check for velocity
            I = INDEX1( SPDATNAM( 5 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                YFLAG = .TRUE.
                MESG = 'NOTE: Velocities hourly input ' //
     &                 'will be used '//CRLF()// BLANK10//
     &                 'to allocate plumes for some sources.'
                CALL M3MSG2( MESG )

                ALLOCATE( HRSTKVE( NHRSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'HRSTKVE', PROGNAME )

            END IF

C.............  Check for flow rate
            I = INDEX1( SPDATNAM( 6 ), NVARS3D, VNAME3D )
            IF ( I .GT. 0 ) THEN
                FFLAG = .TRUE.
                MESG = 'NOTE: Flow rate hourly input ' //
     &                 'will be used '//CRLF()// BLANK10//
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

C.........  If using hourly fire data, open hourly and daily files
C           Will read heat flux from PTMP and area burned from PDAY
        IF( HOURFIRE ) THEN
            MNAME = PROMPTSET(
     &          'Enter logical name for the HOURLY EMISSIONS file',
     &          FSREAD3, CRL // 'TMP', PROGNAME )

C.............  Check to see if appropriate variable list exists
            CALL RETRIEVE_SET_HEADER( MNAME )

            I = INDEX1( 'HFLUX', NVARSET, VNAMESET )
            IF( I <= 0 ) THEN
                MESG = 'ERROR: Cannot find heat flux variable ' //
     &                 '"HFLUX" in hourly emissions file.'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            DAYNAME = PROMPTMFILE(
     &               'Enter logical name for DAY-SPECIFIC file',
     &               FSREAD3, CRL // 'DAY', PROGNAME )

C.............  Check to see if appropriate variable list exists
            CALL RETRIEVE_IOAPI_HEADER( DAYNAME )

            I = INDEX1( 'ACRESBURNED', NVARS3D, VNAME3D )
            J = INDEX1( 'AREA', NVARS3D, VNAME3D )

            IF( I <= 0 .OR. J > 0 ) WILDFLAG = .FALSE.

            IF( I <= 0 .AND. J <= 0 .AND. .NOT. WILDFLAG ) THEN
                MESG = 'ERROR: Cannot find acres burned ' //
     &                 'variable "ACRESBURNED" or "AREA" in daily ' //
     &                  CRLF() // BLANK10 // 'inventory file '
     &                  // TRIM( DAYNAME )
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            DAY_NSRC = NROWS3D
            WRITE( MESG,94010 )'NOTE: Number of Sources in Daily File',
     &                         DAY_NSRC
            CALL M3MSG2( MESG )

            IF( EFLAG ) THEN
                MESG = 'Problem with hourly fire data inputs'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

C.........  If not explicit plume rise only, open and process other met files
        IF ( .NOT. XFLAG ) THEN

            SNAME = PROMPTMFILE(
     &          'Enter name for CROSS-POINT SURFACE MET file',
     &          FSREAD3, 'MET_CRO_2D', PROGNAME )

            XNAME = PROMPTMFILE(
     &          'Enter name for CROSS-POINT LAYERED MET file',
     &          FSREAD3, 'MET_CRO_3D', PROGNAME )

            DNAME = PROMPTMFILE(
     &          'Enter name for DOT-POINT LAYERED MET file',
     &          FSREAD3, 'MET_DOT_3D', PROGNAME )

            IF( GFLAG ) THEN
                CNAME = PROMPTMFILE(
     &              'Enter name for CROSS-POINT SURFACE GRID file',
     &              FSREAD3, 'GRID_CRO_2D', PROGNAME )

                TNAME = PROMPTMFILE(
     &              'Enter name for DOT-POINT SURFACE GRID file',
     &              FSREAD3, 'GRID_DOT_2D', PROGNAME )
            END IF

C.............  Check multiple met files for consistency
            IF( GFLAG ) THEN
                EFLAG = ( .NOT. CHKMETEM( CNAME, SNAME, TNAME,
     &                                    'NONE', XNAME, DNAME ) )
            ELSE
                EFLAG = ( .NOT. CHKMETEM( 'NONE', SNAME, 'NONE',
     &                                    'NONE', XNAME, DNAME ) )
            END IF

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

            IF( VGTYP == -9999 ) VGTYP = 7   ! Reset WRF hybrid to WRF sigma layers

            SELECT CASE( VGTYP )
            CASE ( VGSGPH3, VGHVAL3, VGWRFEM )
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
C           If using a variable grid, do not allow subgrids
        IF( GFLAG ) THEN
            CALL CHKGRID( 'GRIDDESC', 'GRIDDESC', 0, EFLAG )

            IF( EFLAG ) THEN
                MESG = 'Problem with variable grid input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        ELSE
            CALL CHKGRID( 'GRIDDESC', 'GRIDDESC' , 1, EFLAG )

            IF( EFLAG ) THEN
                MESG = 'Problem with gridded input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

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
        CALL GETM3EPI( -1, SDATE, STIME, TSTEP, NSTEPS )
        TSTEP = 10000   ! Only 1-hour time step supported

C.........  Set up and open output file, which will primarily using I/O API
C           settings from the cross-point met file (XNAME), which are
C           already retrieved
        CALL OPENLAYOUT( SDATE, STIME, TSTEP, EMLAYS, REP_LAYR, XFLAG,
     &                   IFDESC2, IFDESC3, METSCEN, CLOUDSHM, VGLVSXG,
     &                   GFLAG, GRDNM, LNAME, RDEV )

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( 'POINT', ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  For fire data, set stack height to zero regardless of inventory
C           Inventory values may have been "corrected" by Smkinven
        IF( FIREFLAG ) THEN
            STKHT = 0.
        END IF

C.........  Call subroutine to convert grid coordinates from lat-lon to
C           coordinate system of the destination grid
        CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
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
        ALLOCATE( PRSFC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRSFC', PROGNAME )

        IF( HOURFIRE ) THEN
            ALLOCATE( BFLX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BFLX', PROGNAME )
            ALLOCATE( DAY_ACRES( DAY_NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DAY_ACRES', PROGNAME )
            ALLOCATE( DAY_INDEX( DAY_NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DAY_INDEX', PROGNAME )
            ALLOCATE( ACRES( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ACRES', PROGNAME )
        END IF

C.........  Allocate per-source and per-layer arrays
        ALLOCATE( DDZH( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DDZH', PROGNAME )
        ALLOCATE( DDZF( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DDZF', PROGNAME )
        ALLOCATE( PRES( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )
        ALLOCATE( DENS( EMLAYS,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DENS', PROGNAME )
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
            LFRAC = 0.0
            TFRAC = 0.0

C.........  If computing plume rise...
        ELSE

C.............  Layer fractions for all sources
            ALLOCATE( LFRAC( NSRC,EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )
            LFRAC = 0.0

C.............  Allocate ungridding arrays
            ALLOCATE( GD( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GD', PROGNAME )
            ALLOCATE( GDS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GDS', PROGNAME )
            ALLOCATE( GX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GX', PROGNAME )

C.............  Allocate per-layer arrays from 1:EMLAYS
            ALLOCATE( WSPD( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WSPD', PROGNAME )
            ALLOCATE( DTHDZ( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DTHDZ', PROGNAME )
            ALLOCATE( TFRAC( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TFRAC', PROGNAME )
            ALLOCATE( PFRAC( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PFRAC', PROGNAME )
            TFRAC = 0.0
            PFRAC = 0.0

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

C.........  Allocate memory so that we can use the GENPTCEL
            ALLOCATE( NX( METNGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NX', PROGNAME )
            ALLOCATE( ND( NDOTS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ND', PROGNAME )
            ALLOCATE( INDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )
            ALLOCATE( SN( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SN', PROGNAME )

C.............  Compute un-gridding matrices for dot and cross point met data
            IF( GFLAG ) THEN
C.............  If variable grid, allocate latitude and longitude arrays
                ALLOCATE( XXVALS( METNCOLS,METNROWS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'XXVALS', PROGNAME )
                ALLOCATE( XYVALS( METNCOLS,METNROWS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'XYVALS', PROGNAME )
                ALLOCATE( DXVALS( METNCOLS+1,METNROWS+1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DXVALS', PROGNAME )
                ALLOCATE( DYVALS( METNCOLS+1,METNROWS+1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DYVALS', PROGNAME )

C.............  Determine grid cells for these coordinate locations
C.............  Determine the X & Y values from the lat and lon
                CALL SAFE_READ3( TNAME, 'LON', 1, 0, 0, DXVALS )
                CALL SAFE_READ3( TNAME, 'LAT', 1, 0, 0, DYVALS )

                CALL CONVRTXY( (METNCOLS+1), (METNROWS+1), GDTYP, GRDNM,
     &                         P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                         DXVALS, DYVALS )

                CALL SAFE_READ3( CNAME, 'LON', 1, 0, 0, XXVALS )
                CALL SAFE_READ3( CNAME, 'LAT', 1, 0, 0, XYVALS )

                CALL CONVRTXY( METNCOLS, METNROWS, GDTYP, GRDNM,
     &                         P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                         XXVALS, XYVALS )

                CALL GENPTVCEL( NSRC, NDOTS, DXVALS, DYVALS, NEXCLD, ND,
     &                      INDX, GD, SN )
                CALL GENPTVCEL( NSRC, METNGRID, XXVALS, XYVALS, NEXCLD, NX,
     &                      INDX, GX, SN )
            ELSE
                CALL GENPTCEL( NSRC, METNGRID, XLOCA, YLOCA, NEXCLD, NX,
     &                      INDX, GX, SN )
C.............  Calculate grid dot reference cell index
                DO S = 1, NSRC
                    NCEL = GX( S )
                    ROW = NCEL / NCOLS
                    IF( MOD( NCEL, NCOLS ) .GT. 0 ) ROW = ROW + 1
                    COL = NCEL - ( ROW - 1 ) * NCOLS
                    GD( S ) = COL + ( NCOLS + 1 ) * ( ROW - 1 )
                END DO

            END IF

C.............  Read time-independent ZF and ZH for non-hydrostatic Met data
C.............  Compute per-source heights
            IF( ZSTATIC ) THEN

C.................  Open GRIDCRO3D file and check that it matches other met files
                GNAME = PROMPTMFILE(
     &              'Enter name for CROSS-POINT LAYERED GRID file',
     &              FSREAD3, 'GRID_CRO_3D', PROGNAME )

                IF( GFLAG ) THEN
                    EFLAG = ( .NOT. CHKMETEM( CNAME, SNAME, TNAME,
     &                                        GNAME, XNAME, DNAME ) )
                ELSE
                    EFLAG = ( .NOT. CHKMETEM( 'NONE', SNAME, 'NONE',
     &                                         GNAME, XNAME, DNAME ) )
                END IF

                IF( EFLAG ) THEN
                    MESG = 'GRID_CRO_3D file inconsistent with ' //
     &                     'other met files.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                CALL RETRIEVE_IOAPI_HEADER( GNAME )
                CALL GET_VARIABLE_NAME( 'ZH', VNAME )
                CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D,
     &                           STIME3D, XBUF )
                CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, ZH )

                CALL GET_VARIABLE_NAME( 'ZF', VNAME )
                CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D,
     &                           STIME3D, XBUF )
                CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, ZF )

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
     &                     '"INDXH" from file "' //
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

C.............  If using hourly fire data, read in acres burned and
C               heat flux
            IF( HOURFIRE ) THEN

C.................  Can't use SAFE_READ3 because data is stored in
C                   a fileset
                IF( .NOT. READSET( MNAME, 'HFLUX', ALLAYS3, ALLFILES,
     &                             JDATE, JTIME, BFLX ) ) THEN     ! PTMP file (all sources)

                    MESG = 'Could not read "HFLUX" from file "' //
     &                     TRIM( MNAME ) // '".'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Convert BTU/hr to Briggs bouyancy
                BFLX = BFLX * BTU2M4PS3  ! array

C.................  Day-specific file replicates data at each hour
C                   of the day
                IF( WILDFLAG ) THEN

                    CALL SAFE_READ3( DAYNAME, 'ACRESBURNED', ALLAYS3,
     &                           JDATE, JTIME, DAY_ACRES )   ! Wildfire inventory format

                     IF ( .NOT. READ3( DAYNAME, 'INDXD', ALLAYS3,
     &                        JDATE, JTIME, DAY_INDEX ) ) THEN

                           MESG = 'Could not read "INDXD" from file "'//
     &                     TRIM( DAYNAME ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                     END IF
                ELSE
                    CALL SAFE_READ3( DAYNAME, 'AREA', ALLAYS3,
     &                           JDATE, JTIME, DAY_ACRES )   ! Bluesky2inv format

                    IF ( .NOT. READ3( DAYNAME, 'INDXD', ALLAYS3,
     &                        JDATE, JTIME, DAY_INDEX ) ) THEN
                       MESG = 'Could not read "INDXD" from file "' //
     &                     TRIM( DAYNAME ) // '".'
                       CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                END IF

            END IF

C.............  Read time-dependent ZF and ZH for hydrostatic Met data
C.............  Compute per-source heights
            IF( .NOT. XFLAG .AND. .NOT. ZSTATIC ) THEN

                CALL SAFE_READ3( XNAME,'ZH',ALLAYS3,JDATE,JTIME,XBUF )
                CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, ZH )

                CALL SAFE_READ3( XNAME,'ZF',ALLAYS3,JDATE,JTIME,XBUF )
                CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, ZF )

C.................  Pre-process ZF and ZH to compute DDZH and DDZF
                CALL COMPUTE_DELTA_ZS

            END IF

C.............  Read and transform meteorology:
            IF ( .NOT. XFLAG ) THEN
            CALL SAFE_READ3( SNAME, 'HFX', ALLAYS3, JDATE, JTIME, XBUF )
            CALL METCELL_2D( METNGRID, NSRC, GX, XBUF, HFX)

            CALL SAFE_READ3( SNAME, 'PBL', ALLAYS3, JDATE, JTIME, XBUF )
            CALL METCELL_2D( METNGRID, NSRC, GX, XBUF, HMIX)

            CALL GET_VARIABLE_NAME( 'TGD', VNAME )
            CALL SAFE_READ3( SNAME, VNAME, ALLAYS3, JDATE, JTIME, XBUF )
            CALL METCELL_2D( METNGRID, NSRC, GX, XBUF, TSFC)

            CALL SAFE_READ3( SNAME, 'USTAR', ALLAYS3, JDATE,JTIME,XBUF )
            CALL METCELL_2D( METNGRID, NSRC, GX, XBUF, USTAR)

            CALL SAFE_READ3( SNAME, 'PRSFC', ALLAYS3, JDATE,JTIME,XBUF )
            CALL METCELL_2D( METNGRID, NSRC, GX, XBUF, PRSFC)

            CALL SAFE_READ3( XNAME, 'TA', ALLAYS3, JDATE, JTIME, XBUF )
            CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, TA )

            CALL SAFE_READ3( XNAME, 'QV', ALLAYS3, JDATE, JTIME, XBUF )
            CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, QV )

            CALL SAFE_READ3( XNAME, 'PRES', ALLAYS3, JDATE, JTIME,XBUF )
            CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, PRES )

            CALL SAFE_READ3( XNAME, 'DENS', ALLAYS3, JDATE, JTIME,XBUF )
            CALL METCELL_3D( METNGRID, NSRC, EMLAYS, GX, XBUF, DENS )

C.............  Read the low left dot and interpolate with the upper right
C.............  Shift over one column
            NCEL = 0
            IF ( .NOT. GFLAG) THEN
                DO S = 1, NSRC
                    NCEL = GD( S )
                    GDS( S ) = NCEL + 1 
                END DO
            END IF
            CALL SAFE_READ3( DNAME, 'UWINDC', ALLAYS3, JDATE,JTIME,DBUF )
            CALL METDOT_3D( NDOTS, NSRC, EMLAYS, GD, GDS, DBUF, UWIND )

C.............  Shift over one row on the dot grid
            NCEL = 0
            IF ( .NOT. GFLAG) THEN
                DO S = 1, NSRC
                    NCEL = GD( S )
                    GDS( S ) = NCEL + (NCOLS + 1)
                END DO
            END IF
            CALL SAFE_READ3( DNAME, 'VWINDC', ALLAYS3, JDATE,JTIME,DBUF )
            CALL METDOT_3D( NDOTS, NSRC, EMLAYS, GD, GDS, DBUF, VWIND )
            END IF

C.............  Precompute constants before starting source loop
            P  = ( SIGH(0) - VGLVSXG(0) ) / ( SIGH( 1 ) - SIGH( 0 ) )

C.............  Loop through sources and compute plume rise
            K = 0

            DO S = 1, NSRC

                IF( .NOT. FIREFLAG .OR. .NOT. AIRFLAG ) THEN
                    DM = STKDM( S )
                    HT = STKHT( S )
                    TK = STKTK( S )
                    VE = STKVE( S )
                    FL = 0.          ! initialize flow
                END IF

                XL = XLOCA( S )
                YL = YLOCA( S )

C.................  Skip fire source if not in day-specific file
                IF( FIREFLAG .AND. HOURFIRE ) THEN

                   MY_INDEX = -1
                   DO MY_LOOP = 1, DAY_NSRC
                       IF(S .EQ. DAY_INDEX(MY_LOOP)) THEN
                           MY_INDEX = MY_LOOP
                       ENDIF
                   ENDDO
                   IF(MY_INDEX .LE. 0) CYCLE
                ENDIF

C.................  Put smoldering fire sources in layer 1 and skip plume rise
                IF(( FIREFLAG ) ) THEN
                    IF( CSCC(S)(9:9) .EQ. 'S' ) THEN
                        LFRAC( S,2:EMLAYS) = 0.0
                        LFRAC( S,1 ) = 1.0
                        CYCLE
                    ENDIF
                ENDIF

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

C.................  Processing EDMS aircraft emission inventory
                ELSE IF( AIRFLAG ) THEN

C.....................  Compute pressures, use sigma values and surface pressures
                    DO L = EMLAYS, 0, -1
                        PRESF( L ) = VGLVSXG( L ) *
     &                               ( PRSFC( S ) - VGTOP ) * CONVPA +
     &                               VGTOP * CONVPA
                    END DO

C.....................  Convert surface pressure from Pa to mb
                    PSFC = PRSFC( S ) * CONVPA

C.....................  Compute derived met vars needed before layer assignments
                    CALL PREPLM( EMLAYS, HMIX( S ), STKHT( S ), PSFC,
     &                       TSFC( S ), DDZF( 1,S ), QV( 1,S ),
     &                       TA( 1,S ), UWIND( 1,S ), VWIND( 1,S ),
     &                       ZH( 1,S ), ZF( 1,S ), ZSTK( 1,S ),
     &                       PRESF( 1 ), LSTK, LPBL, TSTK, WSTK,
     &                       DTHDZ, WSPD, ZZF )

C.....................  Retreive altitude from source characterisitcs and
C                       convert ft to meter
                    ALTITUDE = STKHT( S )
                    ZBOT = ALTITUDE
                    ZTOP = ALTITUDE

C.....................  Setup for computing plume fractions, assuming uniform
C                       distribution in pressure (~mass concentration -- minor
C                       hydrostatic assumption) from bottom to top.
                    SURFACE = PSFC
                    LFULLHT( 0:EMLAYS ) = ZZF ( 0:EMLAYS   )
                    LHALFHT( 1:EMLAYS ) = ZH  ( 1:EMLAYS,S )
                    WEIGHTS( 1:EMLAYS ) = PRES( 1:EMLAYS,S ) * CONVPA
                    TEMPS  ( 1:EMLAYS ) = TA  ( 1:EMLAYS,S )

C.................  For non-explicit plume rise, preprocess met data...
                ELSE

C.....................  Compute pressures, use sigma values and surface pressures
                    DO L = EMLAYS, 0, -1
                        PRESF( L ) = VGLVSXG( L ) *
     &                               ( PRSFC( S ) - VGTOP ) * CONVPA +
     &                               VGTOP * CONVPA
                    END DO

C.....................  Convert surface pressure from Pa to mb
                    PSFC = PRSFC( S ) * CONVPA

C.....................  Compute derived met vars needed before layer assignments
                    IF( FIREFLAG ) THEN
                        CALL FIRE_PREPLM( EMLAYS, HMIX( S ), ZERO, PSFC,
     &                           TSFC( S ), DDZF( 1,S ), QV( 1,S ),
     &                           TA( 1,S ), UWIND( 1,S ), VWIND( 1,S ),
     &                           ZH( 1,S ), ZF( 1,S ), ZSTK( 1,S ),
     &                           PRESF( 1 ), LSTK, LPBL, TSTK, WSTK,
     &                           DTHDZ, WSPD, ZZF )
                    ELSE
                        CALL PREPLM( EMLAYS, HMIX( S ), STKHT(S), PSFC,
     &                           TSFC( S ), DDZF( 1,S ), QV( 1,S ),
     &                           TA( 1,S ), UWIND( 1,S ), VWIND( 1,S ),
     &                           ZH( 1,S ), ZF( 1,S ), ZSTK( 1,S ),
     &                           PRESF( 1 ), LSTK, LPBL, TSTK, WSTK,
     &                           DTHDZ, WSPD, ZZF )
                    END IF

C.....................  Trap USTAR at a minimum realistic value
                    USTMP = MAX( USTAR( S ), USTARMIN )

C.....................  Convert heat flux from watts/m^2 to m K/s
                    HFX( S ) = HFX( S ) / ( CP * DENS( 1,S ) )

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

                        IF( FIREFLAG ) THEN
                            IF( HOURFIRE ) THEN
                                TMPBFLX = BFLX( S )
                            END IF

                            CALL FIRE_PLMRIS( EMLAYS, LPBL, LSTK,
     &                           HFX(S), HMIX(S), TMPBFLX, TSTK, USTMP,
     &                           DTHDZ, TA(1,S), WSPD, ZZF(0), ZH(1,S),
     &                           ZSTK(1,S), WSTK, ZTOP, ZBOT, ZPLM )
                        ELSE
                            CALL PLMRIS( EMLAYS, LPBL, LSTK, HFX(S),
     &                           HMIX(S), DM, HT, TK, VE, TSTK, USTMP,
     &                           DTHDZ, TA(1,S), WSPD, ZZF(0), ZH(1,S),
     &                           ZSTK(1,S), WSTK, ZPLM )
                        END IF

C.........................  Determine bottom and top heights of the plume
                        IF( IPVERT == 0 ) THEN

C.............................  Default Turner approach: plume thickness =
C                               amount of plume rise
                            IF( FIREFLAG ) THEN
                                ZTOP = 1.5 * ZPLM
                                ZBOT = 0.5 * ZPLM
                            ELSE
                                ZTOP = STKHT( S ) +
     &                                 1.5 * ( ZPLM - STKHT( S ) )
                                ZBOT = STKHT( S ) +
     &                                 0.5 * ( ZPLM - STKHT( S ) )
                            END IF
                        ELSE

C.............................  Otherwise, compute plume top and bottom heights
                            CALL PLSPRD( DTHDZ, ZZF, EMLAYS, ZPLM,
     &                                   HMIX(S), ZTOP, ZBOT )
                        END IF
                    END IF

C.....................  Setup for computing plume fractions, assuming uniform
C                       distribution in pressure (~mass concentration -- minor
C                       hydrostatic assumption) from bottom to top.
                    SURFACE = PSFC
                    LFULLHT( 0:EMLAYS ) = ZZF ( 0:EMLAYS   )
                    LHALFHT( 1:EMLAYS ) = ZH  ( 1:EMLAYS,S )
                    WEIGHTS( 1:EMLAYS ) = PRES( 1:EMLAYS,S ) * CONVPA
                    TEMPS  ( 1:EMLAYS ) = TA  ( 1:EMLAYS,S )

                END IF  ! if computing plume rise

C.................  Check plume rise for nonsense values
                CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                IF( ZTOP .LT. STKHT( S ) .AND. K .LE. 0 ) THEN
                    MESG = 'WARNING: Top of plume found to be ' //
     &                     'less than top of stack for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                END IF

C.................  Allocate plume to layers
                IF( FIREFLAG ) THEN

                    IF( HOURFIRE ) THEN

                      MY_INDEX = -1
                      DO MY_LOOP = 1, DAY_NSRC
                          IF (S .EQ. DAY_INDEX(MY_LOOP)) THEN
                              MY_INDEX = MY_LOOP
                          ENDIF
                      ENDDO

                      IF( MY_INDEX .GT. 0) THEN
                          IF( DAY_ACRES( MY_INDEX ) .GT. 0.0) THEN
                            TMPACRE = DAY_ACRES( MY_INDEX )

                          ELSE
                            WRITE( MESG,94010 )'NOTE: Fire Has Zero ' //
     &                      'Acres: all emissions assigned to layer 1'//
     &                      CRLF() // BLANK10 // 'SOURCE NUMBER = ', S,
     &                      ' : ' // 'FIPS CODE = ' // CIFIP( S ) // ' ' //
     &                     CRLF()//BLANK10//BUFFER( 1:L2 )
                           CALL M3MSG2( MESG )

                           TMPACRE = 0.0
                         ENDIF

                      ELSE
                         MESG = 'ERROR: Fire Source not found '//
     &                       CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                      ENDIF
                    END IF

                    CALL FIRE_POSTPLM( EMLAYS, S, ZBOT, ZTOP, PRESF,
     &                                 LFULLHT, TEMPS, LHALFHT, TMPACRE,
     &                                 FPB1FLAG, LTOP, TFRAC )

C.....................  Calculate layer fraction for fire
C                       First, layer fractions for LAY1F under PBL
C                       Second, the rest of (1-LAY1F) gets distributed to above PBL
                    IF( LFLAG .AND. K .GT. 0)  THEN
                    IF( LAY1F( K ) .GT. 0. ) THEN
                        TFRAC = 0.0
                        PFRAC = 0.0
                        PBL = HMIX( S )
                        IF( PBL  .LE. 100.0 ) PBL = 100.0    ! reset PBL if < 100.
                        IF( ZTOP .LE. ZBOT ) THEN           ! Error: ZTOP < ZBOT
                            MESG = 'ERROR: PTOP can not be lower than PBOT'//
     &                          CRLF() // BLANK10 // BUFFER( 1:L2 )
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF

C.........................  Calculate layer fraction below PBL height
C                           Renormalize layer fractions based on LAY1F
                        TBOT = 0.0                ! Zbot is groud
                        TTOP = PBL                ! Ztop is PBL height

                        CALL POSTPLM( EMLAYS, S, TBOT, TTOP, PRESF,
     &                      LFULLHT, TEMPS, LHALFHT, LBOT, LTOP, PFRAC )

                        PFRAC( LBOT:LTOP ) = PFRAC( LBOT:LTOP ) * LAY1F( K )
                        TFRAC = TFRAC + PFRAC

C.........................  Calculate layer fraction above PBL height
C                           Renormalize layer fractions based on 1-LAY1F
                        IF( ZBOT < PBL .AND. ZTOP < PBL ) THEN
                            TBOT = ZBOT
                            TTOP = PBL
                        ELSE IF( ZBOT < PBL .AND. ZTOP > PBL ) THEN
                            TBOT = ZBOT
                            TTOP = ZTOP
                        ELSE IF( ZBOT > PBL .AND. ZTOP > PBL ) THEN
                            TBOT = PBL
                            TTOP = ZTOP
                        ELSE
                            TBOT = ZBOT
                            TTOP = ZTOP
                        END IF

                        CALL POSTPLM( EMLAYS, S, TBOT, TTOP, PRESF,
     &                      LFULLHT, TEMPS, LHALFHT, LBOT, LTOP, PFRAC )

                        TDIFF = 1.0 - LAY1F( K )
                        PFRAC( LBOT:LTOP ) = PFRAC( LBOT:LTOP ) * TDIFF
                        TFRAC = TFRAC + PFRAC
                    ELSE
                        MESG = 'ERROR: LAY1F value can not be zero '//
     &                       CRLF() // BLANK10 // BUFFER( 1:L2 )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF
                    END IF

                ELSE
                    CALL POSTPLM( EMLAYS, S, ZBOT, ZTOP, PRESF,
     &                  LFULLHT, TEMPS, LHALFHT, LBOT, LTOP, TFRAC )

C.....................  If hourly layer-1 fraction is present, reset this and re-
C                       normalize
C.....................  Must account for the case where LAY1F value is
C                       missing
                    IF( LFLAG .AND. K .GT. 0 ) THEN
                        IF( LAY1F( K ) .GT. 0. .AND.
     &                      TFRAC( 1 ) .LT. 1.       ) THEN
                            TSUM = SUM( TFRAC( 2:EMLAYS ) )
                            TDIFF = TSUM + TFRAC( 1 ) - LAY1F( K )
                            FAC = TDIFF / TSUM

                            TFRAC( 1 ) = LAY1F( K )
                            TFRAC( 2:EMLAYS ) = TFRAC( 2:EMLAYS ) * FAC
                        END IF
                    END IF

                END IF

C.................  Check if layer fractions are negative and reset
C                   to output in the first layer if they are.
                X = MINVAL( TFRAC( 1:EMLAYS ) )
                IF( X .LT. 0.0 ) THEN

                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                    MESG = 'WARNING: One or more negative plume ' //
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
                    WRITE( RDEV,OUTFMT ) S, CIFIP( S ),
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

93042   FORMAT( '( I6, ",", A, ",", A', I2.2, ', ","', I2.2, '(A',
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

C.............  This internal subprogram tries to retrieve the I/O API header
C               and aborts if it was not successful
            SUBROUTINE RETRIEVE_SET_HEADER( FILNAM )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C----------------------------------------------------------------------

            IF ( .NOT. DESCSET( FILNAM,-1 ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            END SUBROUTINE RETRIEVE_SET_HEADER

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
                OUTNAME = TGNAM
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

                IF( VARNAM == 'TEMP2' .OR. VARNAM == 'TEMP1P5' ) THEN
                    MESG = 'Please reset PLUME_GTEMP_NAME to match ' //
     &                 'to a variable name from file '// FILNAM( 1:L2)
                    CALL M3MSG2( MESG )
                END IF

                MESG = 'Could not read "' // VARNAM( 1:L1 ) //
     &                 '" from file "' // FILNAM( 1:L2 ) // '".'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )


            END IF

            END SUBROUTINE SAFE_READ3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram computes DDZH and DDZF
            SUBROUTINE COMPUTE_DELTA_ZS

C----------------------------------------------------------------------

            DO S = 1, NSRC

                ZZ0 = ZH( 1,S )
                DDZH( 1,S ) = 1.0 / ZZ0
                ZSTK( 1,S ) = ZF( 1,S ) - STKHT( S )
                ZF0 = ZF( 1,S )
                DDZF( 1,S ) = 1.0 / ZF0

                DO L = 2, EMLAYS

                    ZZ1 = ZH( L,S )
                    ZSTK( L,S ) = ZF( L,S ) - STKHT( S )
                    DDZH( L,S ) = 1.0 / ( ZZ1 - ZZ0 )
                    ZZ0 = ZZ1
                    ZF1 = ZF( L,S )
                    DDZF( L,S ) = 1.0 / ( ZF1 - ZF0 )
                    ZF0 = ZF1

                END DO

            END DO  ! End processing for intermediate layer height calcs

            END SUBROUTINE COMPUTE_DELTA_ZS

C----------------------------------------------------------------------

C.............  This internal subprogram computes assigns the grid cell
C                value to the source
            SUBROUTINE METCELL_2D( M, N, IX, V, Y )

C.............  Subprogram arguments
            INTEGER, INTENT (IN)  :: M              ! length of input
            INTEGER, INTENT (IN)  :: N              ! length of output array
            INTEGER, INTENT (IN)  :: IX( N )        ! index array
            REAL,    INTENT (IN)  :: V( M )         ! input vector
            REAL,    INTENT (OUT) :: Y( N )         ! output vector
C.............  Local variables
            INTEGER S, J

C----------------------------------------------------------------------

            DO S = 1, N
               
                J = IX( S ) 
                Y( S ) = V( J )

            END DO

            END SUBROUTINE METCELL_2D

C----------------------------------------------------------------------

C.............  This internal subprogram computes assigns the grid cell
C                value to the source
            SUBROUTINE METCELL_3D( M, N, P, IX, V, Y )

C.............  Subprogram arguments
            INTEGER, INTENT (IN)  :: M              ! length of input
            INTEGER, INTENT (IN)  :: N              ! length of output array
            INTEGER, INTENT (IN)  :: P              ! number of layers
            INTEGER, INTENT (IN)  :: IX( N )        ! index array
            REAL,    INTENT (IN)  :: V( M, P )         ! input vector
            REAL,    INTENT (OUT) :: Y( P, N )         ! output vector
C.............  Local variables
            INTEGER S, J, L

C----------------------------------------------------------------------

            DO S = 1, N
               
                J = IX( S ) 

                DO L = 1, P

                    Y( L, S ) = V( J, L )

                END DO

            END DO

            END SUBROUTINE METCELL_3D

C.............  This internal subprogram assigns the grid cell 
C                value to the source from the interpolated grid dot
            SUBROUTINE METDOT_3D( M, N, P, JX, KX, V, Y )

C.............  Subprogram arguments
            INTEGER, INTENT (IN)  :: M              ! length of input
            INTEGER, INTENT (IN)  :: N              ! length of output array
            INTEGER, INTENT (IN)  :: P              ! number of layers
            INTEGER, INTENT (IN)  :: JX( N )        ! index array dot 1
            INTEGER, INTENT (IN)  :: KX( N )        ! index array dot 2
            REAL,    INTENT (IN)  :: V( M, P )         ! input vector
            REAL,    INTENT (OUT) :: Y( P, N )         ! output vector
C.............  Local variables
            INTEGER S, J, K, L
            REAL LL, UR

C----------------------------------------------------------------------

            DO S = 1, N
               
                J = JX( S )
                K = KX( S ) 

                DO L = 1, P

                    LL = V( J, L )
                    UR = V( K, L )
                    Y( L, S ) = 0.5 * (LL + UR)

                END DO

            END DO

            END SUBROUTINE METDOT_3D

        END PROGRAM LAYPOINT

