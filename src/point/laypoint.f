
        PROGRAM LAYPOINT

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
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
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C env_progs@mcnc.org
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

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKMETEM
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        CHARACTER*50    GETCFDSC
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        CHARACTER*16    VERCHAR
        INTEGER         WKDAY

        EXTERNAL        CHKMETEM, CRLF, ENVINT, GETCFDSC, HHMMSS, INDEX1,
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, VERCHAR, WKDAY

C...........  LOCAL PARAMETERS and their descriptions:

        REAL        , PARAMETER :: USTARMIN  = 0.1  ! Min valid value for USTAR

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

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

C.........  Fixed-dimension arrays
        REAL         SIGH   ( 0:MXLAYS3-1 )   !  half-level sigma values
        REAL         VGLVSXG( 0:MXLAYS3 )     !  vertical coord values.

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  log file
        INTEGER      :: PDEV = 0!  elevated/PinG source file
        INTEGER      :: RDEV = 0!  optional report iff REP_LAYER_MAX is set
        INTEGER         SDEV    !  ASCII part of inventory file

        CHARACTER*16    ANAME   !  ASCII point-source inventory file
        CHARACTER*16    DNAME   !  dot-point layered met file name
        CHARACTER*16    ENAME   !  point-source inventory input file
        CHARACTER*16    GNAME   !  cross-point layered grid file name
        CHARACTER*16    LNAME   !  layer fractions matrix output file
        CHARACTER*16    SNAME   !  cross-point surface met file name
        CHARACTER*16    XNAME   !  cross-point layered met file name

C...........   Other local variables

        INTEGER          I, L, L1, L2, S, T  ! counters and indices

        INTEGER          EMLAYS    ! number of emissions layers
        INTEGER          GDTYP     ! i/o api grid type code
        INTEGER          IOS       ! tmp i/o status
        INTEGER          JDATE     ! Julian date (YYYYDDD)
        INTEGER          JTIME     ! time (HHMMSS)
        INTEGER          LBOT      ! plume bottom layer
        INTEGER          LDATE     ! previous date
        INTEGER          LPBL      ! first L: ZF(L) above mixing layer
        INTEGER          LSTK      ! first L: ZF(L) > STKHT
        INTEGER          LTOP      ! plume top    layer
        INTEGER          NCOLS     ! cross grid number of grid columns
        INTEGER          NDOTS     ! dot grid number of cells
        INTEGER          NGRID     ! cross grid number of cells
        INTEGER          NLAYS     ! number of layers in met files
        INTEGER          NMAJOR    ! no. major sources
        INTEGER          NPING     ! no. plume-in-grid sources
        INTEGER          NROWS     ! cross grid number of grid rows
        INTEGER          NSTEPS    ! mumber of time steps
        INTEGER          REP_LAYR  ! layer for reporting srcs w/ high plumes
        INTEGER          SDATE     ! Julian start date (YYYYDDD)
        INTEGER          STIME     ! start time (HHMMSS)
        INTEGER          TSTEP     ! output time step
        INTEGER          VGTYP     ! vertical type from cross-point file

        REAL             X, Y, P, Q, PP, QQ
        REAL             XBEG, XEND, XL  ! tmp x-coords
        REAL             YBEG, YEND, YL  ! tmp y-coords
        REAL             PSFC      !  surface pressure (Pa)
        REAL             TSTK      !  temperature at top of stack (K)
        REAL             USTMP     !  tmp Ustar
        REAL             VGTOP     !  top value from cross-point file 
        REAL             WSTK      !  wind speed  at top of stack (m/s)
        REAL             ZZ0, ZZ1, ZF0, ZF1
        REAL             ZBOT      !  plume bottom elevation (m)
        REAL             ZTOP      !  plume top    elevation (m)

        REAL*8           P_ALP     ! projection alpha
        REAL*8           P_BET     ! projection beta
        REAL*8           P_GAM     ! projection gamma
        REAL*8           XCELL     ! cross grid X-coordinate cell dimension
        REAL*8           YCELL     ! cross grid Y-coordinate cell dimension
        REAL*8           XCENT     ! x-center of projection
        REAL*8           YCENT     ! y-center of projection
        REAL*8           XORIG     ! cross grid X-cord origin of grid 
        REAL*8           YORIG     ! cross grid Y-coordinate origin of grid
        REAL*8           XCELLDG   ! dot grid X-coordinate cell dimension
        REAL*8           YCELLDG   ! dot grid Y-coordinate cell dimension
        REAL*8           XORIGDG   ! dot grid X-coordinate origin of grid 
        REAL*8           YORIGDG   ! dot grid Y-coordinate origin of grid

        LOGICAL       :: EFLAG = .FALSE.  ! error flag
        LOGICAL       :: VFLAG = .FALSE.  ! true: use elevated/PinG file (PELV)
        LOGICAL          LF( 9 )          ! true: source characteristic is valid

        CHARACTER*50     CHARS( 9 )!  tmp source characeristics 
        CHARACTER*50  :: METSCEN   !  temporary string for met scenario name
        CHARACTER*50  :: CLOUDSHM  !  temporary string for cloud scheme name
        CHARACTER*200    OUTFMT    !  output format for RDEV report
        CHARACTER*200    BUFFER    !  source characteristics buffer
        CHARACTER*300    MESG      !  buffer for M3EXIT() messages

        CHARACTER(LEN=IOVLEN3)  VNAME ! variable name buffer 
        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER*16  :: PROGNAME = 'LAYPOINT'   !  program name

C***********************************************************************
C   begin body of program LAYPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........   Get setting from environment variables
        EMLAYS = ENVINT( 'SMK_EMLAYS', 'Number of emission layers',
     &                   -1, IOS )

        MESG = 'Indicator for create plume-in-grid outputs'
        VFLAG = ENVYN( 'SMK_PING_YN', MESG, .FALSE., IOS )

        MESG = 'Indicator for defining major/minor sources'
        VFLAG = ( VFLAG .OR. 
     &            ENVYN( 'SMK_SPECELEV_YN', MESG, .FALSE., IOS ) )


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

        IF( VFLAG ) THEN
            PDEV = PROMPTFFILE( 
     &          'Enter logical name for the ELEVATED POINT SOURCE file',
     &          .TRUE., .TRUE., CRL // 'ELV', PROGNAME )
        END IF

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

C.........  Check multiple met files for consistency
        EFLAG = ( .NOT. CHKMETEM( 'NONE', SNAME, GNAME, XNAME, DNAME ) )

C.........  Get grid parameters from 3-d cross-point met file and store needed
C           header information.  Use time parameters for time defaults.
        CALL RETRIEVE_IOAPI_HEADER( XNAME )
        SDATE  = SDATE3D
        STIME  = STIME3D
        NSTEPS = MXREC3D
        NCOLS  = NCOLS3D
        NROWS  = NROWS3D
        NLAYS  = NLAYS3D
        GDTYP  = GDTYP3D
        P_ALP  = P_ALP3D
        P_BET  = P_BET3D
        P_GAM  = P_GAM3D
        XCENT  = XCENT3D
        YCENT  = YCENT3D
        XORIG  = XORIG3D
        YORIG  = YORIG3D
        XCELL  = XCELL3D
        YCELL  = YCELL3D
        VGTYP  = VGTYP3D
        VGTOP  = VGTOP3D

        NGRID = NCOLS * NROWS
        NDOTS = ( NCOLS + 1 ) * ( NROWS + 1 )

        VGLVSXG( 0 ) = VGLVS3D( 0 )
        DO I = 1, NLAYS3D
            VGLVSXG( I ) = VGLVS3D( I )
            SIGH   ( I-1 ) = 0.5 * ( VGLVS3D( I ) + VGLVS3D( I-1 ) )
        END DO

        METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. ) 
        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. ) 

        CALL RETRIEVE_IOAPI_HEADER( DNAME )
        XCELLDG = XCELL3D
        YCELLDG = YCELL3D 
        XORIGDG = XORIG3D
        YORIGDG = YORIG3D

C.........  Compare number of meteorology layers to number of emissions layers
        IF( EMLAYS .LE. NLAYS ) THEN
            WRITE( MESG,94010 ) 'NOTE: The number of emission layers '//
     &             'is', EMLAYS, ', and the number of '// CRLF()// 
     &             BLANK10//'layers in the meteorology file is', NLAYS
            CALL M3MSG2( MESG )

        ELSE
            WRITE( MESG,94010 ) 'Resetting number of emission layers '//
     &             'from', EMLAYS, 'to number of '// CRLF()// BLANK10 //
     &             'layers in the meteorology file,', NLAYS
            CALL M3WARN( PROGNAME, 0, 0, MESG )

            EMLAYS = NLAYS

        END IF

C.........  Update start date/time and duration from the environment
        CALL GETM3EPI( -1, SDATE, STIME, NSTEPS )
        TSTEP = 10000

C.........  Set up and open output file, which will primarily using I/O API 
C           settings from the cross-point met file (XNAME), which are 
C           already retrieved

        CALL HDRMISS3 

        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NROWS3D = NSRC
        NLAYS3D = EMLAYS
        NVARS3D = 1
        VGLVS3D( 1:EMLAYS ) = VGLVSXG( 1:EMLAYS )  ! array
        VGTYP3D = VGTYP
        VGTOP3D = VGTOP

        VNAME3D( 1 ) = 'LFRAC'
        VTYPE3D( 1 ) = M3REAL
        UNITS3D( 1 ) = 'none'
        VDESC3D( 1 ) = 'Fraction of plume emitted into layer'

        FDESC3D = ' '  ! array

        FDESC3D( 1 ) = 'Source level hourly plume rise layer fractions'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
        FDESC3D( 4 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 5 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

        FDESC3D( 11 ) = '/PNTS FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/PNTS VERSION/ ' // IFDESC3

        MESG = 'Enter logical name for LAYER FRACTIONS MATRIX'
        LNAME = PROMPTMFILE( MESG, FSUNKN3, 'PLAY', PROGNAME )

C.........  Get file name of report of plume exceeding specified layer
C.........  Write header to the report
        IF( REP_LAYR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'Enter logical name for report of ' //
     &                          'plumes exceeding layer', REP_LAYR
            RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                          'REPRTLAY', PROGNAME )

        ENDIF

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( 'POINT', ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Call subroutine to convert grid coordinates from lat-lon to
C           coordinate system of the destination grid
        CALL CONVRTXY( NSRC, GDTYP, P_ALP, P_BET, P_GAM, 
     &                 XCENT, YCENT, XLOCA, YLOCA )

C.........  Call elevated sources indicator file, even thought it might not
C           be opened - routine will initialize LMAJOR and LPING regardless
C           of whether the file is available.
        CALL RDPELV( PDEV, NSRC, NMAJOR, NPING )

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

C.........  Allocate layer fractions array
        ALLOCATE( LFRAC( NSRC,EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )

C.........  Allocate ungridding arrays
        ALLOCATE( ND( 4,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ND', PROGNAME )
        ALLOCATE( NX( 4,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NX', PROGNAME )
        ALLOCATE( CD( 4,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CD', PROGNAME )
        ALLOCATE( CX( 4,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CX', PROGNAME )

C.........  Allocate per-layer arrays from 1:EMLAYS
        ALLOCATE( WSPD( EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WSPD', PROGNAME )
        ALLOCATE( DTHDZ( EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DTHDZ', PROGNAME )
        ALLOCATE( TFRAC( EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFRAC', PROGNAME )

C.........  Allocate per-layer arrays from 0:EMLAYS
        ALLOCATE( PRESF( 0:EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRESF', PROGNAME )
        ALLOCATE( ZZF( 0:EMLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ZZF', PROGNAME )

C.........  Allocate array for tmp gridded, layered cross-point met data
        ALLOCATE( XBUF( NGRID,NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XBUF', PROGNAME )
        ALLOCATE( DBUF( NDOTS,NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DBUF', PROGNAME )

C.........  Compute un-gridding matrices for dot and cross point met data
        CALL UNGRIDB( NCOLS+1, NROWS+1, 
     &                XORIGDG, YORIGDG, XCELLDG, YCELLDG,
     &                NSRC, XLOCA, YLOCA, ND, CD )

        CALL UNGRIDB( NCOLS, NROWS, XORIG, YORIG, XCELL, YCELL,
     &                NSRC, XLOCA, YLOCA, NX, CX )

C.........  Read time-independent ZF and ZH
C.........  Use BMATVEC to convert from a gridded array to a source-based array

        CALL RETRIEVE_IOAPI_HEADER( GNAME )
        CALL GET_VARIABLE_NAME( 'ZH', VNAME )
        CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D, STIME3D, XBUF )
        CALL BMATVEC( NGRID, NSRC, EMLAYS, NX, CX, XBUF, ZH )

        CALL GET_VARIABLE_NAME( 'ZF', VNAME )
        CALL SAFE_READ3( GNAME, VNAME, ALLAYS3, SDATE3D, STIME3D, XBUF )
        CALL BMATVEC( NGRID, NSRC, EMLAYS, NX, CX, XBUF, ZF )

C.........  Pre-process ZF and ZH to compute DDZH and DDZF
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

        END DO  ! End processing for 

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
        LF( 1:NCHARS ) = .TRUE.   ! array
        IF( NCHARS .LE. 8 ) LF( NCHARS+1:9 ) = .FALSE.  ! array

C.........  Get variable names from surface meteorology file
        CALL RETRIEVE_IOAPI_HEADER( SNAME )

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
                ENDIF

                LDATE = JDATE
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

C.............  Write to report file if report feature is on
            IF( RDEV .GT. 0 ) THEN
                WRITE( RDEV,93020 ) HHMMSS( JTIME )
            END IF

C.............  Read and transform meteorology:

            CALL SAFE_READ3( SNAME, 'HFX', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( NGRID, NSRC, 1, NX, CX, XBUF, HFX )

            CALL SAFE_READ3( SNAME, 'PBL', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( NGRID, NSRC, 1, NX, CX, XBUF, HMIX )

            CALL GET_VARIABLE_NAME( 'TGD', VNAME )
            CALL SAFE_READ3( SNAME, VNAME, ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( NGRID, NSRC, 1, NX, CX, XBUF, TSFC )

            CALL SAFE_READ3( SNAME, 'USTAR', ALLAYS3, JDATE,JTIME,XBUF )
            CALL BMATVEC( NGRID, NSRC, 1, NX, CX, XBUF, USTAR )

            CALL SAFE_READ3( XNAME, 'TA', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( NGRID, NSRC, EMLAYS, NX, CX, XBUF, TA )

            CALL SAFE_READ3( XNAME, 'QV', ALLAYS3, JDATE, JTIME, XBUF )
            CALL BMATVEC( NGRID, NSRC, EMLAYS, NX, CX, XBUF, QV )

            CALL SAFE_READ3( XNAME, 'PRES', ALLAYS3, JDATE, JTIME,XBUF )
            CALL BMATVEC( NGRID, NSRC, EMLAYS, NX, CX, XBUF, PRES )

            CALL SAFE_READ3( DNAME, 'UWIND', ALLAYS3, JDATE,JTIME,DBUF )
            CALL BMATVEC( NDOTS, NSRC, EMLAYS, ND, CD, DBUF, UWIND )

            CALL SAFE_READ3( DNAME, 'VWIND', ALLAYS3, JDATE,JTIME,DBUF )
            CALL BMATVEC( NDOTS, NSRC, EMLAYS, ND, CD, DBUF, VWIND )

C.............  Precompute constants before starting source loop

            P  = ( SIGH(0) - VGLVSXG(0) ) / ( SIGH( 1 ) - SIGH( 0 ) )
            PP = 1.0 + P 
            QQ =     - P 

C.............  Loop through sources and compute plume rise
            DO S = 1, NSRC

C.................  Skip source if it is outside grid
	        XL = XLOCA( S )
	        YL = YLOCA( S )
                IF( XL .LT. XBEG .OR. XL .GT. XEND .OR.
     &              YL .LT. YBEG .OR. YL .GT. YEND     ) CYCLE

C.................  Skip source if it is a PinG source and assign 0 to the
C                   layer fractions
                IF( LPING( S ) ) THEN
                    LFRAC( S, 1:EMLAYS ) = 0.
                    CYCLE
                END IF

C.................  Skip source if it is minor source and assign layer fractions
C                   that put all emissions in layer 1
                IF( .NOT. LMAJOR( S ) ) THEN
                    LFRAC( S,1 )        = 1.
                    LFRAC( S,2:EMLAYS ) = 0.
                    CYCLE
                END IF

C.................  Compute surface pressure (and convert to mb from Pa)
                PSFC = 1.0E-2 * ( PP * PRES( 1,S ) + QQ * PRES( 2,S ) )

C.................  Compute derived met variables needed before PLMRISs
                CALL PREPLM( EMLAYS, HMIX( S ), STKHT( S ), PSFC, 
     &                       TSFC( S ), DDZH( 1,S ), QV( 1,S ), 
     &                       TA( 1,S ), UWIND( 1,S ), VWIND( 1,S ), 
     &                       ZH( 1,S ), ZF( 1,S ), ZSTK( 1,S ), 
     &                       PRES( 1,S ), LSTK, LPBL, TSTK, WSTK,
     &                       DTHDZ, WSPD, ZZF )

C.................  Trap USTAR at a minimum realistic value
                USTMP = MAX( USTAR( S ), USTARMIN )

C.................  Compute plume rise for this source 
                CALL PLMRIS( EMLAYS, LPBL, LSTK, HFX( S ), HMIX( S ),
     &                       STKDM( S ), STKHT( S ), STKTK( S ), 
     &                       STKVE( S ), TSTK, USTMP, DTHDZ, TA( 1,S ), 
     &                       WSPD, ZZF( 0 ), ZH( 1,S ), ZSTK( 1,S ), 
     &                       WSTK, ZTOP, ZBOT )

C.................  Check plume rise for nonsense output
                IF( ZTOP .LT. STKHT( S ) ) THEN

                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 ) 
     &                     'WARNING: Top of plume found to be ' //
     &                     'less than top of stack for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 )

                    CALL M3MESG( MESG )

                END IF

C.................  Compute plume fractions, assuming uniform distribution
C.................  in pressure (~mass concentration -- minor hydrostatic
C.................  assumption) from bottom to top.
                CALL POSTPLM( EMLAYS, S, PSFC, ZBOT, ZTOP, PRES( 1,S ),
     &                        ZZF, ZH( 1,S ), LTOP, TFRAC )

                LFRAC( S,1:EMLAYS ) = TFRAC( 1:EMLAYS )  ! array

C.................  Check if LTOP out of range, and report (will only work
C.................    if REP_LAYR env var has been set b/c default is -1
                IF( REP_LAYR .GT. 0 .AND. LTOP .GT. REP_LAYR ) THEN

                    CALL PARSCSRC( CSOURC( S ), NCHARS, SC_BEGP, 
     &                             SC_ENDP, LF, I, CHARS )

                    WRITE( OUTFMT, 93042 ) PLTLEN3, NCHARS-2, CHRLEN3
                    WRITE( RDEV,OUTFMT ) S, IFIP( S ),
     &                   ( CHARS( I ), I = 2,NCHARS ), STKHT( S ), 
     &                     STKVE ( S ), STKTK( S ), TSTK, 
     &                     WSTK, LPBL, LTOP 

                END IF

            END DO    !  end loop on sources S

            IF( EFLAG ) THEN

                MESG = 'Problem during plume rise calculations.'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            ELSEIF ( .NOT. WRITE3( LNAME, 'LFRAC', 
     &                             JDATE, JTIME, LFRAC ) ) THEN

                MESG = 'Problem writing "LFRAC" to file "' // 
     &                 LNAME( 1:LEN_TRIM( LNAME ) ) // '."'
     
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

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

            ENDIF

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
                OUTNAME = 'X3HT0F'
            CASE( 'ZF' )
                OUTNAME = 'X3HT0M'
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

        END PROGRAM LAYPOINT

