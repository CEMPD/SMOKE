
        PROGRAM LAYPOINT

C***********************************************************************
C  program body starts at line 271
C
C  DESCRIPTION:
C       Produces source-level layer fractions, using the emissions-source 
C       vector attributes produced by RAWPOINT and ambient meteorology.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level point source emissions in canonical order
C       produced by RAWPOINT
C       Nonhydrostatic sigma-P EDSS/M3 meteorology.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND2, FIND3, GETNUM, GETYN, DAYMON, TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  11/95 by CJC using algorithms adapted from RADM
C       subroutine PR()
C       Version 6/96 by CJC allocates plume fraction uniformly in PRES
C       instead of in Z (as RADM does), in order to conserve mixing ratio.
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

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'GRDIMS3.EXT'   !  gridding parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
        INCLUDE 'VGRD3.EXT'     !  vertical allocation parameters


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         DSCGRID
        INTEGER         GETDATE
        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        LOGICAL         LAMBERT
        LOGICAL         LL2LAM
        CHARACTER*14    MMDDYY
        REAL            POLY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY

        EXTERNAL        CRLF, ENVINT, ENVYN, DSCGRID, GETDATE, GETNUM, 
     &                  GETYN, HHMMSS, LAMBERT, LL2LAM, MMDDYY, POLY, 
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, WKDAY 

C...........   PARAMETERS and their descriptions

        CHARACTER*5 BLANK5
        INTEGER     NDOTS
        REAL        GAMA
        REAL        USTARMIN              ! Minimum realistic value for USTAR

        PARAMETER (  BLANK5 = ' '
     &             , NDOTS    = ( NCOLS + 1 ) * ( NROWS + 1 )
     &             , GAMA     = -0.0098 
     &             , USTARMIN = 0.01
     &            )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point source stack parameters:

        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     IPLANT( NPSRC )  !  source plant ID
        INTEGER     ISTACK( NPSRC )  !  source stack ID
        REAL        XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL        YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL        STKHT ( NPSRC )  !  stack height   (m)
        REAL        STKDM ( NPSRC )  !  stack diameter (m)
        REAL        STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE ( NPSRC )  !  exhaust velocity    (m/s)
 
C.......   Point source meteorology profiles and variables

        INTEGER     LSTK            !  first L: ZF(L) > STKHT
        INTEGER     LPBL            !  first L: ZF(L) above mixing layer

        REAL        TSFC ( NPSRC )  !  surface temperature (deg K)
        REAL        HFX  ( NPSRC )  !  sensible heat flux (M K / S )
        REAL        HMIX ( NPSRC )  !  mixing height (m)
        REAL        USTAR( NPSRC )  !  friction velocity (m/s)
        REAL        ZF   ( EMLAYS, NPSRC )  !  layer surface height (m)
        REAL        ZH   ( EMLAYS, NPSRC )  !  layer center  height (m)
        REAL        TA   ( EMLAYS, NPSRC )  !  temperature (K)
        REAL        TV   ( EMLAYS, NPSRC )  !  virtual temperature
        REAL        UWIND( EMLAYS, NPSRC )  !  wind speed (m/s)
        REAL        VWIND( EMLAYS, NPSRC )  !  wind speed (m/s)
        REAL        PRES ( EMLAYS, NPSRC )  !  pressure (Pa)
        REAL        PRESF( 0:EMLAYS )  !  pressure at full-levels
        REAL        ZZF  ( 0:EMLAYS )      !  elevation at full-levels
        REAL        DTHDZ( EMLAYS )        !  gradient of THETV
        REAL        DDZH ( EMLAYS,NPSRC )  !  1/( zh(l) - zh(l-1) )
        REAL        DDZF ( EMLAYS,NPSRC )  !  1/( zh(l) - zh(l-1) )
        REAL        ZSTK ( EMLAYS,NPSRC )  !  zf( l,s )   - stkht(s)   (m)
        REAL        WSPD ( EMLAYS )        !  wind speed (m/s)
        REAL        PSFC            !  surface pressure (Pa)
        REAL        QSFC            !  surface heat flux
        REAL        ES              !  water vapor pressure
        REAL        TVSFC           !  sfc virtual temperature
        REAL        THETG           !  sfc virtual potential temperature
        REAL        THV1, THVK      !  virtual potential temperature

        REAL        TSTK                  !  temperature at top of stack (deg K)
        REAL        WSTK                  !  wind speed  at top of stack (m/s)

C.......   Cross- and dot-point meteorology input buffers

        REAL        XBUF( NGRID, NLAYS )
        REAL        DBUF( NDOTS, NLAYS )

C.......   Dot point un-gridding matrix (uses bilinear interpolation)

        INTEGER         ND( 4,NPSRC )  !  cell indexes
        REAL            CD( 4,NPSRC )  !  coefficients

C.......   Cross point un-gridding matrix

        INTEGER         NX( 4,NPSRC )  !  cell indexes
        REAL            CX( 4,NPSRC )  !  coefficients

C.......   Output layer fractions (in I/O order)

        REAL            LFRAC( NPSRC, EMLAYS )


C.......   Dot and cross point grid description parameters:

        INTEGER      NCOLSDG      ! dot grid number of grid columns
        INTEGER      NROWSDG      ! dot grid number of grid rows
        INTEGER      NLAYSDG      ! dot grid number of layers
        INTEGER      NCOLSXG      ! cross grid number of grid columns
        INTEGER      NROWSXG      ! cross grid number of grid rows
        INTEGER      NLAYSXG      ! cross grid number of layers
        INTEGER      VGTYPXG      ! vertical coordinate type (VGSIGP3, ...)
        INTEGER      VGTYPDG      ! vertical coordinate type (VGSIGP3, ...)

        REAL         SIGH   ( 0 : MXLAYS3 - 1 )    !  half-level sigma values
        REAL         VGLVSDG( 0 : MXLAYS3 )  !  vertical coord values.
        REAL         VGLVSXG( 0 : MXLAYS3 )  !  vertical coord values.
        REAL         VGTOPDG      ! model-top, for sigma coord types.
        REAL         VGTOPXG      ! model-top, for sigma coord types.

        REAL*8       P_ALPDG      ! dot grid first, second, third map
        REAL*8       P_BETDG      ! dot grid projection descriptive
        REAL*8       P_GAMDG      ! dot grid parameters.
        REAL*8       XCENTDG      ! dot grid lon for coord-system X=0
        REAL*8       YCENTDG      ! dot grid lat for coord-system Y=0
        REAL*8       XORIGDG      ! dot grid X-coordinate origin of grid (map units)
        REAL*8       YORIGDG      ! dot grid Y-coordinate origin of grid
        REAL*8       XCELLDG      ! dot grid X-coordinate cell dimension
        REAL*8       YCELLDG      ! dot grid Y-coordinate cell dimension
        REAL*8       P_ALPXG      ! cross grid first, second, third map
        REAL*8       P_BETXG      ! cross grid projection descriptive
        REAL*8       P_GAMXG      ! cross grid parameters.
        REAL*8       XCENTXG      ! cross grid lon for coord-system X=0
        REAL*8       YCENTXG      ! cross grid lat for coord-system Y=0
        REAL*8       XORIGXG      ! cross grid X-cord origin of grid (map units)
        REAL*8       YORIGXG      ! cross grid Y-coordinate origin of grid
        REAL*8       XCELLXG      ! cross grid X-coordinate cell dimension
        REAL*8       YCELLXG      ! cross grid Y-coordinate cell dimension

        CHARACTER*14    DTBUF     ! Date buffer
        CHARACTER*16    COONMDG   ! dot grid coordinate system name
        CHARACTER*16    COONMXG   ! cross grid coordinate system name
        CHARACTER*16    GRDNMDG   ! dot grid grid name
        CHARACTER*16    GRDNMXG   ! cross grid grid name

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  unit number for log file
        INTEGER         RDEV    !  optional report iff REP_LAYER_MAX is set

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    DNAME   !  dot-point   met file name
        CHARACTER*16    LNAME   !  " timestepped layer fractions output file
        CHARACTER*16    SNAME   !  cross-point sfc met file name
        CHARACTER*16    XGRID   !  cross-point lay grid file name
        CHARACTER*16    XNAME   !  cross-point lay met file name

C...........   Other local variables

        REAL            DDP       !  1/(plume-thickness)
        REAL            PBOT      !  plume bottom sigma
        REAL            PDIFF     !  PBOT - PTOP
        REAL            PTOP      !  plume top    sigma
        REAL            X, Y, P, Q, PP, QQ
        REAL            XBEG, XEND, XL
        REAL            YBEG, YEND, YL
        REAL            ZZ0, ZZ1, ZF0, ZF1
        REAL            ZBOT      !  plume bottom elevation (m)
        REAL            ZTOP      !  plume top    elevation (m)

        INTEGER         I, L, M, S
        INTEGER         IOS
        INTEGER         JDATE, JTIME
        INTEGER         LBOT      !  plume bottom layer
        INTEGER         LDATE     !  previous date
        INTEGER         LTOP      !  plume top    layer
        INTEGER         NSTEPS
        INTEGER         REP_LAYR  !  Layer for reporting srcs w/ high plumes
        INTEGER         STATUS    !  dummy status for envint
        INTEGER         TSTEP     !  output time step

        LOGICAL         EFLAG     !  error flag
        LOGICAL         PROMPTF   !  iff PROMPTFLAG E.V. is true or not defined

        CHARACTER*80    SCEN_STR  !  temporary string for scenario name
        CHARACTER*80    SCHE_STR  !  temporary string for cloud scheme name
        CHARACTER*256   MESG      !  buffer for M3EXIT() messages


C.......   STATEMENT FUNCTION for evaluating relative error: TRUE
C.......   iff (PP,QQ) are significantly different (compares square of
C.......   normalized difference against square of tolerance)
 
        LOGICAL         DBLERR
        LOGICAL         FLTERR
        REAL*8          PD, QD
        REAL            RR, SS
 
        DBLERR( PD,QD ) = 
     &    ( ( PD - QD )**2 .GT. 1.0D-12*( PD*PD + QD*QD + 1.0D-5 ) ) 

        FLTERR( RR,SS ) = 
     &    ( ( RR - SS )**2 .GT. 1.0E-12*( RR*RR + SS*SS + 1.0E-5 ) ) 


C***********************************************************************
C   begin body of program LAYPOINT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &' ',
     &'Program LAYPOINT to take the point source emissions file from',
     &'RAWPOINT and meteorology from the EDSS dot-point and cross-',
     &'point meteorology files and construct the layered point source',
     &'plume fractions file for use by program CSGLMERGE.',
     &' '
        WRITE( *,92000 ) 
     &'You will need to enter the logical names for the input and',
     &'output files (and to have set them prior to program launch,',
     &'using "setenv <logicalname> <pathname>").  You will be ',
     &'prompted for which species to allocate -- these should only ',
     &'be the chemical species, and not the FIP-related or ASC-related',
     &'variables.)', 
     &' '
        WRITE( *,92000 ) 
     &'You will also be prompted for the TIME PERIOD to be covered',
     &'by this run of the program, as specified by STARTING DATE',
     &'in format YYYYDDD (= 1000*year + Julian day n0 (1...365,6) )',
     &'STARTING TIME in format HHMMSS (= 10000*hr + 100*min + sec) )',
     &'and DURATION in format HHMMSS.',
     &  ' '
        WRITE( *,92000 ) 
     &'You may use END_OF-FILE (control-D) to quit the program',
     &'during logical-name entry.  Default responses are indicated',
     &'in brackets [LIKE THIS].',
     &' '

        I = EMLAYS  ! So Alpha compiler will not give errors
        IF ( I .LT. 4 ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'EMLAYS must be at least 4', 2 )
        ELSE IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Ending program', 2 )
        END IF

C.......   Get environment variable logical variables
        REP_LAYR = ENVINT( 'REP_LAYER_MAX', 
     &                     'Layer number for reporting high plume rise',
     &                     9999, STATUS )

        PROMPTF = ENVYN ( 'PROMPTFLAG', 'Prompt for inputs or not',
     &                    .TRUE., IOS )

        IF( STATUS   .EQ. 0        .AND. 
     &    ( REP_LAYR .LT. 1 .OR.
     &      REP_LAYR .GT. EMLAYS )       ) THEN

            WRITE( MESG,94010 ) 'Environment variable REP_LAYR (=', 
     &                          REP_LAYR, ') is out of range 1 to',
     &                          EMLAYS, '.'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ENDIF

C.......   Get file name; open input point sources, temporal cross-reference,
C.......   and temporal profiles files

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'LAYPOINT' )

        DNAME = PROMPTMFILE( 
     &          'Enter name for DOT-POINT   LAYERED  MET file',
     &          FSREAD3, 'METDOT', 'LAYPOINT' )

        SNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT SURFACE  MET file',
     &          FSREAD3, 'SFCCRO', 'LAYPOINT' )

        XNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED  MET file',
     &          FSREAD3, 'LAYCRO', 'LAYPOINT' )

        XGRID = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED GRID file',
     &          FSREAD3, 'LAYGRD', 'LAYPOINT' )

C.........  Check XGRID dimensions - layers don't matter unless too small
        IF ( .NOT. DESC3( XGRID ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                  'Could not get description of file "'
     &                  // XGRID( 1:TRIMLEN( XGRID ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NROWS .OR.
     &          NCOLS3D .NE. NCOLS      ) THEN
 
            WRITE( MESG,94020 )
     &       'LAYGRD grid dimensions inconsistent with program' //
     &       CRLF() // BLANK5 // '(NCOLS,NROWS) file: (', NCOLS3D, ',',
     &       NROWS3D, ') program: (', NCOLS, ',', NROWS, ').'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        ELSEIF( NLAYS3D .LT. EMLAYS ) THEN
            WRITE( MESG,94020 )
     &       'LAYGRD layers less than number of emissions layers!' //
     &       CRLF() // BLANK5 // 'File layers =', NLAYS3D,
     &       ', but program (EMLAYS)=', EMLAYS
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        END IF

C.........  Get surface met files header for scenario and cloud scheme
C.........  and check dimensions
        IF ( .NOT. DESC3( SNAME ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                  'Could not get description of file "'
     &                  // SNAME( 1:TRIMLEN( SNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NROWS .OR.
     &          NCOLS3D .NE. NCOLS      ) THEN
 
            WRITE( MESG,94020 )
     &       'SFCCRO grid dimensions inconsistent with program' //
     &       CRLF() // BLANK5 // '(NCOLS,NROWS) file: (', NCOLS3D, ',',
     &       NROWS3D, ') program: (', NCOLS, ',', NROWS, ').'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        END IF

        SCEN_STR = FDESC3D( 2 )
        SCHE_STR = FDESC3D( 3 )

C.......   Check file grid with program grid
        IF( NCOLS3D * NROWS3D .NE. NGRID ) THEN

            WRITE( MESG,94010 ) 
     &       'Grid cell inconsistency. Met files: ', NCOLS3D * NROWS3D, 
     &       ', but dimensioned (NGRID)=', NGRID
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ELSEIF( GDNAM3D .NE. GRDNM ) THEN

            WRITE( MESG,94010 ) 
     &       'Grid name inconsistency. Met files: "' // GDNAM3D //
     &       '", but compiled is "' // GRDNM // '"'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ENDIF

C.......   Get time period for output file, dot point grid desc parameters:

        IF ( .NOT. DESC3( DNAME ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                  'Could not get description of file "' //
     &                  DNAME( 1:TRIMLEN( DNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NROWS+1 .OR.
     &          NCOLS3D .NE. NCOLS+1      ) THEN
 
            WRITE( MESG,94020 )
     &       'METDOT grid dimensions inconsistent with program' //
     &       CRLF() // BLANK5 // '(NCOLS,NROWS) file: (', NCOLS3D, ',',
     &       NROWS3D, ') program: (', NCOLS+1, ',', NROWS+1, ').'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        ELSE IF ( VGTYP3D .NE. VGSGPN3 ) THEN

            WRITE( MESG,94010 ) 
     &          'VG type', VGTYP3D, 
     &          'in file "' // DNAME( 1:TRIMLEN( DNAME ) ) // 
     &          '" is _not_ nonhydrostatic sigma-P'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        END IF

C.......   Get default time characteristic for output file:
C.......   If we're going to prompt, then set the defaults from the met data
C.......      otherwise, use environment variables to set defaults
        EFLAG = .FALSE.
        IF( PROMPTF ) THEN

            JDATE  = SDATE3D
            JTIME  = STIME3D
            NSTEPS = MXREC3D

        ELSE

            JDATE  = ENVINT( 'ESDATE', 'Start date', 1995182, IOS )
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'Env var PROMPTFLAG is N, but ESDATE is bad'
                CALL M3MSG2( MESG )
            ENDIF

            JTIME  = ENVINT( 'STIME' , 'Start time', 0 , IOS )
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'Env var PROMPTFLAG is N, but STIME is bad'
                CALL M3MSG2( MESG )
            ENDIF
 
            NSTEPS = ENVINT( 'DURATION', 'Duration', 24, IOS )
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'Env var PROMPTFLAG is N, but DURATION is bad'
                CALL M3MSG2( MESG )
            ENDIF

            IF( EFLAG ) THEN
                CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                       'Bad environment variable setting(s)', 2 )
            ENDIF

        ENDIF

        JDATE  = GETDATE( JDATE, 
     &           'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

        JTIME  = GETNUM( 0, 235959, JTIME, 
     &                  'Enter simulation starting time (HHMMSS)' )
        TSTEP  = 10000
        NSTEPS = GETNUM( 1, 999999, NSTEPS,
     &                   'Enter output duration (hours)' )

        DTBUF = MMDDYY( JDATE )
        WRITE( MESG,94050 )
     &  'Output Start Date:', DTBUF( 1:TRIMLEN( DTBUF ) ) //
     &                                         CRLF() // BLANK5 //
     &  '       Start Time:', JTIME,'HHMMSS'// CRLF() // BLANK5 //
     &  '       Time Step :', 1    ,'hour'  // CRLF() // BLANK5 //
     &  '       Duration  :', NSTEPS, 'hours'

        CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

C.........  Set up output file:
        NCOLSDG = NCOLS3D
        NROWSDG = NROWS3D
        NLAYSDG = NLAYS3D
        P_ALPDG = P_ALP3D
        P_BETDG = P_BET3D
        P_GAMDG = P_GAM3D
        XCENTDG = XCENT3D
        YCENTDG = YCENT3D
        XORIGDG = XORIG3D
        YORIGDG = YORIG3D
        XCELLDG = XCELL3D
        YCELLDG = YCELL3D
        GRDNMDG = GDNAM3D
        VGTYPDG = VGTYP3D
        VGTOPDG = VGTOP3D
        VGLVSDG( 0 ) =       VGLVS3D( 0 )
        SIGH   ( 0 ) = 0.5 * VGLVS3D( 0 )
        DO  11  I = 1, NLAYS3D - 1
            VGLVSDG( I ) = VGLVS3D( I )
            SIGH   ( I ) = 0.5 * ( VGLVS3D( I ) + VGLVS3D( I-1 ) )
11      CONTINUE

        IF ( .NOT. DSCGRID( GRDNMDG, COONMDG, GDTYP3D, 
     &                 P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                 XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                 NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get description for grid "' //
     &              GRDNMDG( 1: TRIMLEN( GRDNMDG ) ) //
     &              '"; exiting program'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ELSE IF ( VGTYP3D .NE. VGSGPN3 ) THEN

            WRITE( 94010,MESG ) 
     &          'Unsupported VGTYPE=', VGTYP3D, 
     &          'in "' // DNAME( 1:TRIMLEN( DNAME ) ) //
     &          '"; need VGTYPE=VGSGPN3=2'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        END IF


C.......   Get time period for output file, dot point grid desc parameters:

        IF ( .NOT. DESC3( XNAME ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // XNAME( 1:TRIMLEN( XNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NROWS .OR.
     &          NCOLS3D .NE. NCOLS      ) THEN
 
            WRITE( MESG,94020 )
     &       'LAYCRO grid dimensions inconsistent with program' //
     &       CRLF() // BLANK5 // '(NCOLS,NROWS) file: (', NCOLS3D, ',',
     &       NROWS3D, ') program: (', NCOLS, ',', NROWS, ').'
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        ELSEIF( NLAYS3D .LT. EMLAYS ) THEN

            WRITE( MESG,94020 )
     &       'LAYCRO layers less than number of emissions layers!' //
     &       CRLF() // BLANK5 // 'File layers =', NLAYS3D,
     &       ', but program (EMLAYS)=', EMLAYS
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Bad input file.', 2 )

        END IF

        NCOLSXG = NCOLS3D
        NROWSXG = NROWS3D
        NLAYSXG = NLAYS3D
        P_ALPXG = P_ALP3D
        P_BETXG = P_BET3D
        P_GAMXG = P_GAM3D
        XCENTXG = XCENT3D
        YCENTXG = YCENT3D
        XORIGXG = XORIG3D
        YORIGXG = YORIG3D
        XCELLXG = XCELL3D
        YCELLXG = YCELL3D
        GRDNMXG = GDNAM3D
        VGTYPXG = VGTYP3D
        VGTOPXG = VGTOP3D
        DO  12  I = 0, NLAYS3D - 1
            VGLVSXG( I ) = VGLVS3D( I )
12      CONTINUE

        IF ( .NOT. DSCGRID( GRDNMXG, COONMXG, GDTYP3D, 
     &                 P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                 XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                 NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get description for grid "' //
     &              GRDNMXG( 1: TRIMLEN( GRDNMXG ) ) //
     &              '"; exiting program'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ELSE IF ( COONMXG .NE. COONMDG ) THEN

            MESG = 'Dot and cross MET coordinate system names "' //
     &              COONMDG( 1:TRIMLEN( COONMDG ) ) // '" and "' //
     &              COONMXG( 1:TRIMLEN( COONMXG ) ) // '" do not match'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        ELSE IF ( DBLERR( XCELLXG, XCELLDG ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross XCELL mismatch', 2 )

        ELSE IF ( DBLERR( YCELLXG, YCELLDG ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross YCELL mismatch', 2 )

        ELSE IF ( DBLERR( XORIGXG, XORIGDG + 0.5D0 * XCELLDG ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross XORIG mismatch', 2 )

        ELSE IF ( DBLERR( YORIGXG, YORIGDG + 0.5D0 * YCELLDG ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross YORIG mismatch', 2 )

        ELSE IF ( VGTYPXG .NE. VGTYPDG ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross VG TYPE mismatch', 2 )

        ELSE IF ( FLTERR( VGTOPXG, VGTOPDG ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Dot and cross VGTOP mismatch', 2 )

        ELSE

            DO  13  I = 0, NLAYS3D - 1
                IF ( FLTERR(  VGLVSXG( I ), VGLVSDG( I ) ) ) THEN
                    CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                           'Dot and cross VGLVS mismatch', 2 )
                END IF
13          CONTINUE

        END IF

C.......   Store begining and ending coordinates of grid in X and Y
        XBEG = XORIG3D
        YBEG = YORIG3D
        XEND = XORIG3D + NCOLS * XCELL3D
        YEND = YORIG3D + NROWS * YCELL3D

C.......   Get description of met input file; use most of it 
C.......   for description of output file

        IF ( .NOT. DESC3( ENAME ) ) THEN

            MESG = 'Could not get description of file "' // 
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )
        ELSE IF ( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  EMISSIONS file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )
        END IF

        IF ( GDNAM3D .NE. COONMDG ) THEN

            MESG = 'EMIS and MET coordinate system names "' //
     &              GDNAM3D( 1: TRIMLEN( GDNAM3D ) ) // '" and "' //
     &              COONMDG( 1: TRIMLEN( COONMDG ) ) // '" do not match'
            CALL M3WARN( 'LAYPOINT', 0, 0, MESG )

        END IF

        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = 10000         !  1 hour, no minutes, no seconds
        VGTYP3D = VGTYPXG
        VGTOP3D = VGTOPXG
        DO  21  I = 0, EMLAYS 
            VGLVS3D( I ) = VGLVSXG( I )
21      CONTINUE

        NLAYS3D = EMLAYS
        NVARS3D = 1
        VNAME3D( 1 ) = 'LFRAC'
        VTYPE3D( 1 ) = M3REAL
        UNITS3D( 1 ) = 'none'
        VDESC3D( 1 ) = 
     &  'fraction of point source plume emitted into each model layer'

        FDESC3D( 1 ) =
     &  'Source level hourly plume rise layer fractions file'
        DO  22  I = 2, MXDESC3
            FDESC3D( I ) = ' '
22      CONTINUE

        FDESC3D( 2 ) = SCEN_STR
        FDESC3D( 3 ) = SCHE_STR

        LNAME = PROMPTMFILE( 
     &          'Enter name for LAYER FRACTIONS output file',
     &          FSUNKN3, 'PLAY', 'LAYPOINT' )

        IF( REP_LAYR .NE. 9999 ) THEN
            RDEV = PROMPTFFILE( 'Enter logical name for report',
     &                        .FALSE., .TRUE., 'REPRTLAY', 'LAYPOINT' )
            WRITE( RDEV,93040 ) REP_LAYR

        ENDIF

C...........   Read in emissions stack parameters IFIP, ISCC, XLOCA, YLOCA,
C...........   STKHT, STKDM, STKTK, STKVE

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading IFIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading ISCC from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading IPLANT from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading ISTACK from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3,0,0, XLOCA ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading XLOCA from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3,0,0, YLOCA ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading YLOCA from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading STKHT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading STKDM from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading STKTK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                   'Error reading STKVE from file' // ENAME, 2 )
        END IF


C...........   Compute un-gridding matrices for dot and cross point
C...........   met data:

        CALL UNGRIDB( NCOLSDG, NROWSDG, 
     &                XORIGDG, YORIGDG, XCELLDG, YCELLDG,
     &                NPSRC, XLOCA, YLOCA, ND, CD )

        CALL UNGRIDB( NCOLSXG, NROWSXG, 
     &                XORIGXG, YORIGXG, XCELLXG, YCELLXG,
     &                NPSRC, XLOCA, YLOCA, NX, CX )


C.......   Read time-independent ZF, ZH:

        IF ( .NOT. READ3( XGRID, 'ZH', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            MESG = 'Could not read ZH from "' //
     &             XGRID( 1 : TRIMLEN( XGRID ) ) // '".'
            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

        END IF      !  if read3() failed on xgrid
        CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, ZH )

        IF ( .NOT. READ3( XGRID, 'ZF', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            MESG = 'Could not read ZF from "' //
     &             XGRID( 1 : TRIMLEN( XGRID ) ) // '".'
            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

        END IF      !  if read3() failed on xgrid
        CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, ZF )

        DO  45  S = 1, NPSRC

            ZZ0 = ZF( 1,S )
            ZSTK ( 1,S ) = ZZ0 - STKHT( S )
            ZF0 = ZF( 1,S )
            DDZF( 1,S ) = 1.0 / ZF0

            DO  44  L = 2, EMLAYS

                ZZ1 = ZF( L,S )
                ZSTK( L  ,S ) = ZZ1 - STKHT( S )
                DDZH( L-1,S ) = 1.0 / ( ZZ1 - ZZ0 )
                ZZ0 = ZZ1
                ZF1 = ZF( L,S )
                DDZF( L,S ) = 1.0 / ( ZF1 - ZF0 )
                ZF0 = ZF1

44          CONTINUE

45      CONTINUE


C...........   For each time step, compute the layer fractions:


        CALL M3MSG2( 'Calculating hourly layer fractions...' )

        LDATE = 0
        EFLAG = .FALSE.
        DO  144  I = 1, NSTEPS

            IF ( LDATE .NE. JDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

                IF( RDEV .GT. 0 ) THEN
                    WRITE( RDEV,93000 ) MESG
                ENDIF

                LDATE = JDATE
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

C.............  Write to report file if report feature is on
            IF( RDEV .GT. 0 ) THEN
                WRITE( RDEV,93020 ) HHMMSS( JTIME )
            ENDIF

C...........   Read and transform meteorology:

            IF ( .NOT. READ3( SNAME, 'HFX', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read HFX from "' // 
     &               SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, HFX )

            IF ( .NOT. READ3( SNAME, 'PBL', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read PBL from "' // 
     &               SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, HMIX )

            IF ( .NOT. READ3( SNAME, 'TGD', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read TGD from "' // 
     &               SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, TSFC )

            IF ( .NOT. READ3( SNAME, 'USTAR', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read USTAR from "' // 
     &               SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, USTAR )

            IF ( .NOT. READ3( XNAME, 'TA', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read TA from "' // 
     &               XNAME( 1 : TRIMLEN( XNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, TA )

            IF ( .NOT. READ3( XNAME, 'TV', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read TV from "' // 
     &               XNAME( 1 : TRIMLEN( XNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, TV )

            IF ( .NOT. READ3( XNAME, 'PRES', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read PRES from "' // 
     &               XNAME( 1 : TRIMLEN( XNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, PRES )

            IF ( .NOT. READ3( DNAME, 'UWIND', ALLAYS3,
     &                        JDATE, JTIME, DBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read UWIND from "' // 
     &               DNAME( 1 : TRIMLEN( DNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NDOTS, NPSRC, EMLAYS, ND, CD, DBUF, UWIND )

            IF ( .NOT. READ3( DNAME, 'VWIND', ALLAYS3,
     &                        JDATE, JTIME, DBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read VWIND from "' // 
     &               DNAME( 1 : TRIMLEN( DNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NDOTS, NPSRC, EMLAYS, ND, CD, DBUF, VWIND )


            P  =       1.0  / ( SIGH( 1 ) - SIGH( 0 ) )
            PP =     - P * ( VGLVSXG( 0 ) - SIGH( 0 ) )
            QQ = 2.0 - P * ( VGLVSXG( 1 ) - SIGH( 1 ) ) 

            DO 133  S = 1, NPSRC

C...........   Skip source if it is outside grid

	       XL = XLOCA( S )
	       YL = YLOCA( S )
               IF( XL .LT. XBEG .OR. XL .GT. XEND .OR.
     &             YL .LT. YBEG .OR. YL .GT. YEND     ) THEN 
                   GO TO 133

	       ENDIF

C...........   Convert pressure to millibars from pascals; compute wind speed:

                DO 101  L = 1, EMLAYS

                    PRES( L,S ) = 1.0E-2 * PRES( L,S )

                    P = UWIND( L,S )
                    Q = VWIND( L,S )
                    WSPD( L ) = SQRT( P * P  +  Q * Q )

101             CONTINUE

C...........   Compute other derived meteorology variables:

                X     = PRES( 1,S )
                Y     = PRES( 2,S )

                PSFC  = PP * Y + QQ * X
                ES    = 6.1078 * EXP( 5384.21 / CTOK - 
     &                                5384.21 / TSFC( S ) )
                QSFC  = 0.622 * ES / ( PSFC - ES )
                TVSFC = TSFC( S ) * ( 1.0 + 0.6077 * QSFC )
                THETG = TVSFC * ( 1000.0 / PSFC )**0.286
                THV1  = TV( 1,S )*( 1000.0 / PRES( 1,S ) )**0.286
                DTHDZ( 1 ) = ( THV1 - THETG ) / ZH( 1,S )
                IF ( HMIX ( S ) .LE. ZF( 1,S ) ) THEN
                    LPBL = 1
                END IF
                IF ( STKHT( S ) .LE. ZF( 1,S ) ) THEN
                    LSTK = 1
                END IF
                ZZF( 0 ) = 0.0
                ZZF( 1 ) = ZF( 1,S )

                DO 111  L = 2, EMLAYS

                    IF ( HMIX ( S ) .GT. ZF( L-1,S ) )  LPBL = L
                    IF ( STKHT( S ) .GT. ZF( L-1,S ) )  LSTK = L

                    THVK = TV( L,S ) * ( 1000.0 / PRES( L,S ) )**0.286
                    DTHDZ( L ) = DDZH( L-1,S ) * ( THVK - THV1 )
                    THV1 = THVK

                    ZZF( L ) = ZF( L,S )

111             CONTINUE

                M    = MAX( 1, LSTK - 2 )
                TSTK = POLY( ZSTK( LSTK,S ),
     &                       ZH( M,S ),
     &                       TA( M,S ), 3 )
                WSTK = MAX( POLY( ZSTK( LSTK,S ),
     &                            ZH( M,S ),
     &                            WSPD( M ), 3 ),
     &                      0.1 )

C...............   Compute plume rise for this source (noting inversion of MOL):
C...............   Trap USTAR at a minimum realistic value 

                CALL PLMRIS( ZH( 1,S ), ZZF( 0 ), WSPD, 
     &               TA( 1,S ), DTHDZ, MAX( USTAR( S ), USTARMIN ), 
     &               HMIX( S ), HFX( S ), STKHT( S ), STKDM( S ), 
     &               STKTK( S ), STKVE( S ), ZSTK( 1,S ), TSTK, 
     &               WSTK, LSTK, LPBL, ZTOP, ZBOT )

C...............   Check plume rise for nonsense output
                IF( ZTOP .LT. STKHT( S ) ) THEN
                    WRITE( MESG,94010 ) 
     &                     'WARNING: Top of plume found to be ' //
     &                     '< top of stack!' // CRLF() // BLANK5 //
     &                     'for source: ', S, '( FIP:', IFIP( S ), 
     &                     'SCC:', ISCC( S ), 'PLT:', IPLANT( S ), 
     &                      'STK:', ISTACK( S ), ' )'
                    CALL M3MESG( MESG )

                ENDIF

C...............   Compute plume fractions, assuming uniform distribution
C...............   in pressure (~mass concentration -- minor hydrostatic
C...............   assumption) from bottom to top.  

C...................   Compute LBOT, LTOP so that 
C...................   ZZF( LBOT-1 ) <= ZBOT < ZZF( LBOT ) and
C...................   ZZF( LTOP-1 ) <= ZTOP < ZZF( LTOP )

                DO  121  L = 1, EMLAYS - 1
                    IF ( ZBOT .LE. ZZF( L ) ) THEN
                        LBOT = L
                        GO TO  122
                    ELSE
                        LFRAC( S,L ) = 0.0	! fractions below plume
                    END IF
121             CONTINUE
                LBOT = EMLAYS           !  fallback
122             CONTINUE                !  loop exit:  bottom found at LBOT

                IF ( ZTOP .LE. ZZF( LBOT ) ) THEN  !  plume in this layer

                    LFRAC( S,LBOT ) = 1.0
                    LTOP = LBOT

                    DO  123  L = LBOT + 1, EMLAYS  ! fractions above plume
                        LFRAC( S,L ) = 0.0
123                 CONTINUE

C.................  Note: this check not in original algorithm, but w/o it,
C                         can end up with fractions > 1.0
                ELSEIF( LBOT .EQ. EMLAYS ) THEN    ! plume above top layer

                    LFRAC( S,LBOT ) = 1.0

                    DO  124  L = 1, EMLAYS-1       ! fractions below plume
                        LFRAC( S,L ) = 0.0
124                 CONTINUE

                ELSE                               ! plume crosses layers

                    DO  125  L = LBOT + 1, EMLAYS
                        IF ( ZTOP .LE. ZZF( L ) ) THEN
                            LTOP = L
                            GO TO 126
                        END IF
125                 CONTINUE
                    LTOP = EMLAYS
126                 CONTINUE

C...................   Compute corresponding PBOT,PTOP so that
C...................   PRESF( LBOT-1 ) <= PBOT < PRESF( LBOT ) and
C...................   PRESF( LTOP-1 ) <= PTOP < PRESF( LTOP )
C...................   (Use 3rd order polynomial via POLY() )

                    PRESF( 0 ) = PSFC
                    PRESF( 1 ) = POLY( ZZF ( 1 ), 
     &                                 ZH  ( 1,S ), 
     &                                 PRES( 1,S ), 3 )
                    PRESF( 2 ) = POLY( ZZF ( 2 ), 
     &                                 ZH  ( 1,S ), 
     &                                 PRES( 1,S ), 3 )

                    DO  127  L = 3, EMLAYS-2
                        M = L - 2
                        PRESF( L ) = POLY( ZZF ( L ), 
     &                                     ZH  ( M,S ), 
     &                                     PRES( M,S ), 3 )
127                 CONTINUE

                    PRESF( EMLAYS-1 ) = POLY( ZZF ( EMLAYS-1 ),
     &                                        ZH  ( EMLAYS-3,S ),
     &                                        PRES( EMLAYS-3,S ), 3 )
                    PRESF( EMLAYS   ) = POLY( ZZF ( EMLAYS ),
     &                                        ZH  ( EMLAYS-1,S ),
     &                                        PRES( EMLAYS-1,S ), 1 )

                    M    = MIN( MAX( 0, LBOT-2 ), EMLAYS-3 )
                    PBOT = POLY( ZBOT, ZZF( M ), PRESF( M ), 3 )

                    M    = MIN( MAX( 0, LTOP-2 ), EMLAYS-3 )
                    PTOP = POLY( ZTOP, ZZF( M ), PRESF( M ), 3 )

                    PDIFF = PBOT - PTOP
                    IF( PDIFF .GT. 0 ) THEN
                        DDP  = 1.0 / ( PDIFF )  !  = d(plumefrac)/d(pres)

                    ELSE
                        WRITE( LDEV,93030 ) 
     &                         ZBOT, PBOT, ZTOP, PTOP, PDIFF, VGTOPXG
                        WRITE( MESG,94010 ) 
     &                   'Infinitely small plume created for source ', S
                        CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

                    ENDIF

                    LFRAC( S,LBOT ) = DDP * ( PBOT - PRESF( LBOT ) )
                    LFRAC( S,LTOP ) = DDP * ( PRESF( LTOP-1 ) - PTOP )

                    DO  128  L = LBOT+1, LTOP-1	!  layers in plume
                        LFRAC( S,L ) = DDP*( PRESF( L-1 ) - PRESF( L ) )
128                 CONTINUE

                    DO  129  L = LTOP+1, EMLAYS	!  fractions above plume
                        LFRAC( S,L ) = 0.0
129                 CONTINUE

                END IF          !  if ztop in same layer as zbot, or not

C.................  Check if LTOP out of range, and report (will only work
C.................    if REP_LAYR env var has been set b/c default is 9999
                IF( LTOP .GT. REP_LAYR ) THEN

                    WRITE( RDEV,93042 ) S, IFIP( S ), ISCC( S ), 
     &                     IPLANT( S ), ISTACK( S ), STKHT( S ),
     &                     STKVE ( S ), STKTK( S ), TSTK, WSTK, LPBL, 
     &                     LTOP 

                ENDIF

133         CONTINUE	!  end loop on sources S

            IF( EFLAG ) THEN

                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &                       'Bad plume rise calculations.' , 2 )

            ELSEIF ( .NOT. WRITE3( LNAME, 'LFRAC', 
     &                         JDATE, JTIME, LFRAC ) ) THEN

                MESG = 'Error writing LFRAC to ' // LNAME
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )

144     CONTINUE	!  end loop on time steps I

999     CONTINUE          !  exit program

        MESG = 'Normal completion of LAYPOINT'
        CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, I10 )

92011   FORMAT (/5X , A, I10 )

92020   FORMAT ( 5X , A, I9 )

92030   FORMAT ( 5X , A, I12, :, 2X, A  )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )

93030   FORMAT( 8X, 'ERROR computing plume rise!' /
     &         10X, 'Plume bottom at ',F7.2,' [m] and ',F7.2, ' sigma' /
     &         10X, 'Plume top    at ',F7.2,' [m] and ',F7.2, ' sigma' /
     &         10X, 'Difference      ',F7.2,' [m]' /
     &         10X, 'Domain top   at ',F7.2, '[cb]' )

93040   FORMAT( 'Sources with top of plume greater than layer', I3, //,
     &          'Src ID', 1X, 'FIP', 3X, 'SCC', 6X, 'Plant' , 2X,
     &          'Stack', 2X, 'H[m]', 3X, 'V[m/s]', 1X, 'Ts[K]', 2X, 
     &          'Ta[K]', 1X, 'U[m/s]', 1X, 'LPBL', 1X, 'LTOP' )

93042   FORMAT( I6, 1X, I5.5, 1X, I8.8, 1X, I6, 1X, I6, 1X, F6.1, 1X,
     &          F6.2, 1X, F6.1, 1X, F5.1, 1X, F6.2, 1X, I3, 1X, I3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 12( A, :, I8, :, 1X ) )

94020   FORMAT( 10 ( A, :, I4, :, 1X ) )

94050   FORMAT( A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )


        END

