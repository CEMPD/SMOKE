
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM LAYPOINT

C***********************************************************************
C  program body starts at line 212
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
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'LADIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'GRDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCGRID
        REAL		ENVREAL
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

      EXTERNAL       DSCGRID, ENVREAL, GETNUM, GETYN, 
     &               PROMPTFFILE, PROMPTMFILE, TRIMLEN
     
                                            
C...........   PARAMETERS and their descriptions

        INTEGER     NDOTS
        REAL        GAMA
        PARAMETER ( NDOTS = ( NCOLS + 1 ) * ( NROWS + 1 ) ,
     &              GAMA  = -0.0098 )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Point source stack parameters:

        INTEGER     IFIP ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC ( NPSRC )  !  source SCC
        REAL        XLOCA( NPSRC )  !  UTM X-location (m)
        REAL        YLOCA( NPSRC )  !  UTM Y-location (m)
        REAL        STKHT( NPSRC )  !  stack height   (m)
        REAL        STKDM( NPSRC )  !  stack diameter (m)
        REAL        STKTK( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE( NPSRC )  !  exhaust velocity    (m/s)

C.......   Point source meteorology profiles and variables

        REAL        PSTR0( NPSRC )  !  reference P*
        REAL        TSFC ( NPSRC )  !  surface temperature (deg K)
        REAL        HFX  ( NPSRC )  !  sensible heat flux (M K / S )
        REAL        HMIX ( NPSRC )  !  mixing height (m)
        REAL        ZRUF ( NPSRC )  !  surface roughness (m)
        REAL        USTAR( NPSRC )  !  friction velocity (m/s)
        REAL        MOL  ( NPSRC )  !  inverse Monin-Obukhov length (1/m)
        REAL        PRES0( EMLAYS, NPSRC )  !  reference pressure
        REAL        ZF   ( EMLAYS, NPSRC )  !  layer surface height (m)
        REAL        ZH   ( EMLAYS, NPSRC )  !  layer center  height (m)
        REAL        TA   ( EMLAYS, NPSRC )  !  temperature (K)
        REAL        TV   ( EMLAYS, NPSRC )  !  virtual temperature
        REAL        QV   ( EMLAYS, NPSRC )  !  mixing ratio
        REAL        UWIND( EMLAYS, NPSRC )  !  wind speed (m/s)
        REAL        VWIND( EMLAYS, NPSRC )  !  wind speed (m/s)
        REAL        PRES ( EMLAYS, NPSRC )  !  pressure (Pa)
        REAL        WSPD ( NLAYS )  !  wind speed (m/s)
        REAL        PSFC            !  surface pressure (Pa)
        REAL        QSFC            !  surface heat flux
        REAL        ES              !  water vapor pressure
        REAL        TVSFC           !  sfc virtual temperature
        REAL        THETG           !  sfc virtual potential temperature
        REAL        THV1, THVK      !  virtual potential temperature
        INTEGER     LSTK            !  first L: ZH(L) > STKHT
        INTEGER     LPBL            !  first L: ZH(L) above mixing layer
        REAL        DTHDZ( NLAYS )  !  gradient of THETV
        REAL        DDZH ( NLAYS,NPSRC )  !  1/( zh(l) - zh(l-1) )
        REAL        ZSTK ( NLAYS,NPSRC )  !  zh( l,s )   - stkht(s)   (m)
        REAL        TSTK            !  temperature at top of stack (deg K)
        REAL        WSTK            !  wind speed  at top of stack (m/s)

C.......   Cross- and dot-point meteorology input buffers

        REAL	    XBUF( NGRID, NLAYS )
        REAL	    DBUF( NDOTS, NLAYS )

C.......   Dot point un-gridding matrix (uses bilinear interpolation)

        INTEGER         ND( 4,NPSRC )  !  cell indexes
        REAL            CD( 4,NPSRC )  !  coefficients

C.......   Cross point un-gridding matrix

        INTEGER         NX( 4,NPSRC )  !  cell indexes
        REAL            CX( 4,NPSRC )  !  coefficients

C.......   Output layer fractions (in I/O order)

        REAL        LFRAC( NPSRC, EMLAYS )

C.......   Empirically developed plume spread factors:

        REAL    SZ0FAC          !  ???
        REAL    SPRFAC          !  ???

C...........   Other local variables

        INTEGER         TSTEP   !  output time step
        INTEGER         JDATE, JTIME
        INTEGER         NSTEPS
        INTEGER         I, L, S
        REAL            X, Y, P, Q, PP, QQ
        REAL            ZZ0, ZZ1
        REAL            BOT 	!  plume bottom elevation (m)
        REAL            TOP 	!  plume top    elevation (m)
        REAL            DDP  	!  1/(plume-thickness)
        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    DNAME   !  dot-point   met file name
        CHARACTER*16    SNAME   !  cross-point sfc met file name
        CHARACTER*16    XNAME   !  cross-point lay met file name
        CHARACTER*16    SGRID   !  cross-point sfc grid file name
        CHARACTER*16    XGRID   !  cross-point lay grid file name
        CHARACTER*16    LNAME   !  " timestepped layer fractions output file
        CHARACTER*256   MESG    !  buffer for M3EXIT() messages
        INTEGER		FDEV    !  unit number for stack parameters file
        INTEGER		LDEV    !  unit number for log file
        INTEGER		ISTAT   !  ENVREAL() status

C.......   Dot and cross point grid description parameters:

        INTEGER      NCOLSDG      ! dot grid number of grid columns
        INTEGER      NROWSDG      ! dot grid number of grid rows
        INTEGER      NLAYSDG      ! dot grid number of layers
        REAL*8       P_ALPDG      ! dot grid first, second, third map
        REAL*8       P_BETDG      ! dot grid projection descriptive
        REAL*8       P_GAMDG      ! dot grid parameters.
        REAL*8       XCENTDG      ! dot grid lon for coord-system X=0
        REAL*8       YCENTDG      ! dot grid lat for coord-system Y=0
        REAL*8       XORIGDG      ! dot grid X-coordinate origin of grid (map units)
        REAL*8       YORIGDG      ! dot grid Y-coordinate origin of grid
        REAL*8       XCELLDG      ! dot grid X-coordinate cell dimension
        REAL*8       YCELLDG      ! dot grid Y-coordinate cell dimension
        INTEGER      VGTYPDG      ! vertical coordinate type (VGSIGP3, ...)
        REAL         VGTOPDG      ! model-top, for sigma coord types.
        REAL         VGLVSDG( MXLAYS3 + 1 )  !  vertical coord values.
        CHARACTER*16 GRDNMDG      ! dot grid grid name
        CHARACTER*16 COONMDG      ! dot grid coordinate system name

        INTEGER      NCOLSXG      ! cross grid number of grid columns
        INTEGER      NROWSXG      ! cross grid number of grid rows
        INTEGER      NLAYSXG      ! cross grid number of layers
        REAL*8       P_ALPXG      ! cross grid first, second, third map
        REAL*8       P_BETXG      ! cross grid projection descriptive
        REAL*8       P_GAMXG      ! cross grid parameters.
        REAL*8       XCENTXG      ! cross grid lon for coord-system X=0
        REAL*8       YCENTXG      ! cross grid lat for coord-system Y=0
        REAL*8       XORIGXG      ! cross grid X-coordinate origin of grid (map units)
        REAL*8       YORIGXG      ! cross grid Y-coordinate origin of grid
        REAL*8       XCELLXG      ! cross grid X-coordinate cell dimension
        REAL*8       YCELLXG      ! cross grid Y-coordinate cell dimension
        CHARACTER*16 GRDNMXG      ! cross grid grid name
        CHARACTER*16 COONMXG      ! cross grid coordinate system name
        INTEGER      VGTYPXG      ! vertical coordinate type (VGSIGP3, ...)
        REAL         VGTOPXG      ! model-top, for sigma coord types.
        REAL         VGLVSXG( MXLAYS3 + 1 )  !  vertical coord values.
        REAL         SIGH   ( MXLAYS3 )      !  half-level sigma values

C.......   STATEMENT FUNCTION for evaluating relative error: TRUE
C.......   iff (PP,QQ) are significantly different (compares square of
C.......   normalized difference against square of tolerance)
 
        LOGICAL         DBLERR
        LOGICAL         FLTERR
        REAL*8          PP,QQ
        REAL            RR, SS
 
        DBLERR( PP,QQ ) = 
     &    ( ( PP - QQ )**2 .GT. 1.0D-12*( PP*PP + QQ*QQ + 1.0D-5 ) ) 

        FLTERR( RR,SS ) = 
     &    ( ( RR - SS )**2 .GT. 1.0E-12*( RR*PP + SS*SS + 1.0E-5 ) ) 


C***********************************************************************
C   begin body of program LAYPOINT

        LDEV = INIT3()

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

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 'Ending program', 2 )
        END IF
        
C.......   Get empirical plume rise/spread constants:

        SZ0FAC = ENVREAL( 'SZ0FAC', 
     &                    'Empirical plume spread factor SZ0FAC',
     &                    2.0, ISTAT )
        SPRFAC = ENVREAL( 'SPRFAC',
     &                    'Empirical plume spread factor SPRFAC',
     &                    10.0, ISTAT )
        

C.......   Get file name; open input point sources, temporal cross-reference,
C.......   and temporal profiles files

        ENAME = PROMPTMFILE( 
     &          'Enter name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'LAYPOINT' )

        FDEV = PROMPTFFILE( 'Enter name for STACK PARAMETER input file',
     &                      .TRUE., .TRUE., 'PSTK', 'LAYPOINT' )
        IF ( FDEV .LE. 0 ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &           'Could not open input STACK PARAMETER FILE', 2 )
        END IF

        DNAME = PROMPTMFILE( 
     &          'Enter name for DOT-POINT   LAYERED  MET file',
     &          FSREAD3, 'METDOT', 'LAYPOINT' )

        SNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT SURFACE  MET file',
     &          FSREAD3, 'SFCCROSS', 'LAYPOINT' )

        XNAME = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED  MET file',
     &          FSREAD3, 'LAYCROSS', 'LAYPOINT' )

        SGRID = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT SURFACE GRID file',
     &          FSREAD3, 'SFCCROSS', 'LAYPOINT' )

        XGRID = PROMPTMFILE( 
     &          'Enter name for CROSS-POINT LAYERED GRID file',
     &          FSREAD3, 'LAYCROSS', 'LAYPOINT' )


C.......   Get time period for output file, dot point grid desc parameters:

        IF ( .NOT. DESC3( DNAME ) ) THEN

            CALL M3EXIT( 'LAYPOINT', 0, 0, 
     &                  'Could not get description of file "' //
     &                  DNAME( 1:TRIMLEN( DNAME ) ) // '"', 2 )

        ELSE IF ( VGTYP3D .NE. VGSGPN3 ) THEN

            WRITE( MESG,94010 ) 
     &          'VG type', VGTYP3D, 
     &          'in file "' // DNAME( 1:TRIMLEN( DNAME ) ) // 
     &          '" is _not_ nonhydrostatic sigma-P'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        END IF

        JDATE  = SDATE3D
        JTIME  = STIME3D
        NSTEPS = 24 * MXREC3D

        JDATE  = GETNUM( 0, 9999999, JDATE, 
     &                  'Enter starting date (YYYYDDD)' )
        JTIME  = GETNUM( 0, 999999, JTIME, 
     &                  'Enter starting time  (HHMMSS)' )
        TSTEP  = 10000
        NSTEPS = GETNUM( 1, 999999, NSTEPS,
     &                   'Enter output duration (hours)' )

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
        VGLVSDG( 1 ) =       VGLVS3D( 1 )
        SIGH   ( 1 ) = 0.5 * VGLVS3D( 1 )
        DO  11  I = 2, NLAYS3D
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
     &                  // XNAME( 1:TRIMLEN( DNAME ) ) // '"', 2 )
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
        DO  12  I = 1, NLAYS3D
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

            DO  13  I = 1, NLAYS3D
                IF ( FLTERR(  VGLVSXG( L ), VGLVSDG( L ) ) ) THEN
                    CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                           'Dot and cross VGLVS mismatch', 2 )
                END IF
13          CONTINUE

        END IF


C.......   Get description of met input file; use most of it 
C.......   for description of output file

        IF ( .NOT. DESC3( ENAME ) ) THEN

            MESG = 'Could not get description of file "' // 
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        END IF
        
        IF ( GDNAM3D .NE. COONMDG ) THEN

            MESG = 'EMIS and MET coordinate system names "' //
     &              GDNAM3D( 1: TRIMLEN( GDNAM3D ) ) // '" and "' //
     &              COONMDG( 1: TRIMLEN( COONMDG ) ) // '" do not match'
            CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 2 )

        END IF

        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = 10000         !  1 hour, no minutes, no seconds
        VGTYP3D = VGTYPXG
        VGTOP3D = VGTOPXG
        DO  21  I = 1, EMLAYS + 1
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
        DO  22  I = 2, EMLAYS + 1
            FDESC3D( I ) = ' '
22      CONTINUE
        
        LNAME = PROMPTMFILE( 
     &          'Enter name for LAYER FRACTIONS output file',
     &          FSUNKN3, 'PLAY', 'LAYPOINT' )


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

C...........   Fix up missing stack parameters:

        CALL FIXSTK( FDEV, IFIP, ISCC, STKHT, STKDM, STKTK, STKVE )


C...........   Compute un-gridding matrices for dot and cross point
C...........   met data:

        CALL UNGRIDB( NCOLSDG, NROWSDG, 
     &                XORIGDG, YORIGDG, XCELLDG, YCELLDG,
     &                NPSRC, XLOCA, YLOCA, ND, CD )

        CALL UNGRIDB( NCOLSXG, NROWSXG, 
     &                XORIGXG, YORIGXG, XCELLXG, YCELLXG,
     &                NPSRC, XLOCA, YLOCA, NX, CX )


C.......   Calculation of layer fractions:
C.......   Read time-independent ZZERO, ZF, ZH:
        
        WRITE( *,92000 ) 
     &      ' ', 'Calculating hourly layer fractions...', ' '

        IF ( .NOT. READ3( SGRID, 'ZZERO', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            MESG = 'Could not read ZZERO from "' //
     &             SGRID( 1 : TRIMLEN( SGRID ) ) // '".'
            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

        END IF      !  if read3() failed on dname
        CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, ZRUF )

        IF ( .NOT. READ3( SGRID, 'PSTR0', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            MESG = 'Could not read PSTR0 from "' //
     &             SGRID( 1 : TRIMLEN( SGRID ) ) // '".'
            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )

        END IF      !  if read3() failed on dname
        CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, PSTR0 )

        DO  45  S = 1, NPSRC
  
            ZZ0 = ZH( 1,S )
            PRES0( 1,S ) = SIGH( 1 ) * PSTR0( S )  +  VGTOP3D
            ZSTK ( 1,S ) = ZZ0 - STKHT( S )
            
            DO  44  L = 2, EMLAYS

                PRES0( L,S ) = SIGH( L ) * PSTR0( S )  +  VGTOP3D
                ZZ1 = ZH( L,S )
                ZSTK( L  ,S ) = ZZ1 - STKHT( S )
                DDZH( L-1,S ) = 1.0 / ( ZZ1 - ZZ0 )
                ZZ0 = ZZ1

44          CONTINUE

45      CONTINUE

        IF ( .NOT. READ3( XGRID, 'ZF', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &           'Could not read ZF from "' // 
     &           XGRID( 1 : TRIMLEN( XGRID ) ) // '".', 2 )

        END IF      !  if read3() failed on dname
        CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, ZF )

        IF ( .NOT. READ3( XGRID, 'ZH', ALLAYS3,
     &                    JDATE, JTIME, XBUF ) ) THEN

            CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &           'Could not read ZH from "' // 
     &           XGRID( 1 : TRIMLEN( XGRID ) ) // '".', 2 )

        END IF      !  if read3() failed on dname
        CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, ZH )


C...........   For each time step, compute the layer fractions:


        DO  144  I = 1, NSTEPS

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

            IF ( .NOT. READ3( SNAME, 'MOL', ALLAYS3,	! _inverse_ MOL
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read MOL from "' // 
     &               SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, 1, NX, CX, XBUF, MOL )

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

            IF ( .NOT. READ3( XNAME, 'QV', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read QV from "' // 
     &               XNAME( 1 : TRIMLEN( XNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, QV )

            IF ( .NOT. READ3( XNAME, 'PRES', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read PRES from "' // 
     &               XNAME( 1 : TRIMLEN( XNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NGRID, NPSRC, EMLAYS, NX, CX, XBUF, PRES )

            IF ( .NOT. READ3( DNAME, 'UWIND', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read UWIND from "' // 
     &               DNAME( 1 : TRIMLEN( DNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NDOTS, NPSRC, EMLAYS, ND, CD, DBUF, UWIND )

            IF ( .NOT. READ3( DNAME, 'VWIND', ALLAYS3,
     &                        JDATE, JTIME, XBUF ) ) THEN
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME,
     &               'Could not read VWIND from "' // 
     &               DNAME( 1 : TRIMLEN( DNAME ) ) // '".', 2 )
            END IF      !  if read3() failed on dname
            CALL BMATVEC( NDOTS, NPSRC, EMLAYS, ND, CD, DBUF, VWIND )
            

            P  =       1.0  / ( SIGH( 2 ) - SIGH( 1 ) )
            PP =     - P * ( VGLVSXG( 1 ) - SIGH( 1 ) )
            QQ = 2.0 - P * ( VGLVSXG( 2 ) - SIGH( 2 ) ) 

            DO 133  S = 1, NPSRC

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
     &                                5384.21 / TSFC( S ) ) )
                QSFC  = 0.622 * ES / ( PSFC - ES )
                TVSFC = TSFC( S ) * ( 1.0 + 0.6077 * QSFC )
                THETG = TVSFC * ( 1000.0 / PSFC )**0.286
                THV1  = TV( 1,S )*( 1000.0 / PRES( 1,S ) )**0.286
                DTHDZ( 1 ) = ( THV1 - THETG ) / ZH( 1,S )

                DO 111  L = 2, EMLAYS

                    IF ( HMIX ( S ) .LT. ZH( L-1,S ) ) LPBL = L
                    IF ( STKHT( S ) .LT. ZH( L-1,S ) ) LSTK = L

                    THVK = TV( L,S ) * ( 1000.0 / PRES( L,S ) )**0.286
                    DTHDZ( L ) = DDZH( L-1,S ) * ( THVK - THV1 )
                    THV1 = THVK

111             CONTINUE
                
                P    = MIN( ZSTK( LSTK,S ) * DDZH( LSTK,S ) , 1.0 )
                Q    = 1.0 - P
                TSTK = P * TA  ( LSTK-1,S )  +  Q * TA  ( LSTK,S )
                WSTK = P * WSPD( LSTK-1 )    +  Q * WSPD( LSTK )


C...............   Compute plume rise for this source (noting inversion of MOL):

                CALL PLMRIS( ZH( 1,S ), ZF( 1,S ), WSPD, 
     &           TA( 1,S ), QV( 1,S ), PRES( 1,S ), DTHDZ, USTAR( S ), 
     &           ZRUF( S ), HMIX( S ), 1.0 / MOL( S ), HFX( S ),
     &           STKHT( S ), STKDM( S ), STKTK( S ), STKVE( S ),
     &           ZSTK( 1,S ), TSTK, WSTK, LSTK, LPBL, 
     &           SZ0FAC, SPRFAC, TOP, BOT )
      
     
C...............   Compute plume fractions, assuming uniform distribution
C...............   in Z from bottom to top:

                DDP = 1.0 / ( TOP - BOT )
                ZZ0 = 0

                DO 122  L = 1, EMLAYS

                    IF( ZZ0 .GT. TOP ) THEN
                        LFRAC( S,L ) = 0.0
                    ELSE
                        ZZ1 = ZF( L,S )
                        IF ( ZF( L,S ) .LT. BOT ) THEN
                            LFRAC( S,L ) = 0.0
                        ELSE
                            P = MIN( ZZ1, TOP ) - MAX( ZZ0, BOT )
                            LFRAC( S,L ) = P * DDP
                        END IF
                        ZZ0 = ZZ1
                    END IF

122             CONTINUE

133         CONTINUE

            IF ( .NOT. WRITE3( LNAME, 'LFRAC', 
     &                         JDATE, JTIME, LFRAC ) ) THEN
                MESG = 'Error writing LFRAC to ' // LNAME
                CALL M3EXIT( 'LAYPOINT', JDATE, JTIME, MESG, 2 )
            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )

144     CONTINUE

999     CONTINUE          !  exit program

        MESG = 'Normal completion of LAYPOINT'
        CALL M3EXIT( 'LAYPOINT', 0, 0, MESG, 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, I10 )

92011   FORMAT (/5X , A, I10 )

92020   FORMAT ( 5X , A, I9 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93020   FORMAT( I5.5, 1X, I4, 1X, I3, 1X, I10, 1X, I5, 3 I4 )


C...........   Internal buffering formats.............94xxx

94010   FORMAT( 12( A, :, I8, :, 1X ) )

        END

