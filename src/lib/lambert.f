C.. NOTE: Had to create this version because need for Models-3, but GRIDDESC
C         file is not being used :-(

C.........................................................................
C Version "@(#)$Id$ $Source$ $Date$ "
C EDSS/Models-3 I/O API.  Portions copyright (C) 1992-1997 MCNC
C See file "COPYRIGHT.txt" for conditions of use.
C.........................................................................

C Version "@(#)lambert.F	1.1 /pub/storage/edss/framework/src/smoke/emlib/SCCS/s.lambert.F 18 Sep 1996 12:38:32
      
        LOGICAL FUNCTION LAMBERT( P_ALP_IN, P_BET_IN, P_GAM_IN, 
     &                            XCENT_IN, YCENT_IN )

        IMPLICIT NONE

        LOGICAL          SETLAM
        LOGICAL          LL2LAM
        LOGICAL          LAM2LL
        LOGICAL          UTM2LAM
        LOGICAL          LAM2UTM

C***********************************************************************
C  subroutine LAMBERT body starts at line  157
C  entry      SETLAM       starts at line  236
C  entry      LAM2LL       starts at line  303
C  entry      LL2LAM       starts at line  356
C  entry      LAM2UTM      starts at line  414
C  entry      UTM2LAM      starts at line  469
C
C  FUNCTION:
C     LAMBERT:  set up  GTPZ0() for a particular named Lambert.
C               If CNAME is a _grid_ name, returns corresponding
C               _coordinate_system_name_ and coordinate definition 
C               parms A,B,C,X,Y.
C     SETLAM:   set up  GTPZ0() for a argument-specified anonymous Lambert.
C     LL2LAM:   Convert LAT-LON coordinates to Lambert coordinates
C     LAM2LL:   Convert Lambert coordinates to LAT-LON coordinates
C
C  PRECONDITIONS REQUIRED:
C       For LAMBERT(), CNAME must be the name either of a coordinate
C       system or a grid found in file GRIDDESC; furthermore, the 
C       projection-type must be LAMGRD3 (i.e., Lambert)
C       Must call LAMBERT() or SETLAM() before calling conversion
C       functions LL2LAM(), LAM2LL(), UTM2LAM(), or LAM2UTM().
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       USGS National Mapping Division General Cartographic Transformation
C       Package, routine GTPZ0().
C       I/O API routines INIT3(), DSCGRID(), DSCOORD(), M3MESG()
C
C  REVISION  HISTORY:
C       Prototype 11/95 by CJC:  uses GTPZ0()
C***********************************************************************

        INCLUDE 'PARMS3.EXT'

C...........   ARGUMENTS:

        REAL*8        P_ALP_IN     !  first, second, third map
        REAL*8        P_BET_IN     !  projection descriptive
        REAL*8        P_GAM_IN     !  parameters
        REAL*8        XCENT_IN     !  lon for coord-system X=0
        REAL*8        YCENT_IN     !  lat for coord-system Y=0
        REAL*4        A          !  first secant latitude
        REAL*4        B          !  second secant latitude.  B > A
        REAL*4        C          !  central meridian
        REAL*4        X          !  Lambert easting  in meters
        REAL*4        Y          !  Lambert northing in meters
        REAL*4        U          !  UTM easting  in meters
        REAL*4        V          !  UTM northing in meters
        REAL*4        LON        !  East longitude in decimal degrees
        REAL*4        LAT        !  North latitude in decimal degrees
        INTEGER       Z          !  UTM zone (1...36)


C...........   PARAMETERS:

      REAL*8      D60
      REAL*8      PI 
      REAL*8      RPI180
      PARAMETER ( D60    = 1.0D0 / 60.0D0,
     &            PI     = 3.14159 26535 89793 23846 26433 83279D0 ,
     &            RPI180 = 180.0D0 / PI )


C...........   External Functions

        LOGICAL         DSCOORD
        LOGICAL         DSCGRID
        INTEGER         INIT3           !  from M3IO
        EXTERNAL        DSCOORD, DSCGRID, INIT3


C.......   LOCAL VARIABLES:
C.......   Arguments for GTPZ0:

      REAL*8            CRDIN( 2 )      !  input coordinates x,y
      INTEGER*4         INSYS           !  input projection code
      INTEGER*4         INZONE          !  input utm zone, etc.
      REAL*8            TPAIN( 15 )     !  input projection parameters
      INTEGER*4         INUNIT          !  input units code
      INTEGER*4         INSPH           !  spheroid code
      INTEGER*4         IPR             !  error print flag
      INTEGER*4         JPR             !  projection parameter print flag
      INTEGER*4         LEMSG           !  error message unit number
      INTEGER*4         LPARM           !  projection parameter unit number
      REAL*8            CRDIO( 2 )      !  output coordinates x,y
      INTEGER*4         IOSYS           !  output projection code
      INTEGER*4         IOZONE          !  output utm zone, etc.
      REAL*8            TPARO( 15 )     !  output projection parameters
      INTEGER*4         IOUNIT          !  output units code
      INTEGER*4         LN27            !  NAD1927 file unit number
      INTEGER*4         LN83            !  NAD1983 file unit number
      CHARACTER*128     FN27            !  NAD1927 file name
      CHARACTER*128     FN83            !  NAD1983 file name
      INTEGER*4         LENGTH          !  NAD* record-length
      INTEGER*4         IFLG            !  error flag

C.......   Error codes for GTPZ0:

      CHARACTER*64      GMESG( 9 )
      DATA              GMESG /
     &  'Illegal input system code INSYS',
     &  'Illegal output system code IOSYS',
     &  'Illegal input unit code INUNIT',
     &  'Illegal output unit code IOUNIT',
     &  'Inconsistent unit and system codes for input',
     &  'Inconsistent unit and system codes for output',
     &  'Illegal input zone code INZONE',
     &  'Illegal output zone code IOZONE',
     &  'Projection-specific error' /

C.......   Arguments for DSCGRID() and DSCCORD()

        CHARACTER*16  ANAME
        INTEGER       CTYPE
        REAL*8        P_ALP     !  first, second, third map
        REAL*8        P_BET     !  projection descriptive
        REAL*8        P_GAM     !  parameters
        REAL*8        XCENT     !  lon for coord-system X=0
        REAL*8        YCENT     !  lat for coord-system Y=0
        REAL*8        XORIG     !  X-coordinate origin of grid (map units)
        REAL*8        YORIG     !  Y-coordinate origin of grid
        REAL*8        XCELL     !  X-coordinate cell dimension
        REAL*8        YCELL     !  Y-coordinate cell dimension
        INTEGER       NCOLS     !  number of grid columns
        INTEGER       NROWS     !  number of grid rows
        INTEGER       NTHIK     !  BOUNDARY:  perimeter thickness (cells)
        
C.......   Scratch variables:

        CHARACTER*80    MESG
        INTEGER         DEG, MNT

C.......   SAVED Local Variables:

        INTEGER*4         ZONE
        DATA              ZONE / 62 /
        
        SAVE    ZONE, CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT


C***********************************************************************
C.......   LAMBERT():  Set up a particular named Lambert projection:

C.......   Return the projection parameters as REAL*4 A,B,C,X,Y:

        P_ALP = P_ALP_IN
        P_BET = P_BET_IN
        P_GAM = P_GAM_IN
        XCENT = XCENT_IN 
        YCENT = YCENT_IN

c       A = SNGL( P_ALP )
c       B = SNGL( P_BET )
c       C = SNGL( P_GAM )
c       X = SNGL( XCENT )
c       Y = SNGL( YCENT )
        LAMBERT = .TRUE.
        
C.......   Convert from real degrees to GTPZ0() format  dddmmmsss.sssD0

111     CONTINUE

        ZONE  = ZONE + 1
        XCENT = XCENT - P_GAM   !  convert from lon to offset from P_GAM

        DEG   = INT( P_ALP )                            !  int degrees
        P_ALP = 60.0D0 * ( P_ALP - DBLE( DEG ) )        !  minutes
        MNT   = INT( P_ALP )                            !  int minutes
        P_ALP = 60.0D0 * ( P_ALP - DBLE( MNT ) )        !  seconds
        P_ALP = P_ALP + 1000.0D0 * ( MNT + 1000 * DEG ) !  dddmmmsss.sssD0

        DEG   = INT( P_BET )                            !  int degrees
        P_BET = 60.0D0 * ( P_BET - DBLE( DEG ) )        !  minutes
        MNT   = INT( P_BET )                            !  int minutes
        P_BET = 60.0D0 * ( P_BET - DBLE( MNT ) )        !  seconds
        P_BET = P_BET + 1000.0D0 * ( MNT + 1000 * DEG ) !  dddmmmsss.sssD0

        DEG   = INT( P_GAM )                            !  int degrees
        P_GAM = 60.0D0 * ( P_GAM - DBLE( DEG ) )        !  minutes
        MNT   = INT( P_GAM )                            !  int minutes
        P_GAM = 60.0D0 * ( P_GAM - DBLE( MNT ) )        !  seconds
        P_GAM = P_GAM + 1000.0D0 * ( MNT + 1000 * DEG ) !  dddmmmsss.sssD0

        DEG   = INT( XCENT )                            !  int degrees
        XCENT = 60.0D0 * ( XCENT - DBLE( DEG ) )        !  minutes
        MNT   = INT( XCENT )                            !  int minutes
        XCENT = 60.0D0 * ( XCENT - DBLE( MNT ) )        !  seconds
        XCENT = XCENT + 1000.0D0 * ( MNT + 1000 * DEG ) !  dddmmmsss.sssD0

        DEG   = INT( YCENT )                            !  int degrees
        YCENT = 60.0D0 * ( YCENT - DBLE( DEG ) )        !  minutes
        MNT   = INT( YCENT )                            !  int minutes
        YCENT = 60.0D0 * ( YCENT - DBLE( MNT ) )        !  seconds
        YCENT = YCENT + 1000.0D0 * ( MNT + 1000 * DEG ) !  dddmmmsss.sssD0

        RETURN

C.....................................................................
C.......   Set up anonymous lambert from arguments:

       ENTRY SETLAM( A, B, C, X, Y )

C.......   Check validity of input parameters:

        IF ( A .LT. 0.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad first latitude A =', A
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( A .GT. B ) THEN
            WRITE( MESG, 94020 ) 'Bad latitudes A ', A, 'B =', B
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( B .GE.   90.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad second latitude B =', B
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( C .LT. -180.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad central longitude C =', C
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( C .GT.  180.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad central longitude C =', C
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( X .LT. -180.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad origin longitude X =', X
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( X .GT.  180.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad origin longitude X =', X
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( Y .LT.    0.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad origin latitude Y =', Y
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        ELSE IF ( Y .GE.   90.0 ) THEN
            WRITE( MESG, 94020 ) 'Bad origin latitude Y =', Y
            CALL M3WARN( 'LAMBERT/SETLAM', 0, 0, MESG )
            SETLAM = .FALSE.
            RETURN
        END IF

C.......   Convert to double, then go to GTPZ0() format conversion
        
        P_ALP   = DBLE( A )
        P_BET   = DBLE( B )
        P_GAM   = DBLE( C )
        XCENT   = DBLE( X )
        YCENT   = DBLE( Y )
        SETLAM  = .TRUE.

        GO TO  111      !  convert projection parms to dddmmmsss.sssD0 format


C.....................................................................
C.......   convert from Lambert to lat-lon
C.......   Set up input arguments for GTPZ0()

        ENTRY LAM2LL( X, Y, LON, LAT )

        IF ( ZONE .EQ. 0 ) THEN
            CALL M3WARN( 'LAMBERT/LAM2LL', 0, 0, 
     &                   'Projection not initialized' )
            LAM2LL = .FALSE.
            RETURN
        END IF

        CRDIN( 1 ) = DBLE( X )
        CRDIN( 2 ) = DBLE( Y )
        TPAIN( 1 ) = 0.0D0
        TPAIN( 2 ) = 0.0D0
        TPAIN( 3 ) = P_ALP
        TPAIN( 4 ) = P_BET
        TPAIN( 5 ) = P_GAM
        TPAIN( 6 ) = YCENT
        TPAIN( 7 ) = 0.0D0
        TPAIN( 8 ) = 0.0D0
        INSYS  = 4       !  Lambert conformal conic
        INZONE = ZONE
        INUNIT = 2       !  input units:  meters
        INSPH  = 8       !  GRS 1980 spheroid
        IPR    = 0       !  print error messages, if any
        JPR    = 1       !  do NOT print projection parameters
        LEMSG  = INIT3() !  unit number for log file
        LPARM  = LEMSG   !  projection parameters file
        IOSYS  = 0       !  geographic (lat-lon)
        IOUNIT = 4       !  output units: degrees
        
C.......   Call GTPZ0()

        CALL GTPZ0( CRDIN, INSYS, INZONE, TPAIN, INUNIT, INSPH, 
     &              IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &              TPARO, IOUNIT, LN27, LN83, FN27, FN83, LENGTH, 
     &              IFLG )

        IF ( IFLG .NE. 0 ) THEN
            IFLG = MAX( MIN( 9, IFLG ), 1 )     !  trap between 1 and 9
            CALL M3WARN( 'LAMBERT/LAM2LL', 0,0, GMESG( IFLG ) )
        END IF

C.......   Decode output arguments for GTPZ0()

        LON    = SNGL( CRDIO( 1 ) )
        LAT    = SNGL( CRDIO( 2 ) )
        LAM2LL = .TRUE.
        RETURN


C.....................................................................
C.......   Convert from Lat-Lon to Lambert:

        ENTRY  LL2LAM( LON, LAT, X, Y )

C.......   Check initialization:
        
        IF ( ZONE .EQ. 0 ) THEN
            CALL M3WARN( 'LAMBERT/LL2LAM', 0, 0, 
     &                   'Projection not initialized' )
            LL2LAM = .FALSE.
            RETURN
        END IF

C.......   Set up input arguments for GTPZ0()

        CRDIN( 1 ) = DBLE( LON )
        CRDIN( 2 ) = DBLE( LAT )
        INSYS  = 0       !  projection default (lat-lon)
        INUNIT = 4       !  input units:  degrees
        INSPH  = 8       !  GRS 1980 spheroid
        IPR    = 0       !  print error messages, if any
        JPR    = 1       !  do NOT print projection parameters
        LEMSG  = INIT3() !  unit number for log file
        LPARM  = LEMSG   !  projection parameters file
        IOSYS  = 4       !  Lambert conformal conic
        IOZONE = ZONE    !  LAM zone
        IOUNIT = 2       !  output units: meters
        TPARO( 1 ) = 0.0D0
        TPARO( 2 ) = 0.0D0
        TPARO( 3 ) = P_ALP
        TPARO( 4 ) = P_BET
        TPARO( 5 ) = P_GAM
        TPARO( 6 ) = YCENT
        TPARO( 7 ) = 0.0D0
        TPARO( 8 ) = 0.0D0
        

C.......   Call GTPZ0()

        CALL GTPZ0( CRDIN, INSYS, INZONE, TPAIN, INUNIT, INSPH, 
     &              IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &              TPARO, IOUNIT, LN27, LN83, FN27, FN83, LENGTH, 
     &              IFLG )

        IF ( IFLG .NE. 0 ) THEN
            IFLG = MAX( MIN( 9, IFLG ), 1 )     !  between 1 and 9
            CALL M3WARN( 'LAMBERT/LL2LAM', 0,0, GMESG( IFLG ) )
        END IF

C.......   Decode output arguments for GTPZ0()

        X      = SNGL( CRDIO( 1 ) )
        Y      = SNGL( CRDIO( 2 ) )
        LL2LAM = .TRUE.
        RETURN

C.....................................................................
C.......   convert from Lambert to UTM
C.......   Set up input arguments for GTPZ0()

        ENTRY LAM2UTM( X, Y, Z, U, V )

        IF ( ZONE .EQ. 0 ) THEN
            CALL M3WARN( 'LAMBERT/LAM2UTM', 0, 0, 
     &                   'Projection not initialized' )
            LAM2UTM = .FALSE.
            RETURN
        END IF

        CRDIN( 1 ) = DBLE( X )
        CRDIN( 2 ) = DBLE( Y )
        TPAIN( 1 ) = 0.0D0
        TPAIN( 2 ) = 0.0D0
        TPAIN( 3 ) = P_ALP
        TPAIN( 4 ) = P_BET
        TPAIN( 5 ) = P_GAM
        TPAIN( 6 ) = YCENT
        TPAIN( 7 ) = 0.0D0
        TPAIN( 8 ) = 0.0D0
        INSYS  = 4       !  Lambert conformal conic
        INZONE = ZONE
        INUNIT = 2       !  input units:  meters
        INSPH  = 8       !  GRS 1980 spheroid
        IPR    = 0       !  print error messages, if any
        JPR    = 1       !  do NOT print projection parameters
        LEMSG  = INIT3() !  unit number for log file
        LPARM  = LEMSG   !  projection parameters file
        IOSYS  = 1       !  UTM
        IOZONE = Z       !  UTM zone
        IOUNIT = 2       !  output units: meters
        
C.......   Call GTPZ0()

        CALL GTPZ0( CRDIN, INSYS, INZONE, TPAIN, INUNIT, INSPH, 
     &              IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &              TPARO, IOUNIT, LN27, LN83, FN27, FN83, LENGTH, 
     &              IFLG )

        IF ( IFLG .NE. 0 ) THEN
            IFLG = MAX( MIN( 9, IFLG ), 1 )     !  trap between 1 and 9
            CALL M3WARN( 'LAMBERT/LAM2UTM', 0,0, GMESG( IFLG ) )
        END IF

C.......   Decode output arguments for GTPZ0()

        U       = SNGL( CRDIO( 1 ) )
        V       = SNGL( CRDIO( 2 ) )
        LAM2UTM = .TRUE.

        RETURN


C.....................................................................
C.......   Convert from UTM to Lambert:

        ENTRY  UTM2LAM( U, V, Z, X, Y )

C.......   Check initialization:
        
        IF ( ZONE .EQ. 0 ) THEN
            CALL M3WARN( 'LAMBERT/UTM2LAM', 0, 0, 
     &                   'Projection not initialized' )
            UTM2LAM = .FALSE.
            RETURN
        END IF

C.......   Set up input arguments for GTPZ0()

        CRDIN( 1 ) = DBLE( U )
        CRDIN( 2 ) = DBLE( V )
        INSYS  = 1
        INZONE = Z
        INUNIT = 2       !  meters
        INSPH  = 8       !  GRS 1980 spheroid
        IPR    = 0       !  print error messages, if any
        JPR    = 1       !  do NOT print projection parameters
        LEMSG  = INIT3() !  unit number for log file
        LPARM  = LEMSG   !  projection parameters file
        IOSYS  = 4       !  Lambert conformal conic
        IOZONE = ZONE    !  LAM zone
        IOUNIT = 2       !  output units: meters
        TPARO( 1 ) = 0.0D0
        TPARO( 2 ) = 0.0D0
        TPARO( 3 ) = P_ALP
        TPARO( 4 ) = P_BET
        TPARO( 5 ) = P_GAM
        TPARO( 6 ) = YCENT
        TPARO( 7 ) = 0.0D0
        TPARO( 8 ) = 0.0D0
        

C.......   Call GTPZ0()

        CALL GTPZ0( CRDIN, INSYS, INZONE, TPAIN, INUNIT, INSPH, 
     &              IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &              TPARO, IOUNIT, LN27, LN83, FN27, FN83, LENGTH, 
     &              IFLG )

        IF ( IFLG .NE. 0 ) THEN
            IFLG = MAX( MIN( 9, IFLG ), 1 )     !  between 1 and 9
            CALL M3WARN( 'LAMBERT/UTM2LAM', 0,0, GMESG( IFLG ) )
        END IF

C.......   Decode output arguments for GTPZ0()

        X       = SNGL( CRDIO( 1 ) )
        Y       = SNGL( CRDIO( 2 ) )
        UTM2LAM = .TRUE.
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I5, :, 2X )

94020   FORMAT( A, 1PG12.5, :, 2X )

        END


