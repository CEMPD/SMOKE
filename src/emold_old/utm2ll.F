
C Version "@(#)$Id$ $Source$ $Date$ 
      
      SUBROUTINE UTM2LL( X, Y, Z, LON, LAT )

C***********************************************************************
C  subroutine body starts at line  85
C
C  FUNCTION:
C     Convert UTM zone-Z coordinates X-Y to LAT-LON coords
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       USGS National Mapping Division General Cartographic Transformation
C       Package, routine GTPZ0()
C
C  REVISION  HISTORY:
C       Prototype 6/95 by CJC, adapted from UAM EPS MAPUTG() 
C       X,Y are now in meters.
C       Version 10/95 uses GTPZ0()
C***********************************************************************

      IMPLICIT NONE


C...........   ARGUMENTS:

      REAL*4    X       !  UTM easting  in meters
      REAL*4    Y       !  UTM northing in meters
      INTEGER*4 Z       !  UTM zone
      REAL*4    LON     !  East longitude in decimal degrees
      REAL*4    LAT     !  North latitude in decimal degrees

CC...........   PARAMETERS:
C
C      REAL*8     PI 
C      REAL*8      RPI180
C      PARAMETER ( PI     = 3.14159 26535 89793 23846 26433 83279D0 ,
C     &            RPI180 = 180.0D0 / PI )


C...........   External Functions

        INTEGER         INIT3           !  from M3IO
        EXTERNAL        INIT3


C.......   LOCAL VARIABLES:
C.......   Arguments for GTPZ0:

      REAL*8            CRDIN( 2 )      !  input coordinates x,y
      INTEGER*4         INSYS           !  input projection code
      INTEGER*4         INZONE          !  input utm zone, etc.
      REAL*8            TPARIN( 15 )    !  input projection parameters
      INTEGER*4         INUNIT          !  input units code
      INTEGER*4         INSPH           !  spheroid code
      INTEGER*4         IPR             !  error print flag
      INTEGER*4         JPR             !  projection parameter print flag
      INTEGER*4         LEMSG           !  error message unit number
      INTEGER*4         LPARM           !  projection parameter unit number
      REAL*8            CRDIO( 2 )      !  output coordinates x,y
      INTEGER*4         IOSYS           !  output projection code
      INTEGER*4         IOZONE          !  output utm zone, etc.
      REAL*8            TPARIO( 15 )    !  output projection parameters
      INTEGER*4         IOUNIT          !  output units code
      INTEGER*4         LN27            !  NAD1927 file unit number
      INTEGER*4         LN83            !  NAD1983 file unit number
      CHARACTER*128     FN27            !  NAD1927 file name
      CHARACTER*128     FN83            !  NAD1983 file name
      INTEGER*4         LENGTH          !  NAD* record-length
      INTEGER*4         IFLG            !  error flag

C.......   Error codes for GTPZ0:

      CHARACTER*64      MESG( 9 )
      DATA              MESG /
     &  'Illegal input system code INSYS',
     &  'Illegal output system code IOSYS',
     &  'Illegal input unit code INUNIT',
     &  'Illegal output unit code IOUNIT',
     &  'Inconsistent unit and system codes for input',
     &  'Inconsistent unit and system codes for output',
     &  'Illegal input zone code INZONE',
     &  'Illegal output zone code IOZONE',
     &  'Projection-specific error' /



C***********************************************************************
C   begin body of subroutine  UTM2LL()

C.......   Set up input arguments for GTPZ0()

        CRDIN( 1 ) = DBLE( X )
        CRDIN( 2 ) = DBLE( Y )
        INSYS  = 1
        INZONE = Z
        INUNIT = 2       !  meters
        INSPH  = 8       !  GRS 1980 spheroid
        IPR    = 0       !  print error messages, if any
        JPR    = 1       !  do NOT print projection parameters
        LEMSG  = INIT3() !  unit number for log file
        LPARM  = LEMSG   !  projection parameters file
        IOSYS  = 0       !  geographic (lat-lon)
        IOUNIT = 4       !  output units: degrees
        
C.......   Call GTPZ0()

        CALL GTPZ0( CRDIN, INSYS, INZONE, TPARIN, INUNIT, INSPH, 
     &              IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &              TPARIO, IOUNIT, LN27, LN83, FN27, FN83, LENGTH, 
     &              IFLG )

        IF ( IFLG .NE. 0 ) THEN
            IFLG = MAX( MIN( 9, IFLG ), 1 )     !  trap between 1 and 9
            CALL M3WARN( 'UTM2LL', 0,0, MESG( IFLG ) )
        END IF

C.......   Decode output arguments for GTPZ0()

      LON = SNGL( CRDIO( 1 ) )
      LAT = SNGL( CRDIO( 2 ) )


        RETURN
        END

