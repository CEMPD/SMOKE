
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  CZANGLE( JDATE, JTIME, LAT, LON, COSZEN, ZFLAG )

C***********************************************************************
C  subroutine body starts at line  84
C
C  FUNCTION:
C       Computes cosine of zenith angle for routine HRBIO()
C
C  PRECONDITIONS REQUIRED:
C       JDATE:JTIME represented in GMT
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC, adapted from UAM-BEIS subroutines
C	SOLAR() and ZANGLE() for SMOKE-BEIS2:  produces COS( ZENITH )
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'CONST3.EXT'      ! physical and mathematical constants
      INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants
      INCLUDE 'GRDIMS3.EXT'     ! grid-constants


C...........   ARGUMENTS and their descriptions:

        INTEGER     JDATE   !  current simulation date (YYYYDDD)
        INTEGER     JTIME   !  current simulation time (HHMMSS) (GMT)
        REAL        LAT( NCOLS, NROWS )     !  grid lat (deg) -90 <= LAT <= 90
        REAL        LON( NCOLS, NROWS )     !  grid lon (deg) -180 <= LON <= 180 
        REAL        COSZEN( NCOLS, NROWS )  !  cos of zenith angle
        LOGICAL     ZFLAG                   !  iff sun above horizon somewhere


C...........   PARAMETERS and their descriptions:

        REAL        SIGA, D60, D15, D24, 
     &              ROTDAY , ROTHR, ROTLON, SDEC

        PARAMETER ( 
     &            SIGA = 279.9348  ,
     &            SDEC =   0.39784984 , ! SIN (23^26'37.8") the declination angle
     &             D60 = 1.0 / 60.0,
     &             D15 = 1.0 /15.0 ,
     &             D24 = 1.0 /24.0 ,
     &          ROTDAY = 360.0 / 365.242,   ! rotation per day
     &          ROTHR  = D24 * ROTDAY,      ! rotation per hour
     &          ROTLON = PI180 * ROTHR )    ! rot (radians) per degree lon


C...........   SAVED LOCAL VARIABLES and their descriptions:
C...........   NOTE:  the ANSI standard requires the use of SAVE statements
C...........   for variables which must retain their values from call to call.

        LOGICAL     FIRSTIME
        DATA        FIRSTIME / .TRUE. /

        REAL        SINLAT( NCOLS, NROWS )
        REAL        COSLAT( NCOLS, NROWS )

        SAVE   FIRSTIME, SINLAT, COSLAT


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        REAL        TIMLOC, T, DADM, DF, DESIN, DECOS, DESIN2, DECOS2, 
     &              SIG, DECSIN, DECCOS, EQT, TST, HRANGL, CZ

        INTEGER    R, C
 
 
C***********************************************************************
C   begin body of subroutine  CZANGLE

      IF ( FIRSTIME ) THEN

          FIRSTIME = .FALSE.

          DO  12  R = 1, NROWS
          DO  11  C = 1, NCOLS

              T = PI180 * LAT( C, R )
              SINLAT( C, R ) = SIN( T )
              COSLAT( C, R ) = COS( T )

11        CONTINUE
12        CONTINUE

      END IF          !  if firstime

      TIMLOC = ROTHR * ( FLOAT( JTIME / 10000 ) +
     &                   D60 * ( FLOAT( MOD( JTIME / 100 , 100 ) ) +
     &                           D60 * FLOAT( MOD( JTIME , 100 ) ) ) )
      DADM   = ROTDAY * FLOAT( MOD( JDATE , 1000 ) ) + TIMLOC

      ZFLAG  = .FALSE.  !  is sun anywhere above horizon?

      DO  102  R = 1, NROWS
      DO  101  C = 1, NCOLS

            T      = D15 * LON( C, R )
            DF     = ROTLON * T +  DADM
            DESIN  =  SIN( DF )             !  SINE   of this angle
            DECOS  =  COS( DF )             !  COSINE of this angle
            DESIN2 =  2.0 * DESIN * DECOS            !  SINE   of double-angle
            DECOS2 =  DECOS * DECOS - DESIN * DESIN  !  COSINE of double-angle
            SIG    =  SIGA + DF
     &                     + 1.914827 * DESIN  - 0.079525 * DECOS
     &                     + 0.019938 * DESIN2 - 0.00162  * DECOS2

C...........   The equation of time adjustment, uses sine of the declination

            EQT  =  0.123470 * DESIN  - 0.004289 * DECOS
     &            + 0.153809 * DESIN2 + 0.060783 * DECOS2

C...........   The true solar time

            TST    = T + TIMLOC - EQT
            HRANGL = 15.0 * PI180 * ( TST - 12.0 )

C...........   Compute the sine of the solar elevation (= cos of zenith angle)

            DECSIN = SDEC * SIN( PI180 * SIG )
            DECCOS = SQRT( 1.0 - DECSIN )

            CZ = DECSIN * SINLAT( C, R ) + 
     &           DECCOS * COSLAT( C, R ) * COS( HRANGL )

            COSZEN( C, R ) = CZ

            ZFLAG = ( ZFLAG .OR. ( CZ .GT. -0.01 ) )  !  "sun above horizon"

101   CONTINUE
102   CONTINUE

      RETURN
      END

