
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
C       SOLAR() and ZANGLE() for SMOKE-BEIS2:  produces COS( ZENITH )
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


C.......   PARAMETERS:
 
      REAL*4      AA, BB, CC, SIGA, D60, D15, D24, ROTDAY , SDEC
      PARAMETER ( AA =   0.15    ,
     &            BB =   3.885   ,
     &            CC = - 1.253   ,
     &          SIGA = 279.9348  ,
     &          SDEC =   0.39784984 , !  SIN (23^26'37.8") the declination angle
     &           D60 = 1.0 / 60.0,
     &           D15 = 1.0 /15.0 ,
     &           D24 = 1.0 /24.0 ,
     &        ROTDAY = 360.0 / 365.242 ! fraction of a complete rotation per day
     &           )
 

C.........   Scratch Local variables

      INTEGER*4    R, C
      REAL*4       SLA, GMT,  TK, DAD, DF,
     &             DESIN, DECOS, DESIN2, DECOS2, SIG, DECSIN, DECCOS,
     &             EQT, TST, HRANGL
                     
C.......   SAVED local variables

        LOGICAL  FIRSTIME
        DATA     FIRSTIME / .TRUE. /
        REAL     SINLAT( NCOLS,NROWS )
        REAL     COSLAT( NCOLS,NROWS )
        
        SAVE     FIRSTIME, SINLAT, COSLAT


C................  Begin body of program  .........................

      IF ( FIRSTIME ) THEN

          FIRSTIME = .FALSE.

          DO  12  R = 1, NROWS
          DO  11  C = 1, NCOLS
              SLA = PI180 * LAT( C,R )
              SINLAT( C,R ) = SIN( SLA )
              COSLAT( C,R ) = COS( SLA )
11        CONTINUE
12        CONTINUE

      END IF	!  if firstime


C.......   Convert time to hours and add time-zone offset
      
      GMT    =  FLOAT( JTIME / 10000 ) +			!  hr part
     &            D60 * ( FLOAT( MOD( JTIME / 100 , 100 ) ) 	!  min part
     &                  + D60 * FLOAT( MOD( JTIME, 100 ) ) )	!  sec part
      DAD    =  GMT * D24 + MOD( JDATE , 1000 )
      DF     =  ROTDAY * PI180 * DAD      !  The terrestrial-rotation angle
           
      DESIN  =  SIN( DF )          !  SINE   of this angle
      DECOS  =  COS( DF )          !  COSINE of this angle
           
      DESIN2 =  SIN( DF + DF )     !  SINE   of twice the angle
      DECOS2 =  COS( DF + DF )     !  COSINE of twice the angle
           
      SIG  =  DF +
     &        PI180 * ( SIGA +
     &                  1.914827 * DESIN  - 0.079525 * DECOS +
     &                  0.019938 * DESIN2 - 0.00162  * DECOS2 )
           
C  The sine of the declination
      
      DECSIN = SDEC * SIN( SIG )
      DECCOS = SQRT( 1.0 - DECSIN*DECSIN )
           
C  The equation of time adjustment
      
      EQT  =  0.123470 * DESIN  - 0.004289 * DECOS
     &      + 0.153809 * DESIN2 + 0.060783 * DECOS2

      ZFLAG = .FALSE.
      
      DO 101 R = 1, NROWS
      DO 100 C = 1, NCOLS

          TK     =  GMT + LON( C,R ) * D15      !  Distance in hours from LON=0
          TST    =  TK - EQT                    !  true solar time
          HRANGL =  PI180 * 15.0 * ABS(TST - 12.0) !  hour angle
                   
C  Compute the sine of the solar elevation
          
          COSZEN( C,R ) = DECSIN * SINLAT( C,R ) +
     &                    DECCOS * COSLAT( C,R ) * COS( HRANGL )
          ZFLAG = ( ZFLAG .OR. ( COSZEN( C,R ) .GT. -0.01 ) )  ! "sun > horizon"
                   
100   CONTINUE
101   CONTINUE


      RETURN
      END
