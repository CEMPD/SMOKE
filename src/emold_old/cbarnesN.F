
C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CBARNESN( NG, LAT, LON,
     &                     NP, NV, YLAT, XLON, Z, WL50, GRID )

C***********************************************************************
C    subroutine body starts at line  168
C
C  FUNCTION:
C
C    This routine grids multiple variables (or layers) Z( NP,1...NV ) 
C    at a time from locations with lat-lon coordinates 
C    ( YLAT(1...NP), XLON(1...NP) ) and produces output array 
C    GRID( NG,NV ) on the locations ( LAT(1...NG), LON(1...NG) ).
C    For gridding purposes, NG should be NCOLS*NROWS; CBARNESN will
C    regard the corresponding array as singly-indexed.
C
C    This routine grids multiple variables (or layers) from 
C    source data located at latitudes and longitudes YLAT,XLON.
C    to the locations specified by LAT,LON)
C
C    This routine is an extension of a spatial analysis technique
C    with scale dependent filtering that was originally proposed
C    by Stanly Barnes in 1964 and was expanded in 1973.  The scale
C    dependent response function of the filtering effect is
C    analytically calculable and can be adjusted by the choice of
C    two parameters in the Gaussian weighting function used.
C
C    Input data point locations should be specified in
C    latitude and longitude degrees.  The separation distances
C    between these data points and the grid points where the
C    estimates are made is calculated from a function that is
C    determined by spherical geometry under the assumption of a
C    spherical Earth.
C
C  REVISION HISTORY:
C    CGRID1():
C    5/88   Modified for ROMNET
C    9/88   Modified by CJC -- Subexpression elimination and other code
C           improvements; precalculation of scanning radii
C    7/90   Modified by CJC  for ROM 2.2 -- use weights which are Gaussian
C           where exp ( AR**2 ) > 10**-30; 1/R otherwise.  Error exit via EXWST.
C           Get PI and earth-radius related constants from PI.EXT
C    2/91   Adapted by CJC  from BGRID1.FOR -- uses Barnes analysis to
C           data locations in a single corrective iteration, rather than
C           multiple biquadratic back-interpolation from the predicted grid.
C   ?????   Modified by ?Steve Fudge? for UAM BEIS
C   CBARNESN():
C   12/95   Modified by CJC to fit EDSS/Models-3 conventions from 
C           the UAM BEIS CGRID1().  Uses "exact" spherical-geometry
C           distance formula.
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'CONST3.EXT'

C.......   ARGUMENTS:

        INTEGER    NG               !  Number of output values (=NCOLS*NROWS)
        REAL       LAT( NG )        !  single-indexed output latitudes
        REAL       LON( NG )        !  single-indexed output longitudes
        INTEGER    NP               !  number of input sources
        INTEGER    NV               !  number of input variables
        REAL       YLAT( NP )       !  input latitudes
        REAL       XLON( NP )       !  input longitudes
        REAL       Z   ( NP, NV )   !  input values
        REAL       WL50             !  50 Percent filtered wavelength  (KM)
        REAL       GRID( NG, NV )   !  output values


C...........   EXTERNAL FUNCTIONS:

#ifndef _CRAY
        INTEGER         MALLOC
#endif


C...........   PARAMETERS:

        REAL       G, GINV, R0LN50, ALOG10, R2DSQ, AC

      PARAMETER   ( G      =  0.4          ,
     &              GINV   =  2.5          ,
     &              R0LN50 =  2.1058923    ,
     &              ALOG10 =  2.3025850929 ,
     &              R2DSQ  =  RPI180 * RPI180 ,
     &              AC     =  30.0 * ALOG10 / D2KMSQ )


C...........   LOCAL VARIABLES:

        INTEGER    I, J, M, N    !  cell, variable, station counters
        REAL       C             !  scaling-filter constant
        REAL       C4K           !  Gauss constant for first pass
        REAL       GC4K          !  Gauss constant for second pass
        REAL       YG , XG       !  temporaries for lat, lon of current cell
        REAL       W1, W2        !  weights for Gauss-weighted average
        REAL       WTOT1, WTOT2  !  sum of weights
        REAL       FTOT1( NV )   !  accumulator for values
        REAL       FTOT2( NV )   !  accumulator for corrections
        REAL       DELX          !  lon difference
        REAL       DSQ           !  DELX**2  +  DELY**2
        REAL       XA, YA        !  LON, LAT of current station
        REAL       XB, YB        !  LON, LAT of current station
        REAL       CY, SY        !  cos, sin of lat
        REAL       RMAXA         !  maximum scanning radii, for first
        REAL       RMAXB         !  and second passes
        REAL       ANUMA         !  numerator, beyond scanning radius,
        REAL       ANUMB         !  for first and second passes
        REAL       DZ( NP, NV )  !  work array
        CHARACTER*250  MESG

C...........   STATE VARIABLES:

        REAL         WL501

#ifndef _CRAY
        INTEGER         SIZE, LSIZE
        POINTER         ( P, DZ )
        POINTER         ( P1, FTOT1 )
        POINTER         ( P2, FTOT2 )
 
        DATA            LSIZE  /  -1          /  !  last value of N
        DATA            WL501  /  - 1.0       /  !  last value of WL50
        SAVE            P, P1, P2, LSIZE, WL501
#endif

#ifdef _CRAY
        DATA            WL501  /  - 1.0       /  !  last value of WL50
        SAVE            WL501
#endif


C..................................................................
C.......   begin body of  CBARNESN:

#ifndef _CRAY
        SIZE = NP * NV
        IF ( SIZE .GT. LSIZE ) THEN

            IF ( LSIZE .GT. 0 )  THEN
                CALL FREE( P )
                CALL FREE( P1 )
                CALL FREE( P2 )
            END IF

            P = MALLOC( 8 * SIZE )
            IF ( P .EQ. 0 )
     &           CALL M3EXIT( 'TSHIFT', 0, 0,
     &                        'Memory allocation error for P', 2 )
            P1 = MALLOC( 8 * NV )
            IF ( P1 .EQ. 0 )
     &           CALL M3EXIT( 'TSHIFT', 0, 0,
     &                        'Memory allocation error for P1', 2 )
            P2 = MALLOC( 8 * NV )
            IF ( P2 .EQ. 0 )
     &           CALL M3EXIT( 'TSHIFT', 0, 0,
     &                        'Memory allocation error for P2', 2 )
            LSIZE = SIZE

        END IF		!  if size increased
#endif

C...........   Calculate the scaling filter parameter

      IF  ( WL50 .NE. WL501 )  THEN

          C  =  R0LN50 * ( WL50 / PI )**2

          C4K   = -D2KMSQ / C
          RMAXA =  C * AC           !   exp ( c4k * rmaxa ) = 1.0e-30
          ANUMA =  1.0E-30 * RMAXA

          GC4K  =  GINV * C4K
          RMAXB =  RMAXA * G
          ANUMB =  1.0E-30 * RMAXB

          WL501 = WL50                        !  Record "WL50"

      END IF


C.......   Prediction pass:  Scan each input data point and construct
C.......   estimated error  DZ  at that point:

      DO 155   N = 1 , NP

          DO  111  I = 1 , NV
              FTOT1( I ) = 0.0
111       CONTINUE              !  end loop on variables I

          XA    =  PI180 * XLON( N )
          YA    =  PI180 * YLAT( N )
          CY    =  COS( YA )                
          SY    =  SIN( YA )
          WTOT1 =  0.0

          DO  133  M = 1 , NP

              XB = PI180 * XLON( M )
              YB = PI180 * YLAT( M )
              DELX  =  XB  -  XA
              DSQ   =  R2DSQ * 
     &                 ACOS( COS( DELX ) * CY * COS( YB ) +
     &                                     SY * SIN( YB ) )**2

              IF ( DSQ .LE. RMAXA )  THEN
                  W1     =  EXP ( DSQ * C4K )
              ELSE          !  use matching  1/R**2  weight
                  W1     =  ANUMA / DSQ
              END IF
              WTOT1  =  WTOT1  +  W1
              
              DO  122  I = 1, NV
                  FTOT1( I )  =  FTOT1( I )  +  W1 * Z ( M,I )
122           CONTINUE              !  end loop on variables I

133       CONTINUE              !  end loop on sites M

          WTOT1 = 1.0 / WTOT1

          IF ( WTOT1 .EQ. 0.0  ) THEN
              WRITE ( MESG,94010 ) N, YA, XA, NP
              CALL M3EXIT( 'CBARNESN', 0, 0, MESG, 2 )
          END IF

          DO  144  I = 1 , NV
              DZ ( N,I ) = Z ( N,I )  -  FTOT1( I ) * WTOT1
144       CONTINUE              !  end loop on variables I

155   CONTINUE          !  end prediction loop on sites N


C.......   Grid-prediction pass:  generate estimate using first set of
C.......   weights, and correction using error estimates  DZ and second
C.......   set of weights

      DO  255  J = 1 , NG

          DO  211  I = 1 , NV
              FTOT1( I ) = 0.0
              FTOT2( I ) = 0.0
211       CONTINUE              !  end loop on variables I

          YG   =  PI180 * LAT( J )    !  LAT of grid nodes at index J
          XG   =  PI180 * LON( J )    !  LON of grid nodes at index J
          CY   =  COS( YG )
          SY   =  SIN( YG )
           
C.......   Scan each input data point
      
          WTOT1 = 0.0
          WTOT2 = 0.0

          DO  233  N = 1, NP
           
              XA    =  PI180 * XLON( N )
              YA    =  PI180 * YLAT( N )
              DELX  =  XA  -  XG
              DSQ   =  R2DSQ * 
     &                 ACOS( COS( DELX ) * CY * COS( YA ) +
     &                                     SY * SIN( YA ) )**2

              IF ( DSQ .LE. RMAXB )  THEN
                  W1  =  EXP ( DSQ * C4K )
                  W2  =  EXP ( DSQ * GC4K )
              ELSE IF ( DSQ .LE. RMAXA )  THEN
                  W1  =  EXP ( DSQ * C4K )
                  W2  =  ANUMB / DSQ
              ELSE          !  use matching  1/R**2  weight
                  W1  =  ANUMA / DSQ
                  W2  =  G * W1     ! = ANUMB / DSQ
              END IF

              WTOT1  =  WTOT1  +  W1
              WTOT2  =  WTOT2  +  W2

              DO  222  I = 1 , NV
                  FTOT1( I ) = FTOT1( I )  +  W1 * Z ( N,I )
                  FTOT2( I ) = FTOT2( I )  +  W2 * DZ( N,I )
222           CONTINUE              !  end loop on variables I

233       CONTINUE          !  end loop on data sites nn

          IF ( WTOT1 .EQ. 0.0  .OR.  WTOT2 .EQ. 0.0 ) THEN
              WRITE ( MESG,94020 ) J, I, NP
              CALL M3EXIT( 'CBARNESN', 0, 0, MESG, 2 )
          END IF

          WTOT1 = 1.0 / WTOT1
          WTOT2 = 1.0 / WTOT2
          DO  244  I = 1 , NV
              GRID( J,I ) = WTOT1 * FTOT1( I )  +  WTOT2 * FTOT2( I )
244       CONTINUE              !  end loop on variables I

255   CONTINUE          !  end loop on rows J


      RETURN

C*************************  MESSAGE  FORMATS  **************************


C.......   Internal buffering formats     94xxx


94010 FORMAT ( 'Weights total zero at site', I5, 2X, 
     &         'at lat-lon', F7.2, ':', F7.2, 2X,
     &         'input count=', I6 )

94020 FORMAT ( 'Weights total zero at (C,R)=(', I5, ',', I5, 2X,
     &         'input count=', I6 )


      END
