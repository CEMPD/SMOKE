C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION ETPRED(IEVTYP,KEYEQ,ETPRVP,ETPTMP)
C
C  ETPRED computes the Hot Soak/Running Loss reductions for fail vehicles.
C
C  Called by ETPCAL.
C
C  Input on call:
C    parameter list: IEVTYP,KEYEQ,ETPRVP,ETPTMP
C    /ETPPAR/ ETP,NEWETP
C
C  Output on return:
C
C    function: ETPRED
C
C  Local array subscripts:
C
C  RED(5,4)  -  RED ( IRVP, ITEMP )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  DT       R    Difference between input temperature and reference temp of 65F
C  ETPRVP   R    Weathered RVP value
C  ETPTMP   R    Hot Soak or Running Loss temperature value
C  RED      R    Hot Soak/Running Loss emission reductions (%) at
C                6, 7, 8, 9, and 10 RVP and 65, 80, 95 and 110F
C  IRVP1    I    Upper left RVP point index
C  IRVP2    I    Upper right RVP point index
C  ITEMP1   I    Lower left temperature point index
C  ITEMP2   I    Lower right temperature point index
C  I9RVP    I    Index into the RED array at 9.0 RVP
C  I95TMP   I    Index into the RED array at 95F
C  MAXRVP   I    Maximum RVP dimension of RED array
C  MAXTMP   I    Maximum Temperature dimension of RED array
C  RBASE    R    Lowest allowed RVP
C  TBASE    R    Lowest allowed temperature
C  TINC     R    Increment between successive known temperature values
C  RINC     R    Increment between successive known RVP values
C  DR       R    Difference between known and input RVP values
C  DT       R    Difference between known and input temperature values
C  A        R    Upper interpolation point with respect to RVP only
C  B        R    Lower interpolation point with respect to RVP only
C
C  Notes:
C
C  ETPRED was added in MOBILE5v3.
C  Due to its volatility, the RED matrix is initialized here
C  rather than in BLOCK DATA.
C
C
      COMMON /ETPPAR/ ETP(7,3),NEWETP,IETPYR(4),ETPPIV(4),NETPYR
C
      DIMENSION RED(5,4)
      DATA RED/ 0., 0.,  0.,  0.,  0.,
     *         15., 15., 15., 15., 15.,
     *         30., 30., 30., 30., 30.,
     *         30., 30., 30., 30., 30./
      DATA RBASE,TBASE,RINC,TINC/6.,65.,1.,15./
      DATA MAXRVP,MAXTMP/5,4/
C
C  Find initial location in RED matrix.
C
      IRVP1=(ETPRVP-RBASE)/RINC+1
      ITEMP1=(ETPTMP-TBASE)/TINC+1
C
      IF(IRVP1.LT.1) IRVP1=1
      IF(ITEMP1.LT.1) ITEMP1=1
      IF(IRVP1.GT.MAXRVP) IRVP1=MAXRVP
      IF(ITEMP1.GT.MAXTMP) ITEMP1=MAXTMP
C
C  Find interpolation location in RED matrix.
C
      IRVP2=IRVP1+1
      ITEMP2=ITEMP1+1
C
      IF(IRVP2.GT.MAXRVP) IRVP2=MAXRVP
      IF(ITEMP2.GT.MAXTMP) ITEMP2=MAXTMP
C
C  Perform 4-point linear interpolation.
C
      DR=ETPRVP-(RBASE+(IRVP1-1)*RINC)
      DT=ETPTMP-(TBASE+(ITEMP1-1)*TINC)
C
      A=RED(IRVP1,ITEMP1)+DR*(RED(IRVP2,ITEMP1)-RED(IRVP1,ITEMP1))/RINC
      B=RED(IRVP1,ITEMP2)+DR*(RED(IRVP2,ITEMP2)-RED(IRVP1,ITEMP2))/RINC
C
      ETPRED=A+DT*(B-A)/TINC
C
C  Adjust reduction based on user-supplied input at 9 RVP and 95F.
C
      IF(NEWETP.EQ.2) THEN
        I9RVP=(9.0-RBASE)/RINC+1
        I95TMP=(95.0-TBASE)/TINC+1
        ETPRED=ETPRED*ETP(IEVTYP,KEYEQ)/RED(I9RVP,I95TMP)
      ENDIF
C
      RETURN
      END
