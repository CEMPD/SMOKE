      REAL FUNCTION ETPRED(ETPRVP,ETPTMP) 
C 
C  ETPRED computes the Hot Soak/Running Loss reductions for fail vehicles. 
C 
C  Called by ETPCAL. 
C 
C  Changes : ( Last change first ) 
C 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used for user input 
C            of ETP phase in ; removed code asociated with NEWETP and user supplied 
C            phase in 
C  16 May 00 @EPA- BG    bug#88 changed ETPRVP from INTEGER to REAL 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  23 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C    argument list: IEVTYP,KEYEQ,ETPRVP,ETPTMP 
C    /ETPPAR/ ETP 
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
      IMPLICIT NONE  
      INCLUDE 'ETPPAR.I' 
C 
      REAL, INTENT(IN) :: ETPRVP 
      REAL, INTENT(IN) :: ETPTMP 
C 
      INTEGER IRVP1, ITEMP1, IRVP2, ITEMP2 
      INTEGER, PARAMETER :: MAXRVP=5 
      INTEGER, PARAMETER :: MAXTMP=4 
C 
      REAL DR, A, B, DT 
      REAL, PARAMETER :: RBASE= 6.0 
      REAL, PARAMETER :: TBASE=65.0 
      REAL, PARAMETER :: RINC = 1.0 
      REAL, PARAMETER :: TINC =15.0 
C 
      REAL, DIMENSION(5,4) :: RED=RESHAPE((/ 
     *          0., 0.,  0.,  0.,  0., 
     *         15., 15., 15., 15., 15., 
     *         30., 30., 30., 30., 30., 
     *         30., 30., 30., 30., 30./),(/5,4/)) 
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
      RETURN 
      END 
