C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION ENFORC(CMPLRT,IMATP)
C
C  This REAL-VALUED function computes the enforcement rate given the
C  corresponding compliance rate for either I/M or ATP.
C
C  Called by ATPAER and PCLEFT.
C
C  Input on Call:
C
C    parameter list: CMPLRT,IMATP
C
C  Output on return:
C
C    function: ENFORC
C
C  Local array subscripts:
C
C  A(7,2)      -  A ( IRNG, IMATP )
C  B(7,2)      -  B ( IRNG, IMATP )
C  RANGE(8,2)  -  RANGE ( IRNG, IMATP )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  A        R    first term of the FBL linear equation: FBL=A*X+B
C  B        R    last term of the FBL linear equation: FBL=A*X+B
C  CMPLRT   R    compliance rate
C  FBL      R    fraction of benefit loss
C  IMATP    I    flag indicating whether calculating I/M or ATP rates
C  IRNG     I    index into the A, B, and RANGE arrays
C  NCMPRT   R    non-compliance rate
C  RANGE    R    range limits for applying the proper A and B terms
C
C  Notes:
C
C  None.
C
C
      REAL NCMPRT
      REAL A(7,2),B(7,2),RANGE(8,2)
C
      DATA A/
C I/M
     *   2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0,
C ATP
     *   2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0/
C
      DATA B/
C I/M
     *   .000, .100, .175, .245, .310, .370, .425,
C ATP
     *   .000, .100, .175, .245, .310, .370, .425/
C
      DATA RANGE/
C I/M
     *   .00, .05, .10, .15, .20, .25, .30, 999.,
C ATP
     *   .00, .05, .10, .15, .20, .25, .30, 999./
C
C
      NCMPRT=1.-CMPLRT/100.
C
      FBL=0.
C
      DO 10 IRNG=1,7
C
      IF(NCMPRT.LT.RANGE(IRNG,IMATP).OR.NCMPRT.GE.RANGE(IRNG+1,
     *IMATP)) GOTO 10
C
      FBL=A(IRNG,IMATP)*(NCMPRT-RANGE(IRNG,IMATP))+B(IRNG,IMATP)
      GOTO 20
C
   10 CONTINUE
C
   20 IF(FBL.GT.1.0) FBL=1.0
C
      ENFORC=1.0-FBL
C
      RETURN
      END
