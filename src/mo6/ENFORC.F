      REAL FUNCTION ENFORC(CMPLRT,IMATP) 
C 
C  This REAL-VALUED function computes the enforcement rate given the 
C  corresponding compliance rate for either I/M or ATP. 
C 
C  Called by ATPAER and PCLEFT. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  22 Jan  1999 @DynTel-zk 2-000 Explicit Typing 
C  22 Sept 1998 @DynTel-zk 2-000 Explicit Typing 
C 
C  Input on Call: 
C 
C    argument list: CMPLRT,IMATP 
C 
C  Output on return: 
C 
C    function: ENFORC 
C 
C  Local array subscripts: 
C 
C  A(7,2)      -  A ( IRNG, IMATP ) 
C  B(7,2)      -  B ( IRNG, IMATP ) 
C  BL_RANGE(8,2)  -  BL_RANGE ( IRNG, IMATP ) 
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
C  IRNG     I    index into the A, B, and BL_RANGE arrays 
C  NCMPRT   R    non-compliance rate 
C  BL_RANGE R    range limits for applying the proper A and B terms 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
C 
      INTEGER IMATP, IRNG 
      REAL NCMPRT, CMPLRT, FBL 
      REAL, DIMENSION(7,2) :: A=RESHAPE((/ 
     *   2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0,   !I/M 
     *   2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0/), !ATP 
     *   (/7,2/)) 
C 
      REAL, DIMENSION(7,2) :: B=RESHAPE((/ 
     *   .000, .100, .175, .245, .310, .370, .425,  !I/M 
     *   .000, .100, .175, .245, .310, .370, .425/), !ATP 
     *   (/7,2/)) 
C 
      REAL, DIMENSION(8,2) :: BL_RANGE=RESHAPE((/ 
     *   .00, .05, .10, .15, .20, .25, .30, 999.,   !I/M 
     *   .00, .05, .10, .15, .20, .25, .30, 999./), !ATP 
     *   (/8,2/)) 
C 
C 
      NCMPRT=1.-CMPLRT/100. 
C 
      FBL=0. 
C 
      DO 10 IRNG=1,7 
C 
      IF(NCMPRT.LT.BL_RANGE(IRNG,IMATP).OR.NCMPRT.GE.BL_RANGE(IRNG+1, 
     *IMATP)) GOTO 10 
C 
      FBL=A(IRNG,IMATP)*(NCMPRT-BL_RANGE(IRNG,IMATP))+B(IRNG,IMATP) 
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
