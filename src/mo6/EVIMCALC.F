      SUBROUTINE EVIMCALC(MY,IMPGM,EVIMPGM,JDX) 
C 
C  This subroutine computes the I/M Evap Weighting, without and with OBD. 
C 
C  Called by EVAPCALC. 
C 
C  Calls None. 
C 
C  Changes: (Last change first) 
C 
C  15 Sep 00 AIR Task 08: New routine. 
C 
C  
C  Input on call: 
C 
C    argument list : MY, IMPGM, EVIMPG, JDX 
C    common blocks: 
C     /CHKLST/ IMCHECK 
C     /EVIM/   EVNOIM, FAIL_IL, FAIL_AZFP, FAIL_AZGC 
C     /IMPAR1/ EVFBEN 
C     /IMPAR2/ ITEST 
C     /LDGIMC/ WAVRDC 
C     /LDGOBD/ WMIL, OBDF 
C     /MAXIMA/ MAXYRS 
C 
C  Output on return: 
C 
C    argument list : None 
C    common blocks: 
C     /EVIM/ EVIM 
C 
C  Local variable / array dictionary:] 
C 
C   Name     Type                      Description 
C  ------    ----  ---------------------------------------------------------- 
C  ADJIM      R    I/M Adjustment (see M6IMFIXC.XLS spreadsheet for element origins.) 
C  AGE        R    Vehicle Age with fractional part 
C  APPL       R    Applicability rates, 1=fill pipe, 2=gas cap, 3=purge, 4=OBD. 
C  COMBO      R    Combination of fraction repaired and case 
C  COMP       R    Compliance rates, 1=fill pipe -or- OBD, 2=gas cap, 3=purge 
C  EVIMPGM    I    Applicable evap I/M program number 
C  FBEN       R    Fractional benefit, 1=fill pipe -or- OBD, 2=gas cap, 3=purge 
C  FCASE      R    Fraction repaired by test case  
C  FPREREPAIR R    Fraction pre-repaired 
C  FREPAIRED  R    Fraction repaired 
C  I          I    Temporary index 
C  IAGE       I    Vehicle Age, integer portion only 
C  IMPGM      I    Applicable exhaust I/M program number 
C  IOBD       I    1=no OBD, 2=with OBD 
C  IPF        I    Index for 1=fail pressure, 2=fail either, 3=fail both, 4=leaker 
C  J          I    Temporary index 
C  JDX        I    MOBILE6 JDX 
C  MY         I    Model year 
C  SUMPROD    R    Sum of the products of case and combo 
C  TEST       R    Testability rate 
C  WAIVE      R    Waiver rates, 1=fill pipe -or- OBD, 2=gas cap, 3=purge 
C  WBEN       R    Waiver benefit, 1=fill pipe -or- OBD, 2=gas cap, 3=purge 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BDEVIM.I' 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'EVIM.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'MAXIMA.I' 
C 
      INTEGER EVIMPGM 
      INTEGER I 
      INTEGER IMPGM 
      INTEGER IOBD 
      INTEGER IPF 
      INTEGER J 
      INTEGER JDX 
      INTEGER MY 
      REAL    APPL(4) 
      REAL    COMBO(10,8) 
      REAL    COMP(3) 
      REAL    FBEN(3) 
      REAL    FCASE(8) 
      REAL    FPREREPAIR(10) 
      REAL    FREPAIRED(3) 
      REAL    SUMPROD(10) 
      REAL    TEST(4) 
      REAL    WAIVE(3) 
      REAL    WBEN(3) 
C 
C  If Exh I/M is active, assign those values if no entered evap I/M values 
C  are available. 
C 
        IF(IMPGM.NE.0) THEN 
C 
        IF(IMCHECK(11,EVIMPGM).EQ.1.AND.IMCHECK(11,IMPGM).EQ.2) 
     *    CRIM(EVIMPGM)=CRIM(IMPGM) 
C 
        IF(IMCHECK(12,EVIMPGM).EQ.1.AND.IMCHECK(12,IMPGM).EQ.2) THEN 
          WAVPRE81(EVIMPGM)=WAVPRE81(IMPGM) 
          WAV81PLUS(EVIMPGM)=WAV81PLUS(IMPGM) 
        ENDIF 
C 
        IF(IMCHECK(25,EVIMPGM).EQ.1.AND.IMCHECK(25,IMPGM).EQ.2)  
     *    IMGRACE(EVIMPGM)=IMGRACE(IMPGM) 
C 
        IF(IMCHECK(26,EVIMPGM).EQ.1.AND.IMCHECK(26,IMPGM).EQ.2) THEN 
          DO IPF=1,4 
            EVFBEN(IPF,EVIMPGM)=EVFBEN(IPF,IMPGM) 
          END DO 
        END IF 
C 
        IF(IMCHECK(27,EVIMPGM).EQ.1) WAVRDC(EVIMPGM)=WAVRDC(IMPGM) 
C 
        IF(IMCHECK(28,EVIMPGM).EQ.1) EXEMPTAGE(EVIMPGM)=EXEMPTAGE(IMPGM) 
C 
      ENDIF     
C 
C  Initialize the 4 I/M applicability rates to zero 
C 
        APPL=0. 
C 
C  If the Evap test includes Gas Cap, set the Applicability to 100% 
C 
        IF(ITEST(EVIMPGM).GE.21) APPL(2)=1. 
C 
C  If the Evap test includes Fill Pipe, set the Applicability to 100% 
C 
        IF(ITEST(EVIMPGM).EQ.22) APPL(1)=1. 
C 
C  If Evap tests include OBD, set the Applicability to 100% 
C 
        IF(ITEST(EVIMPGM).LE.21) APPL(4)=1. 
C 
C  Set up the variables/arrays for waiver, compliance, waiver benefit and 
C  fractional benefit. 
C  No OBD Order is: Fill Pipe, Gas Cap, Purge. 
C  OBD Order is: OBD, Gas Cap, Purge 
C  Programming Note: Currently, the waiver, compliance and waiver benefit 
C  values for all test types are equal. However, like fractional benefits, 
C  these value may become test dependent. As a result, the code is already 
C  able to handle by-test type values, hence the use of arrays. Furthermore, 
C  although not used at this time, the code can also handle purge tests. 
C 
        IF(MY.LT.1981) THEN 
          WAIVE(1)=WAVPRE81(EVIMPGM) 
        ELSE 
          WAIVE(1)=WAV81PLUS(EVIMPGM) 
        ENDIF 
        WAIVE(2)=WAIVE(1) 
        WAIVE(3)=WAIVE(1) 
C 
        COMP(1)=CRIM(EVIMPGM)/100. 
        COMP(2)=COMP(1) 
        COMP(3)=COMP(1) 
C 
        WBEN(1)=1.0-WAVRDC(EVIMPGM) 
        WBEN(2)=WBEN(1) 
        WBEN(3)=WBEN(1) 
C 
        IF(ITEST(EVIMPGM).GE.22) THEN 
          FBEN(1)=EVFBEN(1,EVIMPGM)   ! Without OBD 
        ELSE 
          FBEN(1)=EVFBEN(4,EVIMPGM)   ! With OBD 
        ENDIF 
        DO IPF=2,3 
          FBEN(IPF)=EVFBEN(IPF,EVIMPGM) 
        END DO 
C 
C  Loop for without and with OBD 
C 
      DO IOBD=1,2 
C 
C  Compute the fraction repaired. 
C  No OBD order is: FILL PIPE, GAS CAP, PURGE 
C  OBD order is: OBD, GAS CAP, PURGE 
C 
      DO I=1,3 
        FREPAIRED(I)=FBEN(I)*((1.0-WAIVE(I))+WAIVE(I)*WBEN(I)) 
        IF(I.EQ.1.AND.IOBD.EQ.2) FREPAIRED(I)=FREPAIRED(I)*OBDIM 
      END DO 
C 
C  Compute the testabioity values. 
C  No OBD order is: FILL PIPE, GAS CAP, PURGE 
C  OBD order is: OBD, GAS CAP, PURGE, OBD/P OV 
C 
      IF(IOBD.EQ.1) THEN 
        I=MY-MYB_FP+1 
        IF(I.LT.1) I=1 
        IF(I.GT.18) I=18 
        TEST(1)=TEST_FP(I) 
      ELSE 
        TEST(1)=WMIL 
      ENDIF 
C 
      I=MY-MYB_GC+1 
      IF(I.LT.1) I=1 
      IF(I.GT.16) I=16 
      TEST(2)=TEST_GC(I) 
C 
      IF(IOBD.EQ.1) THEN 
        TEST(3)=TEST(1) 
      ELSE 
        I=MY-MYB_FP+1 
        IF(I.LT.1) I=1 
        IF(I.GT.18) I=18 
        TEST(3)=TEST_FP(I) 
      ENDIF 
C 
      TEST(4)=0. 
      IF(IOBD.EQ.2) TEST(4)=TEST(1)*TEST(3) 
C 
C  Compute the pre-repair fractions. Order is: 
C  PP-PPP, FP-PPP, FP-FPP, FP-FFP, FP-PFP, PF-PPF, FF-PPF, FF-FPF,FF-FFF, FF-PFF 
C 
      FPREREPAIR( 1)=EVNOIM(1,JDX,IOBD) 
C 
      FPREREPAIR( 2)=EVNOIM(3,JDX,IOBD)-((EVNOIM(3,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD)))* 
     *     (FAIL_AZFP(IOBD,JDX)+FAIL_IL(IOBD,JDX)*FAIL_AZGC(IOBD,JDX))) 
C 
      FPREREPAIR( 3)=EVNOIM(3,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_IL(IOBD,JDX)*FAIL_AZGC(IOBD,JDX)) 
C 
      FPREREPAIR( 4)=EVNOIM(3,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_IL(IOBD,JDX)*(1.0-FAIL_AZGC(IOBD,JDX))) 
C 
      FPREREPAIR( 5)=EVNOIM(3,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_AZFP(IOBD,JDX)-(FAIL_IL(IOBD,JDX)* 
     *     (1.0-FAIL_AZGC(IOBD,JDX)))) 
C 
      FPREREPAIR( 6)=EVNOIM(2,JDX,IOBD) 
C 
      FPREREPAIR( 7)=EVNOIM(4,JDX,IOBD)-((EVNOIM(4,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD)))* 
     *     (FAIL_AZFP(IOBD,JDX)+FAIL_IL(IOBD,JDX)*FAIL_AZGC(IOBD,JDX))) 
C 
      FPREREPAIR( 8)=EVNOIM(4,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_IL(IOBD,JDX)*FAIL_AZGC(IOBD,JDX)) 
C 
      FPREREPAIR( 9)=EVNOIM(4,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_IL(IOBD,JDX)*(1.0-FAIL_AZGC(IOBD,JDX))) 
C 
      FPREREPAIR(10)=EVNOIM(4,JDX,IOBD)/ 
     *     (EVNOIM(3,JDX,IOBD)+EVNOIM(4,JDX,IOBD))* 
     *     (FAIL_AZFP(IOBD,JDX)-FAIL_IL(IOBD,JDX)* 
     *     (1.0-FAIL_AZGC(IOBD,JDX))) 
C 
C  For no OBD: 
C  Compute the fractions getting repair due to test combinations, order is: 
C  GC+FP+PURGE, GC+FP, GC+P, FP+P, GC, P, FP, NONE 
C 
      IF(IOBD.EQ.1) THEN 
C 
      FCASE(1)=    MIN(     APPL(1),     APPL(2),     APPL(3))* 
     *             MIN(     TEST(1),     TEST(2),     TEST(3))* 
     *             MIN(     COMP(1),     COMP(2),     COMP(3))* 
     *             MIN(FREPAIRED(1),FREPAIRED(2),FREPAIRED(3)) 
C 
      FCASE(2)=MAX(MIN(     APPL(1),     APPL(2))* 
     *             MIN(     TEST(1),     TEST(2))* 
     *             MIN(     COMP(1),     COMP(2))* 
     *             MIN(FREPAIRED(1),FREPAIRED(2))- 
     *                 FCASE(1),0.) 
C 
      FCASE(3)=MAX(MIN(     APPL(2),     APPL(3))* 
     *             MIN(     TEST(2),     TEST(3))* 
     *             MIN(     COMP(2),     COMP(3))* 
     *             MIN(FREPAIRED(2),FREPAIRED(3))- 
     *                 FCASE(1),0.) 
C 
      FCASE(4)=MAX(MIN(     APPL(1),     APPL(2),     APPL(3))* 
     *             MIN(     TEST(1),     TEST(2),     TEST(3))* 
     *             MIN(     COMP(1),     COMP(2),     COMP(3))* 
     *             MIN(FREPAIRED(1),FREPAIRED(2),FREPAIRED(3))- 
     *                 FCASE(1),0.) 
C 
      FCASE(5)=MAX((APPL(2)*TEST(2)*COMP(2)*FREPAIRED(2))- 
     *             (FCASE(1)+FCASE(2)+FCASE(3)),0.) 
C 
      FCASE(6)=MAX((APPL(3)*TEST(3)*COMP(3)*FREPAIRED(3))- 
     *             (FCASE(1)+FCASE(3)+FCASE(4)),0.) 
C 
      FCASE(7)=MAX((APPL(1)*TEST(1)*COMP(1)*FREPAIRED(1))- 
     *             (FCASE(1)+FCASE(2)+FCASE(4)),0.) 
C 
      FCASE(8)=1.0-(FCASE(1)+FCASE(2)+FCASE(3)+FCASE(4)+ 
     *              FCASE(5)+FCASE(6)+FCASE(7)) 
C 
C  For OBD: 
C  Compute the fractions getting repair due to test combinations, order is: 
C  OBD+GC+PURGE, OBD+GC, GC+P, OBD+P, GC, P, OBD, NONE 
C 
      ELSE 
C 
      FCASE(1)=    MIN(     APPL(4),     APPL(2),     APPL(3))* 
     *             MIN(     TEST(2),     TEST(4)             )* 
     *             MIN(     COMP(1),     COMP(2),     COMP(3))* 
     *             MIN(FREPAIRED(1),FREPAIRED(2),FREPAIRED(3)) 
C 
      FCASE(2)=MAX(MIN(     APPL(4),     APPL(2))* 
     *             MIN(     TEST(1),     TEST(2))* 
     *             MIN(     COMP(1),     COMP(2))* 
     *             MIN(FREPAIRED(1),FREPAIRED(2))- 
     *                 FCASE(1),0.) 
C 
      FCASE(3)=MAX(MIN(     APPL(2),     APPL(3))* 
     *             MIN(     TEST(2),     TEST(3))* 
     *             MIN(     COMP(2),     COMP(3))* 
     *             MIN(FREPAIRED(2),FREPAIRED(3))- 
     *                 FCASE(1),0.) 
C 
      FCASE(4)=MAX(MIN(     APPL(4),     APPL(3))* 
     *                      TEST(4)              * 
     *             MIN(     COMP(1),     COMP(3))* 
     *             MIN(       OBDIM,FREPAIRED(3))- 
     *                 FCASE(1),0.) 
C 
      FCASE(5)=MAX((APPL(2)*TEST(2)*COMP(2)*FREPAIRED(2))- 
     *             (FCASE(1)+FCASE(2)+FCASE(3)+FCASE(4)),0.) 
C 
      FCASE(6)=MAX((APPL(3)*TEST(3)*COMP(3)*FREPAIRED(3))- 
     *             (FCASE(1)+FCASE(3)+FCASE(4)),0.) 
C 
      FCASE(7)=MAX((APPL(4)*TEST(1)*COMP(1)*OBDIM)- 
     *             (FCASE(1)+FCASE(2)+FCASE(4)),0.) 
C 
      FCASE(8)=1.0-(FCASE(1)+FCASE(2)+FCASE(3)+FCASE(4)+ 
     *              FCASE(5)+FCASE(6)+FCASE(7)) 
C 
      ENDIF 
C 
C  Compute the by-case by-pre-repair combination fractions 
C 
      COMBO=0. 
C 
C  For no OBD: 
C 
      IF(IOBD.EQ.1) THEN 
C  PP-PPP 
      COMBO(1,1)=FPREREPAIR(1)+FPREREPAIR(3)+FPREREPAIR(4)+ 
     *           FPREREPAIR(5)+FPREREPAIR(6)+FPREREPAIR(8)+ 
     *           FPREREPAIR(9)+FPREREPAIR(10) 
      COMBO(1,2)=FPREREPAIR(1)+FPREREPAIR(3)+FPREREPAIR(4)+ 
     *           FPREREPAIR(5) 
      COMBO(1,3)=FPREREPAIR(1)+FPREREPAIR(3)+FPREREPAIR(6)+ 
     *           FPREREPAIR(8) 
      COMBO(1,4)=FPREREPAIR(1)+FPREREPAIR(5)+FPREREPAIR(6)+ 
     *           FPREREPAIR(10) 
      COMBO(1,5)=FPREREPAIR(1)+FPREREPAIR(3) 
      COMBO(1,6)=FPREREPAIR(1)+FPREREPAIR(6) 
      COMBO(1,7)=FPREREPAIR(1)+FPREREPAIR(5) 
      COMBO(1,8)=FPREREPAIR(1) 
C  FP-PPP 
      COMBO(2,1)=FPREREPAIR(2)+FPREREPAIR(7) 
      COMBO(2,2)=FPREREPAIR(2) 
      COMBO(2,3)=FPREREPAIR(2)+FPREREPAIR(7) 
      COMBO(2,4)=FPREREPAIR(2)+FPREREPAIR(7) 
      COMBO(2,5)=FPREREPAIR(2) 
      COMBO(2,6)=FPREREPAIR(2)+FPREREPAIR(7) 
      COMBO(2,7)=FPREREPAIR(2) 
      COMBO(2,8)=FPREREPAIR(2) 
C  FP-FPP 
      COMBO(3,4)=FPREREPAIR(3)+FPREREPAIR(4)+FPREREPAIR(8)+ 
     *           FPREREPAIR(9) 
      COMBO(3,6)=FPREREPAIR(3)+FPREREPAIR(8) 
      COMBO(3,7)=FPREREPAIR(3)+FPREREPAIR(4) 
      COMBO(3,8)=FPREREPAIR(3) 
C  FP-FFP      
      COMBO(4,6)=FPREREPAIR(4)+FPREREPAIR(9) 
      COMBO(4,8)=FPREREPAIR(4) 
C  FP-PFP 
      COMBO(5,3)=FPREREPAIR(5)+FPREREPAIR(4)+FPREREPAIR(9)+ 
     *           FPREREPAIR(10) 
      COMBO(5,5)=FPREREPAIR(4)+FPREREPAIR(5) 
      COMBO(5,6)=FPREREPAIR(5)+FPREREPAIR(10) 
      COMBO(5,8)=FPREREPAIR(5) 
C  PF-PPF 
      COMBO(6,2)=FPREREPAIR(6)+FPREREPAIR(8)+FPREREPAIR(9)+ 
     *           FPREREPAIR(10) 
      COMBO(6,5)=FPREREPAIR(6)+FPREREPAIR(8) 
      COMBO(6,7)=FPREREPAIR(6)+FPREREPAIR(10) 
      COMBO(6,8)=FPREREPAIR(6) 
C  FF-PPF 
      COMBO(7,2)=FPREREPAIR(7) 
      COMBO(7,5)=FPREREPAIR(7) 
      COMBO(7,7)=FPREREPAIR(7) 
      COMBO(7,8)=FPREREPAIR(7) 
C  FF-FPF 
      COMBO(8,7)=FPREREPAIR(8)+FPREREPAIR(9) 
      COMBO(8,8)=FPREREPAIR(8) 
C  FF-FFF 
      COMBO(9,8)=FPREREPAIR(9) 
C  FF-PFF 
      COMBO(10,5)=FPREREPAIR(9)+FPREREPAIR(10) 
      COMBO(10,8)=FPREREPAIR(10) 
C 
C  For OBD: 
C 
      ELSE 
C 
C  PP-PPP 
      COMBO(1,1)=FPREREPAIR(1)+FPREREPAIR(2)+FPREREPAIR(3)+ 
     *           FPREREPAIR(4)+FPREREPAIR(5)+FPREREPAIR(6)+ 
     *           FPREREPAIR(7)+FPREREPAIR(8)+FPREREPAIR(9)+ 
     *           FPREREPAIR(10) 
      COMBO(1,2)=COMBO(1,1) 
      COMBO(1,3)=FPREREPAIR(1)+FPREREPAIR(3)+FPREREPAIR(6)+ 
     *           FPREREPAIR(8) 
      COMBO(1,4)=COMBO(1,1) 
      COMBO(1,5)=FPREREPAIR(1)+FPREREPAIR(3) 
      COMBO(1,6)=FPREREPAIR(1)+FPREREPAIR(6) 
      COMBO(1,7)=COMBO(1,1) 
      COMBO(1,8)=FPREREPAIR(1) 
C  FP-PPP 
      COMBO(2,3)=FPREREPAIR(2)+FPREREPAIR(7) 
      COMBO(2,5)=FPREREPAIR(2) 
      COMBO(2,6)=COMBO(2,3) 
      COMBO(2,8)=FPREREPAIR(2) 
C  FP-FPP 
      COMBO(3,6)=FPREREPAIR(3)+FPREREPAIR(8) 
      COMBO(3,8)=FPREREPAIR(3) 
C  FP-FFP      
      COMBO(4,6)=FPREREPAIR(4)+FPREREPAIR(9) 
      COMBO(4,8)=FPREREPAIR(4) 
C  FP-PFP 
      COMBO(5,3)=FPREREPAIR(4)+FPREREPAIR(5)+FPREREPAIR(9)+ 
     *           FPREREPAIR(10) 
      COMBO(5,5)=FPREREPAIR(4)+FPREREPAIR(5) 
      COMBO(5,6)=FPREREPAIR(5)+FPREREPAIR(10) 
      COMBO(5,8)=FPREREPAIR(5) 
C  PF-PPF 
      COMBO(6,5)=FPREREPAIR(6)+FPREREPAIR(8) 
      COMBO(6,8)=FPREREPAIR(6) 
C  FF-PPF 
      COMBO(7,5)=FPREREPAIR(7) 
      COMBO(7,8)=FPREREPAIR(7) 
C  FF-FPF 
      COMBO(8,8)=FPREREPAIR(8) 
C  FF-FFF 
      COMBO(9,8)=FPREREPAIR(9) 
C  FF-PFF 
      COMBO(10,5)=FPREREPAIR(9)+FPREREPAIR(10) 
      COMBO(10,8)=FPREREPAIR(10) 
C 
      ENDIF 
C 
C  Compute the sum of the products of the case and combos 
C 
      DO I=1,10 
        SUMPROD(I)=0. 
        DO J=1,8 
          SUMPROD(I)=SUMPROD(I)+FCASE(J)*COMBO(I,J) 
        END DO 
      END DO 
C 
C  Compute the with-im weightings 
C  Order is: PP PF FP FF 
C 
      EVIM(1,JDX,IOBD)=SUMPROD(1) 
      EVIM(2,JDX,IOBD)=SUMPROD(6) 
      EVIM(3,JDX,IOBD)=SUMPROD(2)+SUMPROD(3)+SUMPROD(4)+SUMPROD(5) 
      EVIM(4,JDX,IOBD)=SUMPROD(7)+SUMPROD(8)+SUMPROD(9)+SUMPROD(10) 
C 
C  For age = 0 (JDX=1), set the with I/M values equal to the no I/M values    
C 
      IF(JDX.EQ.1) THEN 
        DO I=1,4 
          EVIM(I,1,IOBD)=EVNOIM(I,1,IOBD) 
        END DO 
      END IF 
C 
      END DO 
C 
      RETURN 
      END 
       
