      SUBROUTINE EVMAIN(IH,MY,KEYEQ,IVGAS,ISTUSE,RIVAL,EVRATE,IFDS) 
C 
C  EVMAIN is the driver "main" for the block of functions that calculate 
C  hot soak given the vehicle type, emission standard, base/in-use/certification 
C  RVP/Index and the untampered/uncontrolled/cert equation group specified in 
C  the argument list. 
C 
C  Called by HOTSOK. 
C 
C  Calls EVPRE, EVSTD, EV81, EV85, EVUNC and EVHI. 
C 
C  Changes : (Last change first) 
C 
C  24 Aug 00 @EPA-djb Bug 233, Passed IFDS from HOTSOK (instead of 
C            calculating it locally.  Removed ISTD (use ISTUSE instead). 
C            Added IFDS to EV85 call. 
C  20 Jun 00 @EPA-ddj Removed IVGAS from parameter list of EVUNC. 
C   2 Jun 00 @EPA-ddj Bug 86, Expanded Computed GO TO for LDGT3/4. 
C  28 Sep  1999 @EPA-EG          Updated EVUNC function to reflect Arcadis Hot Soak Analysis. 
C   8 Mar  1999 @Dyntel-BG 2-664 added hourly index to parameter list 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  23 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C  11 Jun  1997 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C               code to use the symbolic vehicle type names and vehicle type 
C               vectors that were introduced in Mobile6. 
C  20 Sep  1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year 
C               year.  removed IY and ICY from the argument list and thus 
C               from the calls to EV81 and EVSTD. 
C 
C  Input on call: 
C 
C    argument list: MY,KEYEQ,IVGAS,ISTUSE,RIVAL,IH,EVRATE,IFDS 
C    common blocks: 
C    /EVAHIA/ HIADJ 
C    /REGION/ IREJN 
C    /VVNAMS/ VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS, VTLDGV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: EVRATE 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  EVMAX    R    Maximum limit (grams) for low altitude hot soak estimates. 
C 
C  Notes: 
C 
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed in MOBILE4.1. 
C  Modified for MOBILE4.1 to add support for pressure/purge fail. 
C 
C  Hot Soak rates are adjusted for high altitude regions (IREJN=2) as follows: 
C     Hot Soak High Alt Rate = Hot Soak Low Alt Rate * HIADJ 
C  EXCEPT: 
C        LDGVs   1972-76my   Actual LDGV Test Data 
C                1977my      Same as Low Alt 
C                1984+my     Same as Low Alt 
C        LDGT12s 1972-76my   Same as LDGVs 
C                1977my      Same as Low Alt 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAHIA.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'REGION.I' 
C 
C     Declare external functions 
C 
      REAL EVPRE 
      REAL EVSTD 
      REAL EV81 
      REAL EV85 
      REAL EVHI 
      REAL EVUNC 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: KEYEQ 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IVGAS 
      INTEGER ISTUSE 
      REAL EVRATE 
      REAL RIVAL 
C 
C     Declare local variables 
C 
      REAL EVMAX 
      INTEGER IV 
C 
C  Function call depends 
C      first on the equation group key, 
C            second on the vehicle type and 
C                   third on the emission standard. 
C 
      IV = IVGAS 
      IF (IVGAS.EQ.14) THEN 
        IV = VTMC 
      ELSEIF (IVGAS.EQ.15) THEN 
        IV = VTGASBUS 
      ENDIF 
C 
      GOTO (10,55,55), KEYEQ 
C 
C  KEYEQ = 1 = passing rates 
C 
   10 IF (IV.EQ.VTLDGV.OR.VVLDGT12(IV).EQ.1)  GOTO 15 
      IF (VVLDGT34(IV).EQ.1)                  GOTO 35 
      IF (VVHDGV(IV).EQ.1)                    GOTO 50 
      IF (IV.EQ.VTMC)                         GOTO 25 
C 
C  Passing LDGV/LDGT12 
C 
   15 GOTO (20,25,25,25,30,30,30), ISTUSE 
   20 EVRATE=EVPRE(IH,KEYEQ,ISTUSE,IVGAS,RIVAL) 
      GOTO 90 
   25 EVRATE=EVSTD(IH,KEYEQ,ISTUSE,IVGAS,RIVAL,MY) 
      IF(IVGAS.EQ.14) GOTO 100 
      IF(IREJN.EQ.1.OR.MY.EQ.1977) GOTO 100 
      IF(ISTUSE.EQ.2.OR.ISTUSE.EQ.4) GOTO 90 
      IF(ISTUSE.EQ.3) EVRATE=EVHI(KEYEQ,RIVAL) 
      GOTO 100 
   30 EVRATE=EV81(IFDS,IVGAS,RIVAL,MY) 
      GOTO 90 
C 
C   Passing LDGT34 
C 
   35 GOTO (20,40,45,45,45,30,30), ISTUSE 
   40 EVRATE=EVSTD(IH,KEYEQ,ISTUSE,IVGAS,RIVAL,MY) 
      GOTO 90 
   45 EVRATE=EV81(IFDS,IVGAS,RIVAL,MY) 
      GOTO 90 
C 
C   Passing HDGV, GAS BUS 
C 
   50 IF(ISTUSE.EQ.1) GOTO 20 
      EVRATE=EV85(IH,KEYEQ,IVGAS,RIVAL,MY,IFDS) 
      GOTO 90 
C 
C  KEYEQ = 2 = failed purge uncontrolled emissions 
C  KEYEQ = 3 = failed pressure uncontrolled emissions 
C 
   55 IF (IV.EQ.VTLDGV.OR.VVLDGT12(IV).EQ.1)  GOTO 60 
      IF (VVLDGT34(IV).EQ.1)                  GOTO 70 
      IF (VVHDGV(IV).EQ.1)                    GOTO 50 
      IF (IV.EQ.VTMC)                         GOTO 25 
C 
C  Failed LDGV/LDGT12 
C 
   60 GOTO (20,25,25,25,80,80,80), ISTUSE 
C 
C  Failed LDGT34 
C 
   70 GOTO (20,40,80,80,80), ISTUSE 
C 
C  Failed 1981+ LDGV/T1/T2 
C 
   80 EVRATE=EVUNC(IH,KEYEQ,RIVAL) 
   90 IF(IREJN.EQ.1) GOTO 100 
      IF((IVGAS.EQ.VTLDGV.OR.VVLDGT12(IVGAS).EQ.1).AND. 
     &     MY.EQ.1977) GOTO 100 
      IF(IVGAS.EQ.VTLDGV.AND.MY.GE.1984) GOTO 100 
      EVRATE=EVRATE*HIADJ 
C 
C  Limit Hot Soaks to 0-45 grams, inclusive, (low alt). 
C  Limit Diurnals to 0-100 grams, inclusive, (low alt). 
C  Account for hi-altitude accordingly. 
C 
  100 EVMAX=45.0 
C 
      IF(IREJN.EQ.2) EVMAX=EVMAX*HIADJ 
C 
      EVRATE=AMIN1(AMAX1(EVRATE,0.0),EVMAX) 
C 
      RETURN 
      END 
