      REAL FUNCTION 
     *     IDRATE(MY,JDX,IP,IV,ISR,IFDS,ISTD,IMPGM,VMTAGE,INERR) 
C 
C  IDRATE returns the high emitter identification rate for different I/M 
C  test types. (The procedures in this code are in EPA report M6.IM.001 pages 
C  9-11 and 26&27 for I&M 240 run and start fractions, respectively. Pages 
C  57-60 for ASM tests, 61-66 for idle tests, and in the last paragraph of 
C  page 67 for the OBD test.) 
C 
C  Called by IMTEST. 
C 
C  Calls EF_HI_FRCN 
C 
C  Changes: 
C 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C  02 Nov 00 @EPA-djb Bug299, initialized IDRATE. 
C  06 Sep 00 @EPA-elg Eliminate potential divide by zero condition in HI_FRCN. 
C  07 Jul 00 @EPA-elg Changed ASM IM Credit logic. 
C  08 Jun 00 EPA-DJB  bug#113 included IFDS in argument list so that it can be 
C                     passed to ef_hi_frcn 
C  22 May 00 EPA-EG       Changed idr limit to a maximum value of 1.0. 
C  25 Apr 00 AIR Task 02: Included INERR in IDRATE and EF_HI_FRCN parameter lists. 
C  09 SEP 99 @DynTel-BAG 2-660  Renamed HIEMIT.FOR to EF_HI_FRCN.FOR 
C  03 Aug 99 @DynTel-BAG 2-000  Replaced IF statements over I&M test type 
C            with CASE statements and introduced parameters for integer 
C            based test types. 
C  20 jUL 99 @DynTel-BAG 2-000  Created new variable IDRATE_MY to 
C            replace local reassignment of passed nonlocal variable MY 
C   4 May 99 @DynTel-ddj 2-694, Changed assignment of HI4 from 
C            HIFB4(JDX,IAM,IP,IV) to HIGH(JDX,IP,IV). 
C  18 Mar 99 @Dyntel-RJD 2-694 Extended the scope of IDRATE to include 
C            1994+ model year LDGV/Ts. All vehicles are equipped with 
C            OBD after 1996 and OBD I/M programs thus become possible. 
C  02 Dec 98 @Dyntel-RJD 2-663 Created IDRATE for 1981-1993 model year 
C            light duty gasoline vehicles and trucks. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,JDX,VAGE_LSTIM,IP,IV,ISR,IFDS,ISTD,IMPGM,INERR 
C 
C    common blocks: 
C    /IMPAR2/ ITEST 
C    /IMPAR4/ CUTNO,IVCUTIDX 
C    /LDGIMC/ ASMCP1,ASMCP2,IDR240,IDRASM,IMPEFF 
C    /LDGOBD/ HIGH3 
C    /LDGSIM/ HCLOSS,COLOSS,NOXLOSS 
C    /MYRSAV/ AMAR 
C    /RSDINF/ RSDCLEAN 
C 
C  Output on return: 
C 
C    function: IDRATE 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name     Type              Description 
C  ------    ----  ----------------------------------------------------- 
C   HI_FRCN    R    Fraction of high emitters in the fleet 
C   HI3        R    Fraction of high emitters after correction for OBD I/M 
C   IAIM       I    The integer part of the age VAGE_LSTIM. 
C   ICT        I    remomte sensing clean screening effectiveness, RSDCSE(ICT,IP,IMSL), 
C                   index based on NOx cutpoint, RSDNOX(IMPGM) : 
C                   RSDNOX(IMPGM) = 0,    ICT = 1 
C                   RSDNOX(IMPGM) = 1000, ICT = 2 
C                   RSDNOX(IMPGM) = 1500, ICT = 3 
C                   RSDNOX(IMPGM) = 2000, ICT = 4 
C   IMSL       I    Flag indicating strictness of I/M program 
C   IMTT       I    IM test type. 
C   IMY        I    MY pointer for ASM Credits.  1981 = 1. 
C   RAIM       R    For modeling purposes, in selective IM programs such as 
C                   RSD and COIM vehicles are assumed to be eligible for an 
C                   IM test every six months following the date of sale of 
C                   the vehicle into the fleet. RAIM is the difference 
C                   between the age of the vehicle VAGE_LSTIM and its integer part 
C                   IAIM. RAIM is always 0 or 0.5. 
C   VMTAIM     R    Vehicle mileage in 10000 mile units at age VAGE_LSTIM. 
C   IDRATE_MY  I    local variable for vehicle model year 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'RSDINF.I' 
C 
C Declare external functions 
C 
      REAL EF_HI_FRCN 
C 
C Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: ISTD 
      INTEGER, INTENT(IN) :: IMPGM 
      INTEGER, INTENT(INOUT) :: INERR 
      REAL, INTENT(IN) :: VMTAGE 
C 
C Declare local variables 
C 
      REAL    HI_FRCN 
      REAL    HI4 
      INTEGER IMTT 
      INTEGER ICT 
      INTEGER IMSL 
      INTEGER IMY 
      INTEGER IDRATE_MY 
      INTEGER IDRATE_IFDS 
      INTEGER, PARAMETER :: HC=1,CO=2,NOX=3 
      INTEGER, PARAMETER :: IDL_TST=1, IDL25_TST=2 
      INTEGER, PARAMETER :: IM240_TST=4, ASM25P_TST=6 
      INTEGER, PARAMETER :: ASM2550F_TST=11, OBD_TST=12 
      INTEGER, PARAMETER :: LENIENT=1, STRICT=2 
C 
C     I&M test types and integer designations (see also GETIPR.FOR) : 
C 
C       1 = idle test 
C       2 = 2500 rpm idle test 
C       3 = loaded idle test 
C       4 = IM240 test 
C       5 = ASM test procedure 
C       6 = ASM 2525 PHASE-IN 
C       7 = ASM 2525 FINAL 
C       8 = ASM 5015 PHASE-IN 
C       9 = ASM 5015 FINAL 
C      10 = ASM 2525/5015 PHASE-IN 
C      11 = ASM 2525/5015 FINAL 
C      12 = OBD I/M 
C 
      IDRATE=0.0 
      IDRATE_IFDS = IFDS 
      IF (MY.GT.1993) IDRATE_IFDS=3 
      IDRATE_MY = MY 
      IF (MY.GT.1998) IDRATE_MY=1998 
C 
C Calculate the high emitter identification rate for the different test 
C types. 
C 
      IMTT = ITEST(IMPGM) 
C 
      IMTEST_TYPE : SELECT CASE (IMTT) 
C 
C     idle test (IDRATE = 0.0 for NOx, see BD_IM.FOR and p.64 M6.IM.001) 
      CASE (IDL_TST) IMTEST_TYPE 
        IDRATE=IDRID1(IDRATE_IFDS,ISR,IP) 
C 
C     2500 RPM idle test (IDRATE = 0.0 for NOx, see BD_IM.FOR and p.64 M6.IM.001) 
      CASE (IDL25_TST) IMTEST_TYPE 
        IDRATE=IDRID2(IDRATE_IFDS,ISR,IP) 
C 
C     I&M 240 test 
      CASE (IM240_TST) IMTEST_TYPE 
        IDRATE=IDR240(JDX,ISR,IP,IVCUTIDX(IV),IMPGM) 
C 
C     ASM tests ( see pages 57-60 of M6.IM.001 ) 
      CASE (ASM25P_TST:ASM2550F_TST) IMTEST_TYPE 
        IMY = MY - 1980 
        IF(MY.GE.1998) IMY = 18 
        IDRATE=ASMRAT(IP,IMY,JDX)*IDR240(JDX,ISR,IP,IVCUTIDX(IV),IMPGM) 
C 
C     OBD test 
      CASE (OBD_TST) IMTEST_TYPE 
        HI4 = OBDHINIM(JDX,IP,IV) !fraction of high emitting vehicles with OBD and without I&M 
                                  !programs this fraction includes owner response and indicator 
                                  !light percentages see OBD_HI_FRCN.FOR 
        HI_FRCN = 
     &      EF_HI_FRCN(IDRATE_MY,IP,IV,VMTAGE,IDRATE_IFDS,ISTD,INERR) 
C 
        IF (HI_FRCN.EQ.0) THEN 
          IDRATE = 0.0 
        ELSE 
          IDRATE = 1.0-(HI4/HI_FRCN) !ratio of non-I&M OBD high emitting vehicles with age JDX to 
                                     !the total number of high emitting vehicles in model year MY 
                                     !IDRATE is the fraction of OBD high emitting vehicles 
                                     !in I&M programs in model year MY and age JDX 
        ENDIF 
C 
      END SELECT IMTEST_TYPE 
C 
C Correct the high emitter identification rate for the effect of an RSD 
C clean screening program. 
C 
      IF (RSDCLEAN(IMPGM).NE.0) THEN 
C 
C       map the RSDNOX cutpoints to an integer index for the remote sensing 
C       clean screening effectiveness fractions RSDCSE. (See M6.IM.001 table 
C       6a&b pages 48&49) 
        ICT = NINT(RSDNOX(IMPGM)/500.0) 
        IF (ICT .EQ. 0) ICT=1 
C 
C       RSD coverage option 1(see p.39 of M6.IM.001, Poisson dist.) 
C       and RSD NOx cutpoints 
        IF (RSDCLEAN(IMPGM).EQ.1) THEN 
           IMSL = IMSTRICT(IMPGM) 
           IF (IMSL.EQ.0.AND.CUTNO(JDX,IVCUTIDX(IV),IMPGM).GE.3.0) 
     *         IMSL = LENIENT  !RSDCSE from table 6a page 48 
           IF (IMSL.EQ.0.AND.CUTNO(JDX,IVCUTIDX(IV),IMPGM).LT.3.0) 
     *         IMSL = STRICT   !RSDCSE from table 6b page 49 
           IDRATE = IDRATE-RSDCVR1(JDX,IV,IMPGM)*RSDCSE(ICT,IP,IMSL) 
C 
C       RSD coverage option 2(see p.40 of M6.IM.001, user supplied) 
C       and RSD NOx cutpoints 
        ELSE IF (RSDCLEAN(IMPGM).EQ.2) THEN 
           IMSL = IMSTRICT(IMPGM) 
           IF (IMSL.EQ.0.AND.CUTNO(JDX,IVCUTIDX(IV),IMPGM).GE.3.0) 
     *         IMSL = LENIENT  !RSDCSE from table 6a page 48 
           IF (IMSL.EQ.0.AND.CUTNO(JDX,IVCUTIDX(IV),IMPGM).LT.3.0) 
     *         IMSL = STRICT   !RSDCSE from table 6b page 49 
           IDRATE = IDRATE - RSDCVR2(JDX,IMPGM) * RSDCSE(ICT,IP,IMSL) 
C 
C       RSD coverage option 3 (user supplied RSD clean screening identification 
C       rate effects/losses) Note: The M6.IM.001 mentions an option 3 which 
C       uses user supplied RSD failure and fleet size numbers to determine 
C       the RSD losses. However, this is not used to determine RSD effects on 
C       the identification rates. (See the second to the last paragraph on page 
C       47.) 
C 
        ELSE 
           IF (IP.EQ.HC) IDRATE = IDRATE - HCLOSS (IMPGM) 
           IF (IP.EQ.CO) IDRATE = IDRATE - COLOSS (IMPGM) 
           IF (IP.EQ.NOX) IDRATE = IDRATE - NOXLOSS(IMPGM) 
        END IF 
      END IF 
C 
C The efficiency parameter has been introduced to reduce the benefit 
C from programs that are not run at full efficiency. It can take any 
C value between 0 and 1. The default value is 1. 
C 
      IDRATE=IMPEFF(IMPGM,IP)*IDRATE 
C 
C The idr is limited to a maximum value of 1.0 (May 22, 2000). 
C 
      IF (IDRATE.GT.1.0) IDRATE=1.0 
C 
      RETURN 
      END 
