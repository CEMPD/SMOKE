      REAL FUNCTION EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,INH,IFDS) 
C 
C  EF_CASE returns the basic emission factor averaged over fuel delivery 
C  system but uncorrected for operating mode, CO temperature offset, 
C  methane offset, temperature, inspection / maintenance, tampering 
C  offset and reformulated fuels. 
C 
C  Called by HRLOOP. 
C 
C  Calls EF_HDV, EF_MCDSLPRE81, EF81_93, EF_POST81 and EF_TIER2. 
C 
C  Changes: (Last change first) 
C 
C  20 Sep 00 Bug262, Emission standards for diesels are no longer averaged. 
C  07 Jul 00 Removed unused JDX from EF_POST93 call. 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C  12 May 00 @EPA-djb Initialized EF_CASE, skipped Tier 1 NOx 
C            standard effect for LDDV/T and added ISR check to HDV call. 
C            Tier 2 standards now revert to Tier1/LEV. 
C   4 Feb 00 @DynTel-ddj 1-032, Added post-1993 corrections to light 
C            duty diesel cars and trucks. 
C  21 Mar 00 @EPA-eg            hrloop revisions 
C  30 Dec 99 @DynTel-ddj 1-021, Added call to EF_TIER2 and changed 
C            intermediate values of EF_CASE to EF_DUMMY. 
C  16 Sep 99 @DynTel-bag 2-660, IGM was deleted from the EF_POST94 
C            argument list 
C  27 Aug 99 @DynTel-bag 2-660, Renamed to EF_DSPTCH.FOR. Arranged calls to 
C            emissions routines by a CASE statement dependent on vehicle 
C            types and IF statements by model years.  
C  22 Apr 99 @DynTel-ddj 2-000, Revised list of "Called by" routines. 
C  11 Mar 99 @Dyntel-RJD 2-694 Added call to EF_POST81 for 1994+ LDGV/T 
C            basic emission factors. 
C  18 Feb 99 @Dyntel-RJD 2-663 Added an additional argument to EF81_93 
C   6 JAN 99 @Dyntel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust 
C             engine start and running 
C  14 Sep 98 @Dyntel-RJD 2-660 Introduced the new ISR, INH and IFDS 
C            dependencies for light duty gasoline trucks. 
C  09 September 1998 &DynTel-ZK 2-000 Explicit Typing 
C  05 Aug 98 @Dyntel-RJD 2-663 Moved averaging over the fuel delivery 
C            systems to BEFWGT. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,JDX,IP,IV,ISR,VMTAGE,ISTD,INH,IFDS 
C    common blocks: 
C    /IVPCOM/ IVPTRS 
C    /VVNAMS/ VVSTARTS  (IVTYPE.I) 
C    /FLAGS5/ TIER2_FLAG 
C    /LDGBER/ 
C    /TIER2/  T2CRT50K(5,3,25) 
C 
C  Output on return: 
C 
C    function: EF_CASE 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C EF_DUMMY  R    Intermediate variable for basic emission values. 
C 
C 
C  Notes: 
C  15 Apr 98 @Dyntel-RJD 2-650 Modularized EF_DSPTCH for Mobile 6 so that 
C            basic emission factors are now calculated in EF_HDV, EF_MCDSLPRE81_S 
C            and EF_MCDSLPRE81_R. The new argument ISR instructs EF_DSPTCH to return 
C            either running/ftp (ISR=1) or starts (ISR=2) emissions. 
C  09 Jul 98 @Dyntel-MMS 2-666 Expanded separation of start and running 
C            modes using VVSTARTS. 
C  14 Jul 98 @Dyntel-RJD 2-660 Added mobile 6 emission factors for 81-93 
C            LDGV's. 
C  04 Aug 98 @Dyntel-RJD 2-663 Added new INH parameter to argument list. 
C            For the particular case of pre 93 LDGVs EF_DSPTCH returns a 
C            normal emission factor if INH=1 or a high emission factor 
C            if INH=2. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LDGBER.I' 
      INCLUDE 'TIER2.I' 
C 
C     declare external functions 
C 
      REAL EF_HDV 
      REAL EF_MCDSLPRE81 
      REAL EF81_93 
      REAL EF_POST93 
C 
C     declare parameter list 
C 
      INTEGER, INTENT(IN)    :: MY 
      INTEGER, INTENT(IN)    :: JDX 
      INTEGER, INTENT(IN)    :: IP 
      INTEGER, INTENT(IN)    :: IV 
      INTEGER, INTENT(IN)    :: ISR 
      INTEGER, INTENT(IN)    :: ISTD 
      INTEGER, INTENT(IN)    :: INH 
      INTEGER, INTENT(IN)    :: IFDS 
C 
      REAL, INTENT(IN)       :: VMTAGE 
C 
      REAL    TEMP_CRT50K 
      REAL    EF_DUMMY 
      REAL    EF_TIER2 
C 
C     declare local variables 
C 
      INTEGER  IVTEMP 
      INTEGER  IV50K 
      INTEGER  VEH_CLASS 
C 
      INTEGER, PARAMETER :: VCLASS0=0,VCLASS1=1,VCLASS2=2 
      INTEGER, PARAMETER :: TIER0=1, TIER2=4 
C 
C  Initialize 
C 
      IVTEMP=1 
      EF_CASE  = 0.0 
      EF_DUMMY = 0.0 
C 
C     choose emissions by vehicle classes 
C     vclass0 = vvstarts(iv)[=0] + vvldgas(iv)[=0] = 0  (heavy duty vehicles; no model year dependence) 
C     vclass1 = vvstarts(iv)[=1] + vvldgas(iv)[=0] = 1  (motorcycles and light duty diesel) 
C     vclass2 = vvstarts(iv)[=1] + vvldgas(iv)[=1] = 2  (light duty cars and trucks; model year dependence) 
C 
      VEH_CLASS = VVSTARTS(IV) + VVLDGAS(IV) 
C 
      SELECT CASE (VEH_CLASS) 
C 
C     LIGHT DUTY GAS VEHICLES 
      CASE(VCLASS2) 
C 
        IF (MY.LT.1981) THEN       ! MY < 1981 
          EF_DUMMY = EF_MCDSLPRE81(MY,IP,IV,VMTAGE,ISR) 
C 
        ELSE IF( (MY.GE.1981) .AND.(MY.LE.1993) ) THEN    ! 1980 < MY < 1994 
          EF_DUMMY = EF81_93(MY,IP,IV,VMTAGE,IFDS,ISR,INH) 
C 
        ELSE IF (MY.GE.1994) THEN  ! MY > 1993 
          EF_DUMMY = EF_POST93(MY,IP,IV,VMTAGE,ISR,INH,ISTD) 
C 
        ENDIF 
C 
C     Vehicles with separated start and running emissions but no fuel delivery 
C     system technology dependence 
C     LIGHT DUTY DIESEL VEHICLES and MOTORCYCLES 
      CASE(VCLASS1) 
        EF_DUMMY = EF_MCDSLPRE81(MY,IP,IV,VMTAGE,ISR) 
C 
        IV50K=1 
        IF(VVLDDSL(IV) .EQ. 1 .AND. MY .GE. 1994 ) THEN 
C 
C         Adjust emissions for post-Tier0 light duty diesels. 
C         LDDT12 use ratio of only LDT2 standards. 
C         LDDT34 use ratio of only LDT4 standards. 
C 
          IV50K = IVPTRS(IV) 
C 
C         Non-Tier 2 Standards (i.e., Tier 1 & LEV) 
C 
          IF(ISTD .NE. TIER2 .OR. TIER2_FLAG .EQ. 0) THEN 
            TEMP_CRT50K = CRT050K(ISTD,IP,IV50K) 
C 
C         Tier 2 Standards apply to LDDV/T. 
C 
          ELSE 
            IVTEMP = IVPTRC(IV) ! == 1,LDDV, 3 LDDT12, 5 LDDT34 
            TEMP_CRT50K = T2CRT50K(IVTEMP, IP, JDX) 
          END IF 
          EF_DUMMY = EF_DUMMY * TEMP_CRT50K / CRT050K(TIER0,IP,IV50K) 
        END IF 
C 
C     There are only HDV emissions for running (i.e., no start emissions). 
C     HEAVY DUTY VEHICLES 
      CASE(VCLASS0) 
        EF_DUMMY=0.0 
        IF(ISR.EQ.1) EF_DUMMY=EF_HDV(MY,IP,IV,VMTAGE) 
C 
      END SELECT 
C 
C     For Tier2 standards, obtain Tier2 emissions, 
C 
      IF(TIER2_FLAG .EQ. 1 .AND. MY .GT. 2003 .AND. VVLDGAS(IV). EQ. 1 
     &    .AND. ISTD .EQ. TIER2) THEN 
        EF_CASE = EF_TIER2(IV,IP,MY,JDX,INH,EF_DUMMY) 
      ELSE  ! assign non-Tier2 emission value. 
        EF_CASE = EF_DUMMY 
      END IF 
C 
      RETURN 
      END 
