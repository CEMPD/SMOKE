      REAL FUNCTION IMPROG(ICY,MY,ACTUAL_MY,JDX,IP,IV,VMTAGE,ISR,IFDS,
     *                     ISTD,INERR,IDX,IH,IROAD,CAT,HIGHEF)
C
C  IMPROG returns the additive correction factor for an inspection and
C  maintenance program. The correction is applicable to 1981+ light
C  gasoline vehicles and trucks. The total additive IM correction factor
C  is then computed as a weighted average over all the IM test groups.
C
C  Called by HRLOOP.
C
C  Calls MIN, MAX, FINDIM, IDRATE, CLCIMPAR, AFTRPREF
C
C  Changes:
C
C
C  02 May 02 AIR Task Bug Fix: #388 modification of constants via QUITER
C   2 Jan 02 @EPA-BAG bug367 Fix double correction of repaired emissions/ remove
c                            call to IMCF
C  21 Dec 01 @EPA-BAG        Fix IDR for post 1995 vehicles
C  17 Dec 01 @EPA-BAG        Remove I&M sawtooth
C  12 Dec 01 @EPA-BAG        Removed REPAIRX in determination of AVEREF for
C                            post1995/OBD equipped vehicles
C  07 Dec 01 @EPA-BAG        Added variable initialization and fixed argument
C                            list for IDRATE
C  03 Dec 01 @EPA-BAG        Replace IMPN(ISEG) with CIMPN/Remove AGEIM and
C                            associated variables/ see also INITIM and FINDIM
C  27 Nov 01 @EPA-BAG bug331 NOx benefits for idle I&M programs
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C  09 May 01 @EPA-elg Added Exhaust I/M Factor for 1996+ Model Years.
C  03 Jan 01 @EPA-djb Bug307, Removed IMPGM loop from REDPROG calculation.
C  04 Dec 00 @EPA-bag Task X6, added ACTUAL_MY to parameter list
C  02 Nov 00 @EPA-djb Bug299, skip calculations for evap I/M cases.
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM
C  15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M.
C  29 Aug 00 @EPA-djb bug241 Removed call to IMPROG_NP
C  25 Jul 00 @EPA-djb Added methane speciation, removed non-I/M IMREMAIN
C            calculations and altered/added some comments.
C  14 Jun 00 @EPA-elg Replaced EF_POST with EF_CASE.
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list
C  08 Jun 00 @EPA-djb Bug113 Added ISTD to EF_HI_FRCN call.
C  17 May 00 @EPA-bag bug#88 replaced HI3 with its intented value, i.e.,
C                         OBDHINIM(JDX,IP,IV,) from obd_hi_frcn.for
C  24 Apr 00 AIR Task 02: Included INERR into IMPROG, IMBENEFIT, IMPROG_NP,
C                         and EF_HI_FRCN parameter lists.
C  11 Apr 00 @EPA-eg            hrloop revisions
C  21 Mar 00 @EPA-eg            hrloop revisions
C  30 Dec 99 @DynTel-ddj 1-021 Added 1.5 x Tier2 Standard cap on
C            after repair emissions.
C  21 Oct 99 @Dyntel - JWR- Introduced scaling factor REDPROG to 
C            correct the January and July credits for the case when 
C            the program start year equals the calendar year, by
C            removing the spurious Oct -Jan contribution. Fixed OBD credits
C            for 1994 and 1995 model years overwriting all other I&M
C            credits.
C  13 SEP 99 @DynTel-BG Added OBD_IM constant parameter; non-I&M OBD credit for
C            1994 and 1995 needed a sign reversal (OBD for these years increased
C            emissions rather than decreasing them).
C  09 SEP 99 @DynTel-BG Renamed functions HIEMIT.FOR and HN94B.FOR
C            to EF_HI_FRCN and EF_POST93, respectively.
C  11 Aug 99 @DynTel-bag 2-663 Variable names changed, logical and arithmetic
C            expressions now evaluated with same variable types, replaced constants
C            with descriptive parameters, CASE statement for separation of COIM and
C            COIM and RSD non-periodic I&M programs from the other periodic programs,
C            re-introduced max. and min. values for post 1993 after repair rates, and
C            AR94B, determined in IMPROG
C            from IMPROG_NP
C   8 Jul 99 @DynTel-ddj 2-705 Removed bypass around I/M calculation
C            when FTP output option selected.
C   1 Jun 99 @DynTel-ddj 2-700 Changed implementation logic of FTP_FLAG
C            to be consistent with other flags.
C  21 May 99 @DynTel-ddj 2-700, Added code for FTP output option
C            (FTP_FLAG).
C  29 Apr 99 @DynTel-ddj 2-694, Added local variable AR94B.
C  16 Mar 99 @Dyntel-RJD 2-694 Modified IMPROG so that it also returns
C            the I/M additive correction factor for 1994+ model year
C            vehicles equipped with OBD.
C  13 Oct 98 @Dyntel-RJD 2-663 Created IMTEST to correct the 1981-1993
C            light duty gasoline vehicles and trucks basic emission
C            factor for the benefits of IM programs.
C
C  Input on call or from calls:
C
C    argument list: MY,JDX,IP,IV,VMTAGE,ISR,IFDS,ISTD,INERR
C    common blocks:
C    /EVAL/   MEVAL
C    /IMPAR1/ ICYIM
C    /IMPAR2/ ITEST
C    /IMPAR5/ DISCNT
C    /LDGIMC/ WAVRDC
C    /LDGOBD/ OBDIM,WMIL
C
C  Output on return:
C
C    function: IMPROG
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name        Type              Description
C  ------       ----  -----------------------------------------------------
C  IMPGM         I    Local I/M program number.
C  HIGH          I    Integer whose value is 2, indicating high emitter class.
C  HIGHAVE       R    Local high emission rate.
C  IMCRD         R    Fractional reduction in basic emissions from I/M (always negative).
C  IMPGM         I    I/M program number. The user may enter more than one
C                     program for each run. Each program is specified for
C                     a range of vehicle types, calender and model years.
C
      IMPLICIT NONE
      INCLUDE 'EVAL.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'LDGIMC.I'
      INCLUDE 'LDGOBD.I'
C
C     Declare external functions.
C
      INTEGER FINDIM
      REAL    AFTRPREF
      REAL    IDRATE
C
C     Declare parameters.
C
      INTEGER, INTENT(IN) :: ICY
      INTEGER, INTENT(IN) :: ACTUAL_MY
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: JDX
      INTEGER, INTENT(IN) :: IH
      INTEGER, INTENT(IN) :: IROAD
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: ISR
      INTEGER, INTENT(IN) :: IFDS
      INTEGER, INTENT(IN) :: ISTD
      INTEGER, INTENT(IN) :: IDX
      REAL,    INTENT(IN) :: HIGHEF
      INTEGER, INTENT(INOUT) :: INERR
      REAL,    INTENT(IN) :: CAT
      REAL,    INTENT(IN) :: VMTAGE
C
C     Declare local variables.
C
      REAL    BIENADJ
      REAL    FRAC_O
      REAL    FRAC_N
      REAL    IDSUM
      REAL    IMCREDIT
      REAL    OBSUM
      REAL    RPREF
      REAL    RPREF_NET
      REAL    IDR
      REAL    PRCNT_RPRD
      REAL    PRCNT_NTIDD
      REAL    PRCNT_NCOMP
      REAL    PRCNT_WVRS
      INTEGER IMPGM
      INTEGER, SAVE :: IEXHAUST = 1
      INTEGER, SAVE :: JAN=1, JUL = 7
      INTEGER, SAVE :: ANNUAL=1, BIENNIAL=2
      INTEGER, SAVE :: ISTD_TIER0=1
C
      PRCNT_RPRD = 0.0
      PRCNT_NTIDD= 0.0
      PRCNT_NCOMP= 0.0
      PRCNT_WVRS = 0.0
      IDSUM      = 1.0
      IMPGM      = 0
      FRAC_O     = 0.0
      FRAC_N     = 0.0
      IMCREDIT   = 0.0
      OBSUM      = 0.0
      RPREF      = 0.0
      RPREF_NET  = 0.0
      IMPROG     = 0.0
C
C         Test type indices greater than or equal to 12 correspond to
C         OBD I&M programs (see also GETIPR.FOR) :
C
C          1 = idle test
C          2 = 2500 rpm idle test
C          3 = loaded idle test
C          4 = IM240 test
C          5 = ASM test procedure
C          6 = ASM 2525 PHASE-IN
C          7 = ASM 2525 FINAL
C          8 = ASM 5015 PHASE-IN
C          9 = ASM 5015 FINAL
C         10 = ASM 2525/5015 PHASE-IN
C         11 = ASM 2525/5015 FINAL
C         12 = OBD I/M
C
C
C Calculate the correction factor corresponding to the program
C number, IMPGM.
C
C Search algorithm to find I/M program numbers for the emissions
C evaluation date by model year and vehicle type
C
      IMPGM = FINDIM(ICY,MY,IV,IEXHAUST)
C
      IF ( IMPGM==0     .OR.
     *     MY<1981      .OR.
     *     JDX==1       .OR.
     *     (ICYIM(1,IMPGM)==ICY.AND.MEVAL==JAN)
     *  ) GOTO 99
C
C  Determine the emission rate of High emitters after successful repairs
C
      RPREF = AFTRPREF(MY,ACTUAL_MY,IDX,ICY,JDX,IP,IV,ISR,IFDS,
     *                           ISTD,IH,IROAD,CAT,VMTAGE,IMPGM,INERR)
C
      CALL CLCIMPAR(MY,JDX,IP,IV,ISR,IFDS,ISTD,IMPGM,VMTAGE,PRCNT_RPRD,
     *              PRCNT_NTIDD,PRCNT_NCOMP,PRCNT_WVRS,IDR,INERR)
C
C  Determine the base after repair average emissions of highs.
C
C   EF_NTIDD=EF_HIGH*(1.0-IDR)*(1.0-NONCOM)                 EMISSIONS FROM UNIDENTIFIED HIGH EMISSION VEHICLES
C   EF_WVRS =EF_HIGH*WAVRDC(IMPGM)*IDR*WAIVER*(1.-NONCOM)   EMISSIONS FROM WAIVERED HIGH EMISSION VEHICLES
C   EF_NCOMP=EF_HIGH*NONCOM                                 EMISSIONS FROM NONCOMPLAINT HIGH EMISSION VEHICLES
C   EF_RPRD =EF_RPRD*IDR*(1.0-WAIVER)*(1.0-NONCOM)          EMISSIONS FROM REPAIRED HIGH EMISSION VEHICLES
C   EF_AFTR_RPR=EF_RPRD + EF_NTIDD + EF_NCOMP + EF_WVRS     NET EMISSIONS FROM REPAIRED AND
C                                                           NON-REPAIRED HIGH EMISSION VEHICLES
C
C     NET EMISSIONS FROM REPAIRED HIGH EMISSION VEHICLES:
C
      RPREF_NET = HIGHEF*PRCNT_NTIDD  !EMISSIONS FROM UNIDENTIFIED HIGH EMISSION VEHICLES
     *          + RPREF *PRCNT_RPRD   !EMISSIONS FROM REPAIRED HIGH EMISSION VEHICLES
     *          + HIGHEF*PRCNT_NCOMP  !EMISSIONS FROM NONCOMPLAINT HIGH EMISSION VEHICLES
     *          + HIGHEF*WAVRDC(IMPGM)*PRCNT_WVRS   !EMISSIONS FROM WAIVERED HIGH EMISSION VEHICLES
C
C     adjust benefit if the program is an biennial program
      BIENADJ=1.0
      IF(IFREQ(IMPGM)==BIENNIAL) BIENADJ=BIENRATIO(JDX,IP)
      IMCREDIT = MAX(0.0,(HIGHEF-RPREF_NET)*BIENADJ)
C
C     BIENNIAL/ANNUAL ADJUSTMENTS WHICH DEPEND UPON VEHICLE AGE
C
      FRAC_O=OLDFRAC
      FRAC_N=NEWFRAC
C
      IF(JDX==2.AND.ICYIM(1,IMPGM)<ICY) FRAC_N=0.0
      IF(JDX>=2.AND.ICYIM(1,IMPGM)==ICY.AND.MEVAL==JUL) THEN
        FRAC_N=0.0
        FRAC_O=(2.0/3.0)*OLDFRAC
      ENDIF
C
      IF(IFREQ(IMPGM)==ANNUAL) THEN
        IF(JDX==2) IMCREDIT = FRAC_O*IMCREDIT
C
      ELSEIF(IFREQ(IMPGM)>ANNUAL) THEN
        IF(JDX==2) IMCREDIT = (FRAC_O/2.0)*IMCREDIT
        IF(JDX==3) IMCREDIT = (FRAC_O+(FRAC_N/2.0))*IMCREDIT
C
      ENDIF
C
C     Allows for 1996+ Exhaust I/M tests.  It ratios the Exhaust test
C     Identification rate to the OBD I/M Identification rate.
C
      IDSUM = 1.0
      IF(IMPGM>0) THEN
        IF(MY>=1996 .AND. ITEST(IMPGM)<12) THEN !Exhaust I/M After 1996 MY.
C
          IF(WMIL==0.0 .OR. OBDIM==0.0) THEN
            IDSUM = 1.0
          ELSE
C           Use istd=1 (tier0 standards)
            IDSUM=IDRATE(MY,JDX,IP,IV,ISR,IFDS,ISTD_TIER0,IMPGM,VMTAGE,
     *                   INERR) !Exhaust I/M ID Rate
            OBSUM= WMIL * OBDIM !OBD I/M ID Rate
            IDSUM = IDSUM / OBSUM
            IF(IDSUM>1.0) IDSUM = 1.0
            IF(IDSUM<0.0) IDSUM = 0.0
          ENDIF
        ENDIF
      ENDIF
C
C  Apply adjustments to the basic I/M benefit.
C
      IMPROG = IMCREDIT*IDSUM*DISCNT(IP,INTYP(IMPGM))
C
      IMPROG=MAX(0.0,IMPROG)
C
 99   RETURN
      END
