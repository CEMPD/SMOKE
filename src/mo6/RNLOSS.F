      SUBROUTINE RNLOSS(IH,IVTL,MY,ACTUAL_MY,JDX,SPD,RNLOSSTEM,INERR)
C
C  RNLOSS computes the an hourly running loss in grams/mile. There is
C  a running loss value for each model year and gas vehicle class at the
C  given RVP and running Loss temperature using 4-point linear interpolation.
C  Note: This routine is not called for any Diesel vehicles or for motorcycles.
C
C  Called by RNLOSSSUM.
C
C  Calls ETPCAL, IRNLPT.
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  18 Jan 01 @EPA-bag  2007 rule for heavy duty vehicles
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch
C 20 Nov 00 @EPA-bag bug#303 Added weekly index to check on user supplied trip length
C           fraction values, i.e., the TLVMTU's
C 12 Sep 00 @EPA-elg ETP calculations corrected for running losses (Bug #249).
C 20 Jul 00 @EPA-bag running loss; New evaporative weights and moved
C           resting loss subtraction
C 18 Jan 00 @DynTel-MLA 1-011  Added ICY to the parameter list.
C 17 May 2000 AIR Task 03: Rmoved hourly dependency of RVPHS.
C 19 Jan 00 @Dyntel-JWR 1-017 Made RSTLSS_EF the hourly rest loss total
C 15 Nov 1999 @DynTel-BG 1-006 Moved RNLOSSTEM out of common to RNLOSS parameter list;
C             reassigned resting loss variable
C  5 May 1999 @DynTel-HQ 2-664 Modify this subroutine to calculate the
C                              running loss for specific speed.
C 13 Jan 99 - DynTel@RAG 2-664  Expanded TLVMT array to include hourly and
C             weekday/end indices. Running losses now have a hourly distribution.
C             Included temperature dependence to running loss parameters.
C
C  Input on call:
C
C    argument list: MY,IVTL
C    common blocks:
C    /CITRV2/ RVPHS
C    /RESTLS1/ RSTLSS_EF
C    /RUNLS1/ TLVMT,TLVMTU,MAXTL
C    /RUNLS2/ TLEMI,TLFAIL
C    /RUNLS3/ RNLRVP,RNLTEM,RNLSPD,MAXRNL
C    /TEMPS/  TEMHRLY
C    /VMXCOM/ VCOUNT
C    /VVNAMS/ VVGAS  (IVTYPE.I)
C    /VVTYPS/ VTGASBUS  (IVTYPE.I)
C    /WEEKLY/ IWEEK
C
C  Output on return:
C
C    argument list: RNLOSSTEM
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  EM11        R    running loss emission (or offset) at lower boundary
C                   temperature and lower boundary RVP
C  EM12        R    running loss emission (or offset) at lower boundary
C                   temperature and upper boundary RVP
C  EM21        R    running loss emission (or offset) at upper boundary
C                   temperature and lower boundary RVP
C  EM22        R    running loss emission (or offset) at upper boundary
C                   temperature and upper boundary RVP
C  EMT1        R    interpolated running loss emission (or offset)
C                   at the input RVP and the lower boundary temperature
C  EMT2        R    interpolated running loss emission (or offset)
C                   at the input RVP and the upper boundary temperature
C  ETPF        R    holds ETP factors for purge/pressure failures
C  ETPP        R    holds ETP factors for purge/pressure passes
C  FACTOR      R    misc. intermediate factor for interpolations and such
C  RSTLS       R    resting loss rate to be subtracted off
C  RUFAIL      R    running loss emission offset for canister pressure/
C                   purge fail vehicles
C  RUPASS      R    running loss emission rate for pass vehicles
C  RVP1        R    lower boundary RVP
C  RVP2        R    upper boundary RVP
C  RVPINT      R    running loss RVP used for interpolation/extrapolation
C  RVPRL       R    running loss RVP
C  TEMRL       R    running loss temperature
C  TMP1        R    lower boundary temperature
C  TMP2        R    upper boundary temperature
C  IH          I    hourly index passed from EVP_HRLY
C  IHTL        I    running loss hourly index IHTL = IH for IH<15
C                                             IHTL = 14 for IH>13
C  ITL         I    6 TRIP LENGTHS
C  ISPD        I    3 SPEEDS
C  ITLR        I    4 RVP'S
C  ITLT        I    4 TEMP.'S
C  IVTL        I    14 VEHICLE CLASSES; IV FOR GAS VEHICLES AND GAS BUS
C  IVGAS       I    gasoline fueled vehicle type index
C  IVRNL       I    running loss vehicle index
C  MY          I    Vehicle model year; MY=ICY+IDX-25
C  IGTL        I    4 MODEL YEARS
C  IG          I    model year pointer (Pre-Controlled groups use TLFAIL only)
C                   IG = IRNLPT(MY,IVRNL,INERR)
C  SPDRL       R    running loss speed
C  SPD         R    speed passed from RNLOSSSUM
C  GLL_TOTAL   R    running loss liquid leak emission units of grams/hour
C  RUGLL       R    running loss liquid leak emission after subtraction of
C                   resting loss liquid emission in units of grams/mile
C  GLL_REST    R    resting loss liquid leak emission units of grams/hour
C  ENHWTF      R    post 1995 evaporative weights for FAIL either the
C                   purge or the pressure tests
C  ENHWTP      R    post 1995 evaporative weights for PASSED both purge
C                   and pressure tests
C  PREWTF      R    pre-1996 evaporative weights for FAIL either the
C                   purge or the pressure tests
C  PREWTP      R    pre-1996 evaporative weights for PASSED both purge
C                   and pressure tests
C  RUMY        I    running loss model year; has a maximum of 1995
C  REAL_MAXTL  R    MAXTL, maximum # of travel intervals, converted to a
C                   REAL, i.e., REAL(MAXTL)
C  RST         R    resting loss less the resting loss liquid leaks used
C                   to subtract from the running loss emission
C
C  Notes:
C
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1.
C  RNLOSS was modified in Version 4.1 to activate the use of running loss
C  values by trip length rather than using trip-length weighted values
C  and to interpolate running loss emission by vehicle speed.
C  The tampering effects calculated in CAPOFF and CANOFF were replaced
C  with pressure/purge failed vehicle information in MOBILE4.1
C
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6. Added IVRNL to map
C            IVTL.
C  16 December 1997 Dyntel@MMS modified multispeed to a single speed.
C                   SPD changed from an array to a single varaible.
C
      USE DIRNLDAT, ONLY : W95,W99
      IMPLICIT NONE
      INCLUDE 'CITRV2.I'
      INCLUDE 'FLAGS5.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'RESTLS1.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RUNLS1.I'
      INCLUDE 'RUNLS2.I'
      INCLUDE 'RUNLS3.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'WEEKLY.I'
C
C     Declare external functions
C
      INTEGER IRNLPT
      REAL    EVWTENH
      REAL    EVWTPRE
      REAL    GLL_WGT
      REAL    ETPCAL
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: IH
      INTEGER, INTENT(IN) :: IVTL
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: ACTUAL_MY
      INTEGER, INTENT(IN) :: JDX
      INTEGER, INTENT(INOUT) :: INERR
      REAL, INTENT(IN)    :: SPD
      REAL, INTENT(OUT)   :: RNLOSSTEM
C
C     Declare local variables
C
      INTEGER IG
      INTEGER IGTL
      INTEGER IHTL
      INTEGER ISPD1
      INTEGER ISPD2
      INTEGER ITL
      INTEGER ITLR
      INTEGER ITLT
      INTEGER IVRNL
      REAL    RST
      REAL    EM11
      REAL    EM12
      REAL    EM21
      REAL    EM22
      REAL    EMT1
      REAL    EMT2
      REAL    ETPF
      REAL    ETPP
      REAL    FACTOR
C
C     For now, emissions from Gross Liquid Leakers are hard coded using
C     the values specified in report M6.EVP.009. All emissions are in
C     units of grams/hour. Total emissions are the sum of running and
C     resting loss emissions. Subtracting the resting loss emissions
C     from the total emissions yields the running loss. Dividing by the
C     speed gives emissions in units of grams/mile.
C
      REAL, SAVE ::  GLL_TOTAL = 345.94
      REAL, SAVE ::  GLL_REST =    9.16
      REAL       ::  RUFAIL
      REAL       ::  ENHWTF
      REAL       ::  ENHWTP
      REAL       ::  PREWTF
      REAL       ::  PREWTP
      REAL       ::  RUGLL
      INTEGER    ::  RUMY
      INTEGER, SAVE :: IYES=1
C 
      REAL    REAL_MAXTL
      REAL    RUPASS
      REAL    RVP1
      REAL    RVP2
      REAL    RVPINT
      REAL    RVPRL
      REAL    SPDRL
      REAL    TEMRL
      REAL    TMP1
      REAL    TMP2
C
C     Initializations
C
      IVRNL  = IVTL
      RNLOSSTEM = 0.0
      RUPASS = 0.0
      RUFAIL = 0.0
      RVPRL  = 0.0
      SPDRL  = SPD
      TEMRL  = 0.0
C
      IHTL = IH
C
C     Use the 14th time intervals' TLVMTU's for hourly intervals
C     14 through 24
C
      IF( IH > 14 ) IHTL = 14
C
      IF ( IVTL == VTGASBUS ) IVRNL = 14
C
C     Initialize TLVMTU array to built-in values from TLVMT
C     unless some value was already read in (or RNLOSSHOUR
C     was already called)
      IF(TLVMTU(1,1,IWEEK) < 0.0)
     *        TLVMTU(1:6,1:14,IWEEK)=TLVMT(1:6,1:14,IWEEK)
C
C     Running losses from non-GLL vehicles are zero if any evaporative
C     temperature is 40F or less, or if minimum temperature is 25F
C     or less. But gross liquid leakers have non-zero emissions even
C     at temperatures that are out of range for ordinary running
C     loss emissions.
C
      IF( TEMMAX      <= 40.0 .OR.
     *    TEMHRLY(IH) <= 40.0 .OR.
     *    TEMMIN      <= 25.0 ) GOTO 85
C
C     Choice of RVP depends on whether or not cert fuel change applies to the
C     current MY. Also, set running loss temperature.
C
      RVPRL = RVPHS
      TEMRL = TEMHRLY(IH)
C
C     Reset RVP and Running Loss Temperature if above running loss limits
C
      IF(RVPRL < RNLRVP(1)) RVPRL=RNLRVP(1)
      IF(RVPRL > RNLRVP(4)) RVPRL=RNLRVP(4)
C
C     Maximum temperature is set at 105 degrees Fahrenheit
C
      IF(TEMRL > RNLTEM(4)) TEMRL = RNLTEM(4)
C
C     Set the interpolation RVP.
C     (Allow for extrapolation 0.5 RVP below RNLRVP lower limit.)
C
      RVPINT = RVPHS
      IF(RVPINT < RNLRVP(1)-0.5) RVPINT=RNLRVP(1)-0.5
      IF(RVPINT > RNLRVP(4)) RVPINT=RNLRVP(4)
C
C     Find and assign the bracketing temperature values.
C
      ITLT=4
      IF( TEMRL <= RNLTEM(3) ) ITLT = 3
      IF( TEMRL <= RNLTEM(2) ) ITLT = 2
      IF( TEMRL <= RNLTEM(1) ) ITLT = 1
C
      TMP1=40.0
C
      IF(ITLT > 1) TMP1=RNLTEM(ITLT-1)
      TMP2=RNLTEM(ITLT)
C
C     Find and assign the bracketing RVP values.
C
      ITLR=4
      IF( RVPRL <= RNLRVP(3) ) ITLR = 3
      IF( RVPRL <= RNLRVP(2) ) ITLR = 2
C
      RVP1=RNLRVP(ITLR-1)
      RVP2=RNLRVP(ITLR)
C
C     Find and assign the bracketing SPD values (TLVMT groups)
C
      SPDRL=SPD
      ISPD1=1
      ISPD2=1
      IF( SPDRL > RNLSPD(1) ) ISPD2=2
      IF( SPDRL >= RNLSPD(2) ) ISPD1=2
      IF( SPDRL > RNLSPD(2) ) ISPD2=3
      IF( SPDRL >= RNLSPD(3) ) ISPD1=3
      FACTOR=0.0
      IF( ISPD1 /= ISPD2 ) THEN
        FACTOR=( SPDRL-RNLSPD(ISPD1) )/( RNLSPD(ISPD2)-RNLSPD(ISPD1) )
      END IF
C
C     Determine model year pointer (Pre-Controlled groups use TLFAIL only)
C     Running losses for 96+ vehicles are based on the 1995 MY rate.
C
      RUMY = MY
      IF (MY > 1995) RUMY = 1995
      IG = IRNLPT(RUMY,IVRNL,INERR)
      IF(INERR.GT.0) GOTO 99
C
C     (For pre-control vehicles do the fail pressure/jump to line 80.)
C
      IF( IG == 1 ) GO TO 80
C
C     Reduce IG by one to use correct TLEMI components for pass IV
C
      IGTL = IG - 1
C
C     Assign bracketing running loss emission values
C     Interpolate on vehicle speed as the values are added up
C     Subtract off resting losses. First, convert resting losses
C     into grams/mile. The TLEMI are in units of grams/mile.
C
      EM11 = 0.0
      EM12 = 0.0
      EM21 = 0.0
      EM22 = 0.0
C
C     Vehicles which have PASSED purge and pressure tests.
C     Sum over MAXTL (different travel lengths) and weight each emission
C     by the appropriate fraction of vehicle miles travelled,
C     TLVMTU(ITL,IHTL,IWEEK), i.e., the hourly per trip activity fraction.
C
      DO ITL = 1, MAXTL
C
        IF( ITLT > 1 ) THEN
          EM11 = EM11 +
     *    (
     *      TLEMI(ITL,ISPD1,ITLR-1,ITLT-1,IGTL,IVRNL) +
     *      FACTOR * ( TLEMI(ITL,ISPD2,ITLR-1,ITLT-1,IGTL,IVRNL) -
     *                 TLEMI(ITL,ISPD1,ITLR-1,ITLT-1,IGTL,IVRNL) )
     *    ) * TLVMTU(ITL,IHTL,IWEEK)
C
          EM12 = EM12 +
     *    (
     *      TLEMI(ITL,ISPD1,ITLR,ITLT-1,IGTL,IVRNL) +
     *      FACTOR * ( TLEMI(ITL,ISPD2,ITLR,ITLT-1,IGTL,IVRNL) -
     *                 TLEMI(ITL,ISPD1,ITLR,ITLT-1,IGTL,IVRNL) )
     *    ) * TLVMTU(ITL,IHTL,IWEEK)
        END IF
C
        EM21 = EM21 +
     *    (
     *      TLEMI(ITL,ISPD1,ITLR-1,ITLT,IGTL,IVRNL) +
     *      FACTOR * ( TLEMI(ITL,ISPD2,ITLR-1,ITLT,IGTL,IVRNL) -
     *                 TLEMI(ITL,ISPD1,ITLR-1,ITLT,IGTL,IVRNL) )
     *    ) * TLVMTU(ITL,IHTL,IWEEK)
C
        EM22 = EM22 +
     *    ( TLEMI(ITL,ISPD1,ITLR,ITLT,IGTL,IVRNL) +
     *      FACTOR * ( TLEMI(ITL,ISPD2,ITLR,ITLT,IGTL,IVRNL) -
     *                 TLEMI(ITL,ISPD1,ITLR,ITLT,IGTL,IVRNL) )
     *    ) * TLVMTU(ITL,IHTL,IWEEK)
C
      END DO
C     End loop over time of travel percentages for vehicles which
C     have passed canister pressure/purge tests.
C
C     Interpolate first on RVP, then temperature
C
      EMT1 = RVPINT*(EM12-EM11)/(RVP2-RVP1) +
     *       (RVP2*EM11-RVP1*EM12)/(RVP2-RVP1)
      EMT2 = RVPINT*(EM22-EM21)/(RVP2-RVP1) +
     *       (RVP2*EM21-RVP1*EM22)/(RVP2-RVP1)
C
      RUPASS = TEMRL*(EMT2-EMT1)/(TMP2-TMP1) +
     *         (TMP2*EMT1-TMP1*EMT2)/(TMP2-TMP1)
C
C     Now, compute the running loss for canister pressure/purge FAIL vehicles.
C     This is also the only calculation needed for vehicles without emissions
C     controls, i.e., pre-control vehicles.
C
C     Assign bracketing running loss pressure/purge fail offset values
C     Don't forget to subtract off fraction of resting loss rate (g/hr).
C     Fail numbers, TLFAIL, are in unit of grams. See (BD34.FOR.) (Resting
C     loss emissions, RSTLSS_EF, are in units of grams/hour.)
C
   80 EM11 = 0.0
      EM12 = 0.0
      EM21 = 0.0
      EM22 = 0.0
C
C     Vehicles which have failed canister pressure/purge tests.
C     The TLFAIL values are in grams by bag.
C     Sum over ITL [MAXTL(=6) different travel lengths] and weight each emission
C     by the appropriate fraction of vehicle miles travelled for that particular
C     travel partition, TLVMTU(ITL,IHTL,IWEEK), i.e., the hourly per trip activity
C     fraction. This sum yields hourly resting loss emissions based on the number
C     of miles for each portion of the trip and the time of day of the trip.
C
      REAL_MAXTL=REAL(MAXTL)
      DO ITL = 1, MAXTL
C
C       Figure conversion factor for g-->g/mi plus the %TLVMT
C
        FACTOR = TLVMTU(ITL,IHTL,IWEEK)*REAL_MAXTL / (SPD*FLOAT(ITL))
C
        IF( ITLT > 1 ) THEN
          EM11 = EM11 + TLFAIL(ITL,ITLR-1,ITLT-1)*FACTOR
          EM12 = EM12 + TLFAIL(ITL,ITLR  ,ITLT-1)*FACTOR
        END IF
        EM21 = EM21 + TLFAIL(ITL,ITLR-1,ITLT)*FACTOR
        EM22 = EM22 + TLFAIL(ITL,ITLR  ,ITLT)*FACTOR
C
      END DO
C     End loop over time of travel percentages for vehicles which
C     have failed canister pressure/purge tests.
C
C     Interpolate first on RVP, then temperature
C
      EMT1 = RVPINT * (EM12-EM11)/(RVP2-RVP1) +
     *      (RVP2*EM11-RVP1*EM12)/(RVP2-RVP1)
      EMT2 = RVPINT*(EM22-EM21)/(RVP2-RVP1) +
     *       (RVP2*EM21-RVP1*EM22)/(RVP2-RVP1)
C
      RUFAIL = TEMRL*(EMT2-EMT1)/(TMP2-TMP1) +
     *         (TMP2*EMT1-TMP1*EMT2)/(TMP2-TMP1)
C
      IF(IG == 1) RUPASS=RUFAIL
C
C     Reset fail if less than pass.
C
      IF( RUFAIL < RUPASS ) RUFAIL = RUPASS
C
C     Compute the rate for GLL, using the hard-coded data points, and
C     divide by speed to convert to units of grams/mile.
C
   85 RUGLL = (GLL_TOTAL - GLL_REST)/SPD
C
C     Initialize the weighting factors for the different pass/fail
C     strata. Phase in the Evap Test Procedure over the model year
C     range 1996-1998 by taking an average of the 1995 and the 1999
C     emision factors. The weighting factors that are used for the 
C     averages are 80/20 for 1996, 60/40 for 1997, and 90/10 for 1998.
C     1999 pass emisisons are defined to be 20% of the 1995 pass
C     emissions and 1999 fail emissions are defined to be 70% of the
C     1995 fail emission levels.
C
C     Years 1999 and later use only the 'enhanced' weights. Years 1995
C     and earlier use only the 'pre-enhanced' weights. Years 1996 to
C     1998 use blended rates.
C
      IF (MY < 1999) THEN
        PREWTP = EVWTPRE(MY,ACTUAL_MY,JDX,1,3)
        PREWTF = EVWTPRE(MY,ACTUAL_MY,JDX,2,3) +
     &           EVWTPRE(MY,ACTUAL_MY,JDX,3,3) +
     &           EVWTPRE(MY,ACTUAL_MY,JDX,4,3)
      ELSE
        PREWTP = 0.0
        PREWTF = 0.0
      ENDIF
C
      IF (MY > 1995) THEN
        ENHWTP = EVWTENH(MY,ACTUAL_MY,JDX,1,3)
        ENHWTF = EVWTENH(MY,ACTUAL_MY,JDX,2,3) +
     &           EVWTENH(MY,ACTUAL_MY,JDX,3,3) +
     &           EVWTENH(MY,ACTUAL_MY,JDX,4,3)
      ELSE
        ENHWTP = 0.0
        ENHWTF = 0.0
      ENDIF
C
C     ETPDUM implements the Enhanced Test Procedure effects for Running losses.
C
      ETPP = 1.0 -  ETPCAL(6,MY,1,RVPRL,TEMRL)
      ETPF = 1.0 - (ETPCAL(6,MY,2,RVPRL,TEMRL) +
     *              ETPCAL(6,MY,3,RVPRL,TEMRL))/2.0
C
      IF (MY <= 1995) THEN
        RNLOSSTEM = RUPASS*PREWTP+RUFAIL*PREWTF+
     *              RUGLL*GLL_WGT(MY,ACTUAL_MY,3,JDX)
      ELSE
        RNLOSSTEM = ETPP*RUPASS*ENHWTP +
     &              ETPF*RUFAIL*ENHWTF +
     &              RUGLL*GLL_WGT(MY,ACTUAL_MY,3,JDX)
      ENDIF
C
      IF(IVTL>=6.AND.MY>2007.AND.HDR2007==IYES
     &   .AND.TIER2_FLAG==IYES) THEN
        RST=(RSTLOS_FR_RNLOSS-GLL_WGT(MY,ACTUAL_MY,2,JDX)*GLL_REST)/SPD
      ELSE
        RST=(RSTLSS_EF-GLL_WGT(MY,ACTUAL_MY,2,JDX)*GLL_REST)/SPD
      ENDIF
      RST=AMAX1(0.0,RST)
      RNLOSSTEM=(RNLOSSTEM-RST)
      IF(RNLOSSTEM < 0.) RNLOSSTEM=0.0
C
C     If speed is below 7 or above 48, extrapolate the running loss rate
C
      IF(SPDRL < RNLSPD(1))
     *           RNLOSSTEM=RNLOSSTEM*RNLSPD(1)/SPDRL
      IF(SPDRL > RNLSPD(3))
     *           RNLOSSTEM=RNLOSSTEM*RNLSPD(3)/SPDRL
C
  99  RETURN
      END
