      REAL FUNCTION HOTSOK(IH,MY,ACTUAL_MY,IVGAS,JDX)
C
C  HOTSOK selects an average hourly value for hot soak HC emissions
C  (given an hourly temperature and given a model year, MY. It
C  then uses this average hourly number to compute hourly hot soak
C  emissions considering soak time distributions. The final result,
C   is in units of grams/hour.
C
C  Called by EVP_HRLY.
C
C  Calls EVMAIN, ETPCAL.
C
C  Changes: (Last change first)
C
C  02 May 02 @EPA-bag Added code for user input ofhot soak liquid leak rates
C  27 Feb 01 @EPA-BAg bug314 revised tier2 effects on heavy duty vehicles
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch
C  24 Aug 00 @EPA-djb Bug 233, Applied Arcadis equations to fuel
C            injected trucks as well as carbureted.  Added IFDS to the
C            call to EVMAIN. Allowed KEYEQ >1 to loop by IFDS.
C  01 Aug 00 @EPA-elg Resting losses are now added as an average rather than by PP BIN.
C  21 Jul 00 @EPA-elg Added code to account for Tier2 Evap on IV=6 vehicles (2b's).
C  11 aUG 00 @EPA-bag Task 06: rearranged if_then conditional for user
C            activity fractions
C  17 Jul 00 AIR Task 06: Added code to apply user supplied
C            Hot Soak Activity Fractions.
C  26 Jul 00 @EPA-bag running loss; removed citcin.i include previously held
C                 purge/pressure pass fail fractions
C  13 Jul 00 @EPA-ddj Bug191 Modified equation for HSPF when TEMHRLY
C            less than 70F.
C  17 May 00 AIR Task 03: Removed hourly dependency of RVPHS.
C   2 Jun 00 @EPA-ddj Bug 86, Changed FDS conditional from VVLDGAS to
C            VVGAS.
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations
C  13 Apr 00 @Dyntel- BG 1-014  Use LEAKYDAT variable, ILIQUID_LEAKS; removed
C            GLL_SWITCH
C  22 Mar 00 @Dyntel- JWR 1-014 Updated to include gross liquid leakers
C            and new evaporative weights as well as Tier2 phase-in. 
C  14 Jan 00 @Dyntel-JWR 1-017 Made RSTLSS_EF the hourly rest loss total
C  29 Nov 99 EPA@ELG Added code to handle ETP/OBD Phase-in for MY 1996-98.
C  30 Sep 99 EPA@ELG Added trip per day distribution to calculate hot soaks
C            in units of grams per hour.
C  29 Sep 99 EPA@ELG Added new temperature correction logic and values for 1981+ MYs
C  04 Nov 99 @DynTel-BAG 1-006 Deleted the HSK(IH) array and HOTSOK is now a
C            function
C  03 Dec 98 @DynTel-BAG 2-664 Replaced TEMEVP(IHRD) with TEMHRLY(IH). Deleted
C            TEMFLG. Added hurly hot soak code from EVAPHOUR.
C  06 Nov 98 @ DynTel-ZK 2-000 explicit Typing 
C
C  Input on call:
C
C    argument list: IVGAS,JDX,MY
C    modules :
C     LEAKYDAT - ILIQUID_LEAKS
C    common blocks:
C    /CITRV2/ RVPHS
C    /ETPPAR/ IETPYR,ETPPIV,NETPYR                  ! THIS IS NEW
C    /EVAHS1/ HSCOEF,HSEQ, HS_HRLY_TRP_FRCN
C    /EVAPACT/ HSACTUSER, HSACTFLAG
C    /EVAPAR/ ISTDC,ISTDT,ISTDP,FINJ
C    /EVAPHD/ IVVMAP
C    /EVAPHR/ IH
C    /EVAHS1/ HSIV5
C    /FLAGS4/ HCFLAG
C    /HOTSOK1/ HSK
C    /MAXIMA/ MAXYRS
C    /RESTLS1/ RSTLSS_EF
C    /TAMOUT/ THS
C    /TEMPS/  TEMHRLY,TEMMIN,TEMMAX
C    /VVNAMS/ VVGAS, VVHDGV, VVLDGAS, VVLDGT12, VVLDGT34  (IVTYPE.I)
C    /VVTYPS/ VTGASBUS, VTLDGV, VTMC  (IVTYPE.I)
C    /WEEKLY/ IWEEK
C
C  Output on return:
C
C  Function value HOTSOK
C
C  Local array subscripts
C
C  HSTCF  (3) - HSTCF  (IFDS)
C
C  Local variable / array dictionary
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IH       I    hourly index
C  JDX      I    vehicle age index; JDX=25-IDX+1
C  MY       I    Vehicle model year; MY=ICY-JDX+1=ICY+IDX-25
C  IVGAS    I    gasoline fueled vehicle type index
C  IVTEMP   I    hot soak vehicle type index
C  A        R    Variable used with hot soak temperature correction factor
C                It has values of 1 or HSTCF depending on vehicle type and
C                and carburetor type
C  HSGLTOT  R    Liquid leaker total (summed over fuel systems)
C  HSPASS   R    Pass emission rate to insure FPressure and FPurge
C                never fall below the pass rate
C  HSPF     R    the EVMAIN hot soak loss for Pass, Failed Purge, Failed
C                Pressure and controlled and uncontrolled
C  HSTCF    R    hot soak temperature correction factor for CARB, TBI,
C                and PBI
C  HSTEMP   R    temperature used in the hot soak temperature correction
C                factor equation
C  IVV      I    vehicle (IVV=1) or truck (IVV=2) classification switch
C  RIVAL    R    hot soak RVP adjusted for weathering
C  FULADJ   R    FTP vs In-Use fuel level adjustment factor for fuel
C                injected vehicles
C  HS_RESID R    holds the non Tier2 hotsoak value for IV=6 vehicles.
C  MAXKEQ   I    number of purge/pressure test to loop through
C  MAXFDS   I    number of fuel delivery systems to loop through
C  KEYEQ    I    purge/pressure test index
C  IHST     I    hot soak hourly index
C  IMIN     I    hot soak soak time loop index
C  IV6T2    R    phase-in of IV=6 (MDGV) for Tier2.  First value is for
C                CY2007, then CY2008, and CY2009+
C  IV6HS    R    Hot Soak Correction Factor for IV=6 (MDGV).
C  CUM      R    cumulative hot soak fraction
C  PREV     R    hot soak fraction at earlier soak time
C  DIF      R    difference of hot soak fractions, CUM-PREV
C  HMIN     R    hot soak emission at soak time IMIN
C  HSBASE   R    base MOBILE5 hot soak emission used for hourly hot soak
C                determination from soak time distributions
C  TRIPFR   R    fraction of vehicle fleet with daily trips
C
C  Notes:
C
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed in MOBILE4.1.
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C  Sep-21-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year
C  removed IY from the argument list, and the calculation of ICY is no longer
C  needed and thus removed IY and ICY from the calls to EVMAIN
C
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS, HSUSERGL_INPT, HSGLL_USER
C
      IMPLICIT NONE
C
      INCLUDE 'CITRV2.I'
      INCLUDE 'ETPPAR.I'
      INCLUDE 'EVAPACT.I'
      INCLUDE 'EVAHS1.I'
      INCLUDE 'EVAPAR.I'
      INCLUDE 'EVAPHD.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'HOTSK1.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'RESTLS1.I'
      INCLUDE 'STDIST.I'
      INCLUDE 'TAMOUT.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'WEEKLY.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'FLAGS5.I'
C
C     Declare external functions
C
      REAL TR2FCTR
      REAL ETPCAL
      REAL EVWTPRE
      REAL EVWTENH
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: IH
      INTEGER, INTENT(IN) :: IVGAS
      INTEGER, INTENT(IN) :: JDX
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: ACTUAL_MY
C
C
C     Declare local variables
C
      INTEGER  IVTEMP
      INTEGER  ISTUSE
      INTEGER  MAXKEQ
      INTEGER  MAXFDS
      INTEGER  KEYEQ
      INTEGER  IHST
      INTEGER  IVV
      INTEGER  IMIN
      INTEGER  IETPY
      INTEGER, SAVE :: IYES = 1
C
      REAL  CUM
      REAL  DIF
      REAL  HMIN(60)
      REAL  HSBASE
      REAL  HSPASS
      REAL  HSPF
      REAL  HSGLTOT
      REAL  HSTCF(3)
      REAL  HSTEMP
      REAL  PREV
      REAL  PFRTPRE(3)
      REAL  PFRTENH(3)
      REAL  GLWT
      REAL  RIVAL
      REAL  ETPIV
      REAL, SAVE :: FULADJ = 0.88
      REAL, DIMENSION(3), SAVE :: HSGLL=(/16.95,45.0,57.14/)
C
C
C     Initializations
C
      HSBASE=0.0
      HOTSOK=0.0
      PREV = 0.0
C
C     Hot soak is zero if any evaporative temperature is 40F or less,
C     or if minimum temperature is 25F or less.
C
      IF( (TEMMAX.LE.40.0)       .OR.
     *    (TEMHRLY(IH).LE.40.0) )     RETURN
C
C     Initialize the pressure\purge weighting factors
C
      PFRTPRE = 0.0
      PFRTENH = 0.0
C
C     Get the pressure\purge weighting factors
C
      IF(MY.GE.1996) THEN
        PFRTENH(1)=EVWTENH(MY,ACTUAL_MY,JDX,1,1)
        PFRTENH(2)=EVWTENH(MY,ACTUAL_MY,JDX,2,1)
        PFRTENH(3)=EVWTENH(MY,ACTUAL_MY,JDX,3,1) +
     *             EVWTENH(MY,ACTUAL_MY,JDX,4,1)
        GLWT = EVWTENH(MY,ACTUAL_MY,JDX,5,1)
      END IF
C
      IF(MY.LE.1998) THEN
        PFRTPRE(1)=EVWTPRE(MY,ACTUAL_MY,JDX,1,1)
        PFRTPRE(2)=EVWTPRE(MY,ACTUAL_MY,JDX,2,1)
        PFRTPRE(3)=EVWTPRE(MY,ACTUAL_MY,JDX,3,1) +
     *             EVWTPRE(MY,ACTUAL_MY,JDX,4,1)
        GLWT = EVWTPRE(MY,ACTUAL_MY,JDX,5,1)
      END IF
C
      IVTEMP = IVGAS
      IF (IVGAS.EQ.14) THEN
        IVTEMP = VTMC
      ELSEIF (IVGAS.EQ.15) THEN
        IVTEMP = VTGASBUS
      ENDIF
C
      RIVAL=RVPHS
C
C     Lower bound on Hot Soak RVP is 5.0, reset if necessary.
C
      IF ( RIVAL .LT. 5.0 ) RIVAL=5.0
C
C     Select the base temperature to be used for the hot soak correction
C     factor. Use the calculated hot soak temperature as default. However,
C     if that temperature is below 70F, reset to 70F, then interpolate later.
C
      HSTEMP = TEMHRLY(IH)
      IF ( HSTEMP .LT. 70.0 ) HSTEMP = 70.0
C
C     Compute the Hot Soak Temperature Correction Factor for each fuel
C     delivery system.
C
      IVV=IVVMAP(IVGAS)
C
      DO IFDS=1,3
C
C     All 1981 and newer vehicles get the new MOBILE6 correction.
C     All others keep the MOBILE5 correction.
C
           IF(MY.GE.1981 .AND. MY.LE.1985) THEN
             HSTCF(IFDS)=(HSEQ_81(IFDS,IVV,6) +
     *                    HSEQ_81(IFDS,IVV,7)*HSTEMP +
     *                    HSEQ_81(IFDS,IVV,8)*HSTEMP**2) /
     *                    HSEQ_81(IFDS,IVV,9)
           ELSEIF(MY.GE.1986) THEN
             HSTCF(IFDS)=(HSEQ_86(IFDS,IVV,6) +
     *                    HSEQ_86(IFDS,IVV,7)*HSTEMP +
     *                    HSEQ_86(IFDS,IVV,8)*HSTEMP**2) /
     *                    HSEQ_86(IFDS,IVV,9)
           ELSE
             HSTCF(IFDS)=(HSEQ(IFDS,IVV,6) +
     *                    HSEQ(IFDS,IVV,7)*HSTEMP +
     *                    HSEQ(IFDS,IVV,8)*HSTEMP**2) /
     *                    HSEQ(IFDS,IVV,9)
         ENDIF
C
      END DO      !IFDS Loop
C
C     All HDGV AND pre-1981 LDGV/T are 100% carbureted, but are affected by
C     pass/fail.  For these vehicle classes, figure a carbureted number
C     accordingly and then skip the fuel injected rate case.
C     Motorcycles are 100% carbureted and not affected by pass/fail.
C
C     Loop and combine for Passed, Failed Purge and Failed Pressure
C     and for carbureted and fuel injected vehicles.
C     (MC do not have purge or pressure checks.)
C
      HSPASS=0.0
C
      MAXKEQ=3
C
C     Motocycles do not have pressure or purge failures.
C
      IF(IVTEMP.EQ.VTMC) MAXKEQ=1
C
C     KEYEQ is the pressure/purge index
C     KEYEQ = 1 : passing both purge and pressure
C     KEYEQ = 2 : failed purge uncontrolled emissions
C     KEYEQ = 3 : failed pressure uncontrolled emissions
C
      DO KEYEQ=1,MAXKEQ
        MAXFDS=3
C
C     Motorcycles and pre-1981 model years do not have fuel injection estimates.
C
        IF(IVGAS.EQ.14.OR.MY.LT.1981) MAXFDS=1
C
C       Loop through fuel delivery systems
C
C       IFDS = 1 : Carbureted
C              2 : Throttle-body fuel injection
C              3 : Ported or Multi-point fuel injection
C
        DO IFDS=1,MAXFDS
          IF(IFDS.EQ.1) ISTUSE=ISTDC
          IF(IFDS.EQ.2) ISTUSE=ISTDT
          IF(IFDS.EQ.3) ISTUSE=ISTDP
C
          CALL EVMAIN(IH,MY,KEYEQ,IVGAS,ISTUSE,RIVAL,HSPF,IFDS)
C
C         Skip the TCF for Failed Purge and Failed Pressure and only if
C         the last standard group
C
          IF( KEYEQ.GT.1 ) THEN
            IF(IVTEMP.EQ.VTLDGV.AND.ISTUSE.GE.5)  GOTO 60
            IF(VVLDGT12(IVTEMP).EQ.1.AND.ISTUSE.GE.5) GOTO 60
            IF(VVLDGT34(IVTEMP).EQ.1.AND.ISTUSE.GE.3) GOTO 60
            IF(VVHDGV(IVTEMP).EQ.1.AND.ISTUSE.GE.2) GOTO 60
          ENDIF
          HSPF = HSPF * HSTCF(IFDS)
C
C         Save the composite pass rate for future purge/pressure checks.
C
 60       IF(KEYEQ.EQ.1) HSPASS=HSPF
C
C         Don't allow purge or pressure rates to fall below pass rates.
C
          IF(KEYEQ.GT.1.AND.HSPF.LT.HSPASS) HSPF=HSPASS
C
C         If temperature is below 70F, interpolate with assumption
C         that HS=0 AT 40F.
C
          IF(TEMHRLY(IH).LT.70.0)
     *      HSPF=HSPF*(TEMHRLY(IH)-40.0)/30.0
C
C         Evap. Test Procedure adjustment calculation.
C
          ETPIV = 0.0
          IF(MY.GE.IETPYR(1)) THEN
            IETPY = MY - (IETPYR(1) - 1)
            IF(IETPY.GE.4) IETPY = NETPYR
            ETPIV = ETPPIV(IETPY)/100.
          ENDIF
C
C
C       Applied ETP to gas vehicles except Motorcycles.
C
          IF(IVTEMP.EQ.VTMC) THEN
            HSPF = HSPF * PFRTPRE(KEYEQ)*(1-ETPIV) +
     *             HSPF * PFRTENH(KEYEQ)*ETPIV
          ELSE
            HSPF = HSPF * PFRTPRE(KEYEQ)*(1-ETPIV) +
     *             HSPF * PFRTENH(KEYEQ)*ETPIV*
     *                 (1.0-ETPCAL(1,MY,KEYEQ,RVPHS,HSTEMP))
          ENDIF
C
C        Phase in the Tier 2 vehicles by introducing reduction 
C        factors for the 2004 and bigger years of pass-pass
C        light cars and trucks but only if the 
C        Tier 2 switch is on.
C
          IF(TIER2_FLAG==IYES)
     *               HSPF=HSPF*TR2FCTR(MY,ACTUAL_MY,IVGAS,KEYEQ)
C
C         Weight carbureted and fuel injected rates by corresponding sales fractions
C         For fuel injected vehices, adjust for the difference between in-use
C         and FTP fuel tank levels
C
C         All pre-1981 LDGV/T, HDGV and MC are 100% carbureted
C         so skip the fuel injection weighting.
C
          IF(IFDS.EQ.1) HSPF=(1.0-FINJ(1)-FINJ(2))*HSPF
          IF(IFDS.GT.1) HSPF=HSPF*FINJ(IFDS-1)*FULADJ
C
C         Combine the hot soak components together to form an average hourly
C         fleet value. This will be used to compute hourly hot soak emissions
C         below.
C
          HSBASE=HSBASE+HSPF
        END DO        !   End loop over fuel delivery systems
      END DO          !   End loop over KEYEQ
C
C     Initialize GLL total
C
          HSGLTOT=0.0
C
C     Sum over fuel systems
C
C         USER ENTERED GROSS LIQUID LEAKS
          IF(HSUSERGL_INPT==IYES) THEN
             HSGLTOT = (1.0-FINJ(1)-FINJ(2))*HSGLL_USER(1)*GLWT
             HSGLTOT = HSGLTOT+ FINJ(1)*HSGLL_USER(2)*GLWT
             HSGLTOT = HSGLTOT+ FINJ(2)*HSGLL_USER(3)*GLWT
          ELSE ! STANDARD MOBILE6 VALUES FOR LIQUID LEAKS
            DO IFDS= 1,3
              IF(IFDS.EQ.1) THEN
                HSGLTOT = HSGLTOT+(1.0-FINJ(1)-FINJ(2))*HSGLL(IFDS)*GLWT
              ELSE
                HSGLTOT = HSGLTOT+ FINJ(IFDS-1)*HSGLL(IFDS)*GLWT
              END IF
            END DO
          ENDIF
C
C         Add the gross liquid leaker contribution that is attributable to
C         HOT SOAK.
C
          HSBASE = HSBASE + HSGLTOT
C
C         Subtract off resting losses and the gross liquid leaker
C         contribution that is attributable to resting losses.
C         Both resting losses and HSBASE have
C         units of grams/hour at this point in the calculation.
C
          HSBASE=AMAX1(HSBASE - RSTLSS_EF,0.0)
C
C     Now, use the above base or hourly average to compute hourly hot soak
C     emissions. The remainder of the code uses the above hourly average to
C     produce hourly hot soak emissions from the fractions of vehicles which
C     have been running for a given length of time (soak time). These soak
C     times are in minutes and run from 1 to 60 minutes.
C
C     Go over the soak time distribution for this hour and add the
C     contributions for the soak time during this time interval.
C
C     Use the 14th time intervals' HSCOEF's for hourly intervals
C     14 through 24
      IHST = IH
      IF( IH .GT. 14 ) IHST = 14
C
C
C     If user-supplied hot soak activity rates have been specified, apply them.
C
      IF(HSACTFLAG.EQ.0) THEN   ! check for user supplied activity fractions
        DO IMIN = 1, 60
          CUM = HSCOEF(1,IHST,IWEEK) - HSCOEF(2,IHST,IWEEK) *
     &           EXP(-HSCOEF(3,IHST,IWEEK)*IMIN**HSCOEF(4,IHST,IWEEK))
          IF (CUM.LT.0.0) CUM = 0.0
          IF (IMIN.EQ.60) CUM = 100.0
          DIF = 0.01 * (CUM - PREV)
          PREV = CUM
          HMIN(IMIN) = HSBASE * HSFRACT(IMIN) * DIF
          HOTSOK = HOTSOK + HMIN(IMIN)
        END DO
      ELSE             ! user supplied activity fractions
        DO IMIN = 1, 60
          HMIN(IMIN) = HSBASE * HSFRACT(IMIN) * HSACTUSER(IMIN,IHST)
          HOTSOK = HOTSOK + HMIN(IMIN)
        END DO         !End of loop over IMIN with user supplied activity fractions
      ENDIF
C
      RETURN
      END
