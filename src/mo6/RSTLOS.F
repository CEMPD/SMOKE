      SUBROUTINE RSTLOS(IH, MY, ACTUAL_MY, ICY, IVGAS, INERR)
C
C  RSTLSHOUR computes the Resting Loss rates RSTLSS_EF  in grams/hour
C  for each model year and gas vehicle class at the given (hot soak)
C  temperature.
C
C  Called by EVAPHOUR.
C
C  Calls  EVWTPRE, EVWTENH, MAX and FLOAT.
C
C  Changes : (Last change first)
C
C  05 Apr 02 AIR Task Bug 388: modification of constants via QUITER
C  27 Feb 01 @EPA-BAg bug314 revised tier2 effects on heavy duty vehicles
C  31 Jan 01 @EPA-bag  low temperature resting losses
C  18 Jan 01 @EPA-bag  2007 rule for heavy duty vehicles
C  13 Dec 00 @EPA-bag  t-x6 MY and ACTUAL_MY added to parameter list for CAA switch
C  31 Oct 00 @EPA-bag  bug280 HDGV multiplicative factor applied to gas bus
C                      and ldgt3&4 conditional now includes all years less
C                      than 1979 (expand years for bug 25)
C  22 Aug 00 @EPA-bag  bug25 1972-1978 ldgt3&4 should use hdgv2b&3 parameters
C  21 Jul 00 @EPA-elg Added code to account for Tier2 Evap on IV=6 vehicles (2b's).
C  01 Jul 00 @EPA- BG  removed 40F temperature check at beginning of
C                      code
C  20 Jul 00 @EPA- BG  running loss; no motorcycle resting losses for
C                      Temp.<40F; replaced local phase-2 fractions with
C                      fractions from DIRNLDAT.FOR
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations
C  19 Jan 00 @Dyntel-JWR 1-017 Made RSTLSS_EF the hourly rest loss total
C  17 Jan 00 @Dyntel - JWR 1-017 Tier2 flag and phase-in fractions
C            derived from exhaust user-inputs
C  28 Dec 99 @Dyntel JWR Incorporated new weighting factors, included 
C            high altitude correction and heavy duty correction 
C  19 Nov 99 @Dyntel- JWR 1-017 Removed the 40-50 F condition- now the
C     linear T-dependence comes in at 40F. Also expanded to allow 4
C            pressure/purge test outcomes and extra model year range 
C            96+. Introduced array RSTLS2 for T-coefficient. Changed
C            liquid leaker to 9.16. Phase-in for 96-98 model years
C  13 Dec 99 @Dyntel-BAG 1-006 resting loss results are carried in RESTLS1.I
C  29 Dec 98 @Dyntel-BAG 2-664 Added motorcycle and heavy duty gas vehicle
C            hourly resting losses.
C  06 Nov 98 @Dyntel-BAG 2-000 Fixed missing model year ranges for
C            model years 1996 and later
C
C  Input on call:
C
C    argument list: IH, MY, ICY, IVGAS
C    common blocks:
C    /RESTLS/ RSTLS1,RSTLS2
C    /RESTLS1/ RSTLSS_EF
C    /TEMPS/  TEMHRLY
C    /VVTYPS/ VVLDGAS (IVTYPE.I)
C    /EVAPHD/ HDWGT
C    /EVAPAR/ FINJ
C    /EVAHIA/ HIADJ
C    /REGION/ IREJN
C    /IOUCOM/ IOUREP
C    /FLAGS5/ TIER2_FLAG 
C    /RESTLS1/ NO_RSTLOS
C    /TIER2/  T2_EVAP, T2_RED
C
C  Output on return:
C
C    argument list:
C    common blocks:
C    /RESTLS1/ RSTLSS_EF 
C
C
C  Local variable / array dictionary:
C
C   Name      Type              Description
C  ------     ----  -------------------------------------------------------
C  FDS_WGT     R     Sales fraction weighting factor used to sum resting loss
C                    emissions
C  GLLEF       R     Gross liquid leaker emission factor for resting losses in grams/hr.
C  RSTLSS_PP   R     Hourly resting loss array for pass/fail  pressure/purge 
C                    test vehicles with emissions controls
C  RSTPRE_PP   R     Hourly resting loss array for pass/fail  pressure/purge 
C                    test vehicles with emissions controls but pre-95
C  RSTENH_PP   R     Hourly resting loss array for pass/fail  pressure/purge 
C                    test vehicles with emissions controls but post-99
C  ALSS        R     Holds RSTL1() for calculation
C  BLSS        R     Holds RSTL2() for calculation
C  ALSSPRE     R     Holds RSTL1() for calculation of pre-1995 vhcles
C  BLSSPRE     R     Holds RSTL2() for calculation of pre-1995 vhcles 
C  ALSSENH     R     Holds RSTL1() for calculation of post-1999 vhcles
C  BLSSENH     R     Holds RSTL2() for calculation of post-1999 vhcles
C  RSTCARB     R     Array for the carbureted LD resting losses
C                    for use in calculating heavy duty losses
C  RSTFDS      R     Resting loss before weighting over FDS_WGT
C  RSTFDPRE    R     Resting loss before weighting over FDS_WGT (pre-95)
C  RSTFDENH    R     Resting loss before weighting over FDS_WGT (post-99)
C  LDCARB      R     Dummy variable for use in the HD calculation
C  HDRSTLSS    R     Heavy duty resting loss storage variable
C  HDRSTPH     R     Heavy duty resting loss storage variable (96-98)
C  RSTLSS_PP2  R     Hourly resting loss for pass pressure test vehicles (vapor leaks)
C                    vehicles without emissions controls
C  RSTLSS_FP2  R     Hourly resting loss for fail pressure test vehicles (vapor leaks)
C                    vehicles without emissions controls
C  RSTLSS      R     Sum of hourly resting loss for pass pressure test vehicles (vapor
C                    leaks) vehicles with emissions controls
C  RSTHDV      R     Holds resting losses for heavy-duty vehicles.
C  IVGAS_MC    I     Vehicle type label appropriate for motorcycles
C  HDWT        R     Overall multiplication factor for HD (1.5 or 2.0)
C  ICY         I     Calender year
C  IV          I     Vehicle type
C  IT          I     Index for initialization of weighting factors
C  MY          I     Model year
C  MAGE        I     Age (JDX type i.e. CY-MY +1)
C  PPST        I     Test-status label (1-4)
C  QPST        I     Test-status label (1-4)
C  QQST        I     Test-status label (1-4)
C  RSTFI       R     Stores fuel-injected rest loss for use by HDGV
C  RSTFIPRE    R     Stores FI rest loss for use by HDGV (pre-95)
C  RSTFIENH    R     Stores FI rest loss for use by HDGV (post-99)
C  RSTLSS_PR72 R     Resting loss for pre-control vehicles 
C  RIFDS       I     Fuel delivery system loop variable. IFDS is mapped into
C                    two fuel delivery systems, either carbureted
C                    (RIFDS=1) or  fuel injected (RIFDS=2)
C  MYRNG       I     Model year range  label
C  TMRNG       I     Temperature/model year range used to categorize
C                    calculation method
C  TMPRTURE    R     Temperature for hour, IH, in degrees Fahrenheit
C  EWT         R     Stores pressure/purge weighting factors
C  FAC         R      Overall reduction factor for Tier2 emissions
C
C  Notes:
C
C  RSTLSHOUR was added for MOBILE6.
C
      USE DIRNLDAT, ONLY : W95,W99
      IMPLICIT NONE
C
      INCLUDE 'EVAPAR.I'
      INCLUDE 'EVAPHD.I'
      INCLUDE 'EVAHIA.I'
      INCLUDE 'EVAHS1.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESTLS.I'
      INCLUDE 'RESTLS1.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'FLAGS5.I'
C
C     Declare external functions
C
      REAL TR2FCTR
      REAL EVWTPRE
      REAL EVWTENH
      REAL GLL_WGT
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: ICY
      INTEGER, INTENT(IN) :: IH
      INTEGER, INTENT(IN) :: IVGAS
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: ACTUAL_MY
      INTEGER, INTENT(INOUT) :: INERR
C
C     Declare local variables
C
      INTEGER, SAVE :: PRE_CNTRL = 2
      INTEGER, SAVE :: CNTRL     = 3
      INTEGER, SAVE :: IYES      = 1
      INTEGER, SAVE :: INO       = 0
      INTEGER, SAVE :: IGASBUS   = 15
      INTEGER MYRNG
      INTEGER RIFDS
      INTEGER PPST
      INTEGER QPST
      INTEGER QQST
      INTEGER TMRNG
      INTEGER MAGE
      INTEGER IT
      INTEGER IUSE2B
      INTEGER IVQTR
C
      REAL FDS_WGT
      REAL TMPRTURE
      REAL RSTLSS_PP(4)
      REAL RSTCARB(4)
      REAL RSTLSS
      REAL HDRSTLSS
      REAL RSTLSS_PP2
      REAL RSTLSS_FP2
      REAL RSTLSS_PR72 
      REAL ALSS
      REAL BLSS
      REAL RSTFDS
C
      REAL EWT(5)
      REAL RSTPRE_PP(4)
      REAL RSTENH_PP(4)
      REAL RSTLSSPH
C
      REAL ALSSPRE
      REAL BLSSPRE
      REAL RSTFDPRE
      REAL ALSSENH
      REAL BLSSENH
      REAL RSTFDENH
      REAL HDWT
      REAL LDCARB
      REAL FAC
      REAL RSTHDV(4,2,3)
      REAL RMY
      REAL, SAVE    :: GLLEF  = 9.16
C
C     Initializations and Parameter Determinations
C
      RSTLOS_FR_RNLOSS = 0.0
      RSTLSS_PP    = 0.0
      RSTLSS       = 0.0
      RSTFDS       = 0.0
      RSTLSS_EF    = 0.0
      RSTCARB      = 0.0
      MAGE         = 1
      IUSE2B       = INO
      IVQTR        = 0
      FAC          = 0.0
      ALSS         = 0.0
      BLSS         = 0.0
      RSTFDS       = 0.0
      ALSSPRE      = 0.0
      BLSSPRE      = 0.0
      RSTFDPRE     = 0.0
      ALSSENH      = 0.0
      BLSSENH      = 0.0
      RSTFDENH     = 0.0
      RSTPRE_PP    = 0.0
      RSTENH_PP    = 0.0
      RSTLSSPH     = 0.0
      RMY          = 0.0
C
      TMPRTURE = TEMHRLY(IH)
      IF(NO_RSTLOS==IYES ) GOTO 99
      IF( TMPRTURE .GT. 105. ) TMPRTURE = 105.
C
C     Set vehicle age (JDX type age)
      MAGE=ICY-ACTUAL_MY+1
C
      DO IT=1,5
        EWT(IT)=0.0
C       Set the evaporative weights according to the model year
C       for years outside the phase-in 96-98
        IF(MY.LE.1995) EWT(IT)=EVWTPRE(MY,ACTUAL_MY,MAGE,IT,2)
        IF(MY.GE.1999) EWT(IT)=EVWTENH(MY,ACTUAL_MY,MAGE,IT,2)
      END DO
C
C     Motorcycles
      IF(IVGAS==14) THEN
C
C       Motorcycle resting losses are discussed on page 26
C       of M6.EVP.002.)
C
C         resting losses for temperatures greater than 40F
          RSTLSS=0.0
          IF(TMPRTURE > 40.) RSTLSS=0.04434 +
     *             0.006134*(FLOAT(MAGE)-1.) + 0.000859*TMPRTURE
C
C       Gross liquid  leakers
C
          RSTLSS=RSTLSS*(1.0 - GLL_WGT(MY,ACTUAL_MY,2,MAGE-1))
     &     + GLL_WGT(MY,ACTUAL_MY,2,MAGE-1) *GLLEF
C
C      For MCs set the total
C
          RSTLSS_EF = RSTLSS
C
          GOTO 98
C
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     All other gas vehicles
C     Temperature/Model Year Ranges :
C
      IF( MY<=1971 )  THEN
        TMRNG = PRE_CNTRL !2
      ELSE IF( MY>1971 )  THEN
        TMRNG = CNTRL     !3
      END IF
C
      IF( TMPRTURE<=40.0 ) TMRNG = 1
C
C     Model Year Ranges for parameter determinations
C     (used only for model years after 1971):
C
      IF( MY<1980 ) THEN
        MYRNG = 1
        IF((IVGAS==4.OR.IVGAS==5).AND.(MY<1979)) IUSE2B=IYES !1972-1978 LDGT3&4 use
                                                             !hdgv2b&3 parameters
      ELSE IF( MY>=1980 .AND. MY<=1985 ) THEN
        MYRNG = 2
      ELSE IF( MY>=1986 .AND. MY<=1995 ) THEN
        MYRNG = 3
      ELSE IF( MY>=1996) THEN
        MYRNG = 4
      ENDIF
C
C     RESTING LOSSES :
C     Select resting loss calculation depending on the model year/
C     temperature range
C
      MODEL_YR_TMP_RNG : SELECT CASE(TMRNG)
C
C ******************************************************************
C       Resting loss is zero if any evaporative temperature is 40F or less,
C       (Except for liquid leakers)
C
        CASE(1) MODEL_YR_TMP_RNG
C
C         Temperature too low to get vapour emissions
C
          RSTLSS_PP(1) = 0.0
          RSTLSS_PP(2) = 0.0
          RSTLSS_PP(3) = 0.0
          RSTLSS_PP(4) = 0.0
C
C         Add liquid leakers
C
          RSTLSS_EF= GLL_WGT(MY,ACTUAL_MY,2,MAGE-1) * GLLEF
          GOTO 98
C
C ******************************************************************
C           
C       Vehicles With Emissions Controls, i.e., model years
C       1972 and later
C
        CASE(3) MODEL_YR_TMP_RNG
C
          RSTLSS_PP(1) = 0.0
          RSTLSS_PP(2) = 0.0
          RSTLSS_PP(3) = 0.0
          RSTLSS_PP(4) = 0.0
          RSTPRE_PP(1) = 0.0
          RSTENH_PP(1) = 0.0
          RSTPRE_PP(2) = 0.0
          RSTENH_PP(2) = 0.0
          RSTPRE_PP(3) = 0.0
          RSTENH_PP(3) = 0.0
          RSTPRE_PP(4) = 0.0
          RSTENH_PP(4) = 0.0
C
C         Loop over the four pressure/purge test outcomes
C         1=pass pressure pass purge
C         2=pass pressure fail purge
C         3=fail pressure pass purge
C         4=fail pressure fail purge
C
          DO 500 PPST =1,4
C
C           Since the fit parameters are fuel delivery system dependent
C           need a loop over fuel injected and carbureted vehicles
C           Fuel System Categorization (used only for model years after 1972):
C           (fuel delivery system type)
C
            FDS_WGT = 0.0
            DO RIFDS = 1,2
C
C           Set the weighting factor, FDS_WGT, for carbureted and fuel injected
C           rates by corresponding sales fractions. The weighting factor is used
C           below for adding resting loss results together to obtain the fleet
C           resting loss rate.
C
              IF(RIFDS==1) FDS_WGT = 1.0 - FINJ(1) - FINJ(2)
              IF(RIFDS==2) FDS_WGT = FINJ(1) + FINJ(2)
C
              ALSS = RSTLS1(MYRNG,PPST,RIFDS)
              BLSS = RSTLS2(MYRNG,PPST)
C
C             Phase in the post 95 coefficients gradually
C             according  to the prescriptions
              IF(MY==1996) THEN
                ALSSPRE = RSTLS1(3,PPST,RIFDS)*W95(MY-1995)
                ALSSENH = RSTLS1(4,PPST,RIFDS)*W99(MY-1995)
C
                BLSSPRE = RSTLS2(3,PPST)*W95(MY-1995)
                BLSSENH = RSTLS2(4,PPST)*W99(MY-1995)
C
              ELSE IF(MY==1997) THEN
                ALSSPRE = RSTLS1(3,PPST,RIFDS)*W95(MY-1995)
                ALSSENH = RSTLS1(4,PPST,RIFDS)*W99(MY-1995)
C
                BLSSPRE = RSTLS2(3,PPST)*W95(MY-1995)
                BLSSENH = RSTLS2(4,PPST)*W99(MY-1995)
C
              ELSE IF(MY==1998) THEN
                ALSSPRE = RSTLS1(3,PPST,RIFDS)*W95(MY-1995)
                ALSSENH = RSTLS1(4,PPST,RIFDS)*W99(MY-1995)
C
                BLSSPRE = RSTLS2(3,PPST)*W95(MY-1995)
                BLSSENH = RSTLS2(4,PPST)*W99(MY-1995)
              END IF
C
C        Phase in the Tier 2 vehicles by introducing reduction
C        factors for the 2004 and bigger years of pass-pass
C        if the Tier 2 switch is on.
C
              IF(TIER2_FLAG==IYES) THEN
                FAC = TR2FCTR(MY,ACTUAL_MY,IVGAS,PPST)
                ALSS=ALSS*FAC
                BLSS=BLSS*FAC
              END IF
C                              
C             Calculate resting loss prior to weighting
C             over fuel system
C
              RSTFDS= ALSS + BLSS*TMPRTURE
              RSTFDPRE= ALSSPRE + BLSSPRE*TMPRTURE
              RSTFDENH= ALSSENH + BLSSENH*TMPRTURE
C
C            But if the emission rate becomes negative, make it 0.
C
              RSTFDS = MAX(0.0,RSTFDS)
              RSTFDPRE=MAX(0.0,RSTFDPRE)
              RSTFDENH=MAX(0.0,RSTFDENH)
C
C            Store the carbureted and fuel-injected resting losses
C            for use in the heavy duty calculation later
C
              IF(RIFDS==1) RSTCARB(PPST)=RSTFDS
              IF(RIFDS==2) THEN
                RSTHDV(PPST,RIFDS,2)=RSTFDPRE      !HD Pre Enhanced RSTLOS
                RSTHDV(PPST,RIFDS,3)=RSTFDENH      !HD Enhanced RSTLOS
              END IF
              RSTHDV(PPST,RIFDS,1)=RSTFDS          !HD Regular RSTLOS
C
C            Carry out the weighting over fuel system
C
              RSTLSS_PP(PPST) = RSTLSS_PP(PPST)+FDS_WGT*RSTFDS
              RSTPRE_PP(PPST) = RSTPRE_PP(PPST)+FDS_WGT*RSTFDPRE
              RSTENH_PP(PPST) = RSTENH_PP(PPST)+FDS_WGT*RSTFDENH
C
            END DO  !RIFDS Loop
C
 500      END DO !PPST Loop
C
C         Do the sum over the four pressure/purge pass/fail
C         categories using the externally calculated
C         weighting factors for the diurnal & resting losses
C
          DO QPST=1,4
C
            RSTLSS=RSTLSS +EWT(QPST)*RSTLSS_PP(QPST)
C
C           Calculate the phase-in total
C
            IF(MY<=1998 .AND. MY>=1996) THEN
               RSTLSSPH = RSTLSSPH +
     &             EVWTPRE(MY,ACTUAL_MY,MAGE,QPST,2)*RSTPRE_PP(QPST) +
     &             EVWTENH(MY,ACTUAL_MY,MAGE,QPST,2)*RSTENH_PP(QPST)
            END IF
          END DO
C
C         Add the gross liquid leaker term
          RSTLSS=RSTLSS + GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)*GLLEF
          IF(MY>=1996 .AND. MY<=1998)
     *         RSTLSSPH=RSTLSSPH+GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)*GLLEF
C
C         Pre-Control Vehicles (i.e., model years 1972 and older):
C         EPA parametrization of pre-control vehicles has negative
C         emissions results. Fix is the same as for diurnals. Use
C         parameters from the above parametrization of vehicles with
C         emissions controls as a check. The larger of the pre-control
C         and the above control emissions will be used.
C
              
        CASE(2) MODEL_YR_TMP_RNG
C
C         Calculate what the 1972  emissions would be for
C         pressure pass and pressure failing vehicles
C
          RSTLSS_PP2 = 0.0553 + 0.002812*TMPRTURE
          RSTLSS_FP2 = 0.07454 + 0.002812*TMPRTURE
          RSTLSS_PR72 = 0.0
C
C         for all test-status categs and fuel systems we have
C         to compare with the formula for pre-72 model years
C
          RSTLSS_PR72 = -0.76844+0.002812*TMPRTURE +
     &                   0.040528*FLOAT(ICY - MY)
C
C    Examine the two ways of calculating TMRNG=2
C    vehicles  and take whichever is the larger (avoids negatives)              
C
          RSTLSS_PP2 = MAX(RSTLSS_PP2,RSTLSS_PR72)
          RSTLSS_FP2 = MAX(RSTLSS_FP2,RSTLSS_PR72)
C
C    Store the  pressure fail results for use by heavy dutys
C
          RSTCARB(3)= RSTLSS_FP2
C
C    Weight these together using the new weighting factors
C    and include the gross liquid leakers
C
          RSTLSS=RSTLSS_PP2*(EWT(1) +EWT(2)) +
     &               RSTLSS_FP2*(EWT(3) + EWT(4)) + GLLEF*EWT(5)
C
C ******************************************************************
C
      END SELECT MODEL_YR_TMP_RNG
C
C ******************************************************************
C
      RSTLSS_EF = RSTLSS
C
      IF(MY<=1998 .AND. MY>=1996) RSTLSS_EF = RSTLSSPH
C
C     Heavy duty gas vehicles are about 1.56 times larger than resting
C     losses for light duty gas trucks. (See page 24 of M6.EVP.002.)
C     IVGAS = 1,2,3,4,5 : light duty gas vehicles
C     IVGAS = 6,7,8,9,10,11,12,13,15 : heavy duty gas vehicles and gas buses
C     IVGAS = 14 : motorcycles
C     VVHDGV(6)=VVHDGV(7)=VVHDGV(8)=VVHDGV(9)=VVHDGV(10)=VVHDGV(11)=VVHDGV(12)=VVHDGV(13)=1
C     VVHDGV(1-5,14,15)=0
C
      IF(VVHDGV(IVGAS)==IYES .OR. IUSE2B==IYES .OR. IVGAS==IGASBUS) THEN !heavy duty vehicles, buses
C                                                                         and 1972-1978 ldgt3&4
C     Initialize the temporary storage variables for HD resting loss
C           
        HDRSTLSS=0.0
        LDCARB=0.0
        RSTLOS_FR_RNLOSS=0.0
C
C     (Actually HDGV are 1.5 or 2.0 times the relevant LD emissions)
        IF(IVGAS<8) THEN !heavy duty vehicles2b&3 and
          HDWT=HDWGT(1)  !1972-1978 ldgt3&4
        ELSE
          HDWT=HDWGT(2)
        END IF
C     Loop over the pressure/purge test-status bins
C
        DO QQST=1,4
C
C       For model years before 1985 the relevant LD carbureted
C       vehicles are the ones that failed the pressyre test,
C       otherwise use the proper test-status. For post 96
C       vehicles use Fuel injected LD emission factors
C
          IF(MY<1985) THEN
            LDCARB=RSTCARB(3)
          ELSE
            LDCARB=RSTHDV(QQST,1,1)*(1.0 - FINJ(1) - FINJ(2)) +
     *             RSTHDV(QQST,2,1)*(FINJ(1) + FINJ(2))
          END IF
C
C         RSTHDV(QQST,RIFDS,IE)  ==  RIFDS - 1:CARB, 2:FI
C                                    IE - 1:STANDARD, 2:PRE-ENHANCED, 3:ENHANCED
C                                    2 and 3 are used only for the 1996-98 MYs.
C
C     For the phase-in years use the fuel-injected press/purge LDs
C
          IF(MY>=1996 .AND. MY<=1998) THEN
            HDRSTLSS=HDRSTLSS+HDWT*(EVWTPRE(MY,ACTUAL_MY,MAGE,QQST,2)*
     &        RSTHDV(QQST,2,2) +
     &        EVWTENH(MY,ACTUAL_MY,MAGE,QQST,2)*RSTHDV(QQST,2,3))
          ELSE
            HDRSTLSS = HDRSTLSS + HDWT*EWT(QQST)*LDCARB
          END IF
C
          IF(IVGAS>=6.AND.MY>2007.AND.HDR2007==IYES
     *       .AND.TIER2_FLAG==IYES) THEN
            SELECT CASE (QQST)
              CASE (1)     !QQST=IPASS_BTH
                IF(FAC<=0.000001) THEN
                  IVQTR=IVGAS
                  IF(IVGAS==14) IVQTR=24
                  IF(IVGAS==15) IVQTR=25
                  RMY=REAL(MY)
                  CALL QUITER(RMY,IVQTR,53,INERR)
                  GOTO 99
                ENDIF
                RSTLOS_FR_RNLOSS= HDRSTLSS/FAC
              CASE DEFAULT !QQST=2-4
                RSTLOS_FR_RNLOSS=RSTLOS_FR_RNLOSS+HDWT*EWT(QQST)*LDCARB
            END SELECT
          ENDIF
C
        END DO      !End QQST Loop
C
C     Finally add on the unenhanced gross liquid leaker term
C
        RSTLOS_FR_RNLOSS=RSTLOS_FR_RNLOSS +
     *           GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)*GLLEF
        HDRSTLSS=HDRSTLSS  + GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)*GLLEF
C
C     Set it equal to the resting loss total
C
        RSTLSS_EF =HDRSTLSS
C
      END IF
C
C     Adjust the emissions if the region under consideration is a
C     high altitude region. (See page 36 of M6.EVP.001.)
C
  98  IF(IREJN==2) THEN
        RSTLSS_EF=HIADJ*RSTLSS_EF -
     *            0.3*GLLEF*GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)
        RSTLOS_FR_RNLOSS=HIADJ*RSTLOS_FR_RNLOSS -
     *            0.3*GLLEF*GLL_WGT(MY,ACTUAL_MY,2,MAGE-1)
C
        RSTLSS_EF=MAX(0.0,RSTLSS_EF)
        RSTLOS_FR_RNLOSS=MAX(0.0,RSTLOS_FR_RNLOSS)
C
      END IF
C
  99  RETURN
      END
