      SUBROUTINE EVP_HRLY(ICY,IVGAS,IDX,IV,JDX,MY,ACTUAL_MY,INERR)
C
C  EVAP_HRLY calls the subroutines and functions to which determine hydrocarbon
C  evaporative emissions for hot soak, diurnal, crankcase, resting losses, and
C  running losses. (Hourly refueling losses have been determined in REFUEL.FOR
C  called by EFCALX.FOR.) Hourly emissions are determined from the calls to
C  the various individual evaporative routines. Sums for daily emission values
C  are performed and units conversions are made.
C
C  Called by:  HCCALX
C
C  Calls:  CRANKC, DB_EVAP, DIRNL_HR, DIRNLDAY, DIURATIO, FUEL2, HOTSOK, IEVPTR,
C          PCTLEV, QUITER, RNLOSSSUM, RSTLOS, CNVRT_U
C
C  Changes:  (Last change first)
C
C  21 Dec 01 M6011221: Changes to CNVRT_U parameters.
C  15 Nov 01 AIR Task 27: Added INERR processing due to PMFLECON update.
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch
C  01 Nov 00 @EPA-bag BUG296 add IVGAS to DIRNL_HR parameter list
C  06 Oct 00 @EPA-bag BUG137 Remove RN_FLAG
C  20 Sep 00 @EPA-elg Altered the ZEV and NGVFRAC calculations.  Bug 257 fix.
C  08 Sep 00 @EPA-bag bug252&253; removed IV from DIRNL_HR parameter list
C  08 Aug 00 @EPA-elg Added NO_EVAPCALC to shut down evap calculations.
C  02 Aug 00 @EPA-djb Bug 212: Initialized running losses for motorcycles & comments.
C  24 Jul 00 @EPA-bag D15 Moved call to multiday diurnals to EVP_HRLY.FOR; multi-day
C            ratios are calculated once for each MY and IV and then passed through
C            DIRNLDAT.FOR
C  21 Jul 00 @EPA-bag running loss; Added ICY and JDX to RNLOSSSUM  parameter list.
C  28 Jun 00 @EPA-ddj D15   Added MY and IV to DIRNL_HR parameter list.
C 23 Jun 00 @EPA-bag bug52  No refueling losses for motorcycles
C 26 Jun 00 @EPA-bag bug129 Running losses have database output for freeways,
C                    ramps, arterials, and local road types (removed all
C                    roads category)
C   2 Jun 00 @EPA-ddj E11, Subtracted the fraction of Natural Gas
C            Vehicles from the fractions of non-zero emission vehicles
C            (variable Z).  Added HDGVs to calculation of FINJ.
C  30 mar 00 @DynTel-mla  1-008 database output
C  30 Jan 00 @DynTel-bg  1-005 Calls to new diurnal subroutines and functions; added
C            re-initializations of FINJ's
C  25 Feb 00 @DynTel-bg  1-006 divided crankcase gram/hour sum by 24 (this quantity,
C            EVAP_GM_HR(ICRNCASE), has not been requested by EPA to be output;
C            it has only been added for completeness
C  12 Jan 00 @DynTel-bg  1-006 Applied HVMT fractions to average daily emissions (hourly
C            sums) which are in units of g/mi
C  05 Jan 00 @DynTel-bg  1-006 Applied zero emission fractions to crankcase tampering
C            effects, TCC
C  10 Dec 99 @DynTel-bg  1-006 Renamed to EVAP_HRLY and changed evaporative hourly and
C            daily emisions arrays. Added a call to the function, CNVRT_U.FOR, to
C            do unit conversions.
C  10 Jun 99 @DynTel-bg  0-632 Fuel injection vehicle fractions, FINJ, now initialized
C            before DIRNLDAY is called. The reverse order caused incorrect fractions to be
C            used in DIRNLDAY.
C  18 May 99 @DynTel-HQ  2-650 removed conditionals for determining MPD and TPD
C   6 May 99 @DynTel-HQ  2-664 call the rnlosssum instead of rnloss
C  26 Apr 99 @DynTel-HQ  2-698 renamed BDSPEED.I to SPEED9.I
C  06 Apr 99 @DynTel-bag 2-000 MAXIH_EVAP replaced with MAXIH
C  30 Mar 99 @DynTel-bag  Had running loss hourly sum in an IF statement
C            with argument depending on VVGAS(IVGAS) different that of an
C            IF statement for call to RNLOSS with argument VVGAS(IV)
C            resulting with running loss sum missing gas bus running losses
C  3  Mar 99 @DynTel-bag  Renamed to CCEVRT and added hourly hot soak
C            and running loss code.
C  30 Oct 98 @DynTel-bag  Renamed to EVAPHOUR and added hourly diurnal
C                         and resting loss emissions.
C
C  Input on call:
C
C    Constants: 
C    MAXIH from MAXIMA.I
C
C    argument list: ICY,IVGAS,IDX,IV,JDX
C    function : CRANKC
C    common blocks:
C
C    /CEVBMY/        BMYMPD
C    /EVAPAR/        IFDS,ISTDC,ISTDT,ISTDP,FINJ  (to IEVPTR, HOTSOK)
C    /INJECT/        TBI,PFI
C    /IVTYPE/        VTMC,VVGAS
C    /MAXIMA/        MAXIH
C    /NGV/           NGVFRAC
C    /RESTLS1/       RSTLSS_EF, NO_EVAPCALC
C    /RLCOM3/        RLRATE
C    /SPEED9/        HVMT
C    /TAMOUT/        TCC
C
C  Output on return:
C
C    argument list:
C    common blocks:
C    /DAILY_EVAP/    EVAP_GM_MI(6),EVAP_GM_HR(6),EVAP_GM_TRP,EVAP_GM_GA
C    /RESTLS1/       RSTLSS_EF
C  Local variable / array dictionary:
C
C   Name        Type              Description
C  ------       ----  -------------------------------------------------------
C  ICY           I    vehicle calendar year; user input
C  IDX           I    vehicle yearly loop index(1-25 years)
C  JDX           I    vehicle age index; JDX=25-IDX+1
C  MY            I    Vehicle model year; MY=ICY-JDX+1=ICY+IDX-25
C  IV            I    Vehicle type loop control variable
C  IVGAS         I    gasoline fueled vehicle type index
C  CC_GM_HR      R    crankcase emissions in grams per hour
C  CC_GM_MI      R    crankcase emissions in grams per mile
C  DRNL          R    hourly diurnal losses
C  DRNL_GM_HR    R    diurnal losses in grams per hour
C  DRNL_GM_MI    R    diurnal losses in grams per mile
C  FCRNCASE      I    database output indice (7 : crankcase emissions)
C  FDIURNAL      I    database output indice (4 : diurnal emissions)
C  FHOTSOAK      I    database output indice (3 : hot soak emissions)
C  FREFUEL       I    database output indice (8 : refueling emissions)
C  FRESTLOSS     I    database output indice (5 : resting loss emissions)
C  FRUNLOSS      I    database output indice (6 : running loss emissions)
C  HS_GM_HR      R    hot soak emissions in grams per hour
C  HS_GM_MI      R    hot soak emissions in grams per mile
C  HS_GM_TRP     R    hot soak emissions in grams per trip
C  HSK           R    hourly hot soak losses; HOTSOK + tampering hs offset, if any
C  IALLROADS     I    roadway index (5 : all roads weighted together)
C  ICRNCASE      I    evaporative emissions index (4=ICRNCASE) -  crankcase
C  IDIURNAL      I    evaporative emissions index (2=IDIURNAL) -  diurnal
C  IFREEWAY      I    roadway index (1 : freeways)
C  IG            I    pointer to fuel injected sales percentage for given my
C  IH            I    hourly index loop (hours 1-24, 6am-5am)
C  IHOTSOAK      I    evaporative emissions index (1=IHOTSOAK) -  hotsoak
C  IRAMP         I    roadway index (4 : freeway ramps)
C  IREFUEL       I    evaporative emissions index (6=IREFUEL)  -  refuelng loss
C  IRESTLOSS     I    evaporative emissions index (5=IRESTLOSS)-  resting loss
C  IROAD         I    roadway index loop
C                            1 : freeways
C                            2 : arterial
C                            3 : local
C                            4 : ramp
C                            5 : all roads weighted together
C  IRUNLOSS      I    evaporative emissions index (3=IRUNLOSS) -  running loss
C  IVTAM         I    index for vehicles with tampering
C  IVTEMP        I    index for fuel injection
C  MAX_DIURNL_HR I  (=18)diurnal emissions end at midnight
C  MY            I    Vehicle model year; MY=ICY-JDX+1=ICY+IDX-25
C  RFL_GM_GA     R    refueling emissions in grams per gallon
C  RFL_GM_HR     R    refueling emissions in grams per hour
C  RFL_GM_MI     R    refueling emissions in grams per mile
C  RN_GM_HR      R    running loss emissions in grams per hour
C  RN_GM_MI      R    running loss emissions in grams per mile
C  RSTLS_GM_HR   R    resting loss emissions in grams per hour
C  RSTLS_GM_MI   R    resting loss emissions in grams per mile
C  RUNLOSS       R    hourly running losses by road type
C  ZEV           R    zero emission vehicle fractions
C
C
C  Notes:
C
C  Hot soak, diurnal and crankcase without tampering and with tampering HC are
C  figured separately and then added together, according to Evap Study 2.
C  The additions are done before miles per day and trips per day are applied
C  as weights.
C
C  CCEVRT was modified in MOBILE4.1 to transform hot soaks, diurnals
C  and partial diurnals to grams per mile (g/mi), when IV is LDGV/T
C  and MPD is positive.    Also, oxygenated and cert. fuels support
C  from MOBILE4 was removed in MOBILE4.1.
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C
C  (Sep-19-1994) @ CSC ked (for bsg) Request 438 Fix 1990 and earlier calendar
C  year removed IY from the argument list, and from the calls to HOTSOK and
C  DIURNL
C
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C  29 Oct. 1997 Dyntel@MMS 0-632a deleted condition for calling PCTLEV
C
C
      IMPLICIT NONE
C
      INCLUDE 'CEVBMY.I'
      INCLUDE 'DAILY_EVAP.I'
      INCLUDE 'EVAPAR.I'
      INCLUDE 'INJECT.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'NGV.I'
      INCLUDE 'RLCOM3.I'
      INCLUDE 'RESTLS1.I'
      INCLUDE 'SPEED9.I'
      INCLUDE 'TAMOUT.I'
      INCLUDE 'UNIT_PARAM.I'
C
C     Declare parameters.
C
      INTEGER, INTENT(IN) :: ICY       ! calender year
      INTEGER, INTENT(IN) :: IDX       ! Model year index (chronological order)
      INTEGER, INTENT(IN) :: IV        ! vehicle type index (see BD06.FOR)
      INTEGER, INTENT(IN) :: IVGAS     ! gasoline fueled vehicle type index
      INTEGER, INTENT(IN) :: JDX       ! Model year index (reverse chronological order)
      INTEGER, INTENT(IN) :: MY        ! Model year corrected if NO CLEAN AIR ACT (NO CAA) is enabled
      INTEGER, INTENT(IN) :: ACTUAL_MY ! Model year (uncorrected; ACTUAL_MY=ICY+IDX-MAXYRS)
      INTEGER, INTENT(INOUT) :: INERR  ! Errors
C
C     Declare external functions.
C
      INTEGER IEVPTR
      REAL    CRANKC
      REAL    HOTSOK
      REAL    PCTLEV
      REAL    CNVRT_U
      REAL    DIRNL_HR
C
C     Declare local variables.
C
      REAL     DRNL        ! hourly diurnal losses
      REAL     HSK         ! hourly hot soak losses; HOTSOK + tampering hs offset, if any
      REAL     RUNLOSS(5)  ! running losses indexed by 5 roadway types
      REAL     ZEV         ! discount evap emission by ZEV fraction
      INTEGER  IH          ! Hourly loop control variable
      INTEGER  IG5         ! pointer to fuel injected sales percentage for given my   
      INTEGER  IROAD       ! road type index
                           ! 1 : freeways
                           ! 2 : arterial
                           ! 3 : local
                           ! 4 : ramp
                           ! 5 : all roads weighted together
      INTEGER  IVTAM       ! index for vehicles with tampering
      INTEGER  IVTEMP      ! index for fuel injection
C
      INTEGER, PARAMETER :: MAX_DIURNL_HR=18 ! diurnal emissions end at midnight
C                       evap. emission indices
      INTEGER, PARAMETER :: IHOTSOAK=1
      INTEGER, PARAMETER :: IDIURNAL=2
      INTEGER, PARAMETER :: IRUNLOSS=3
C                       evap. emission indices
      INTEGER, PARAMETER :: ICRNCASE=4
      INTEGER, PARAMETER :: IRESTLOSS=5
      INTEGER, PARAMETER :: IREFUEL=6
C                       three of the five roadway type indices
      INTEGER, PARAMETER :: IFREEWAY=1
      INTEGER, PARAMETER :: IRAMP=4
      INTEGER, PARAMETER :: IALLROADS=5
C                      evap. emission indices for database output
      INTEGER, PARAMETER :: FHOTSOAK=3
      INTEGER, PARAMETER :: FDIURNAL=4
      INTEGER, PARAMETER :: FRESTLOSS=5
      INTEGER, PARAMETER :: FRUNLOSS=6
      INTEGER, PARAMETER :: FCRNCASE=7
      INTEGER, PARAMETER :: FREFUEL=8
C
C     Initialize and define the hourly evaporative emissions arrays.
C
      REAL, DIMENSION(24)   :: DRNL_GM_HR
      REAL, DIMENSION(24)   :: DRNL_GM_MI
      REAL, DIMENSION(24)   :: HS_GM_TRP
      REAL, DIMENSION(24)   :: HS_GM_MI
      REAL, DIMENSION(24)   :: HS_GM_HR
      REAL, DIMENSION(24)   :: CC_GM_HR
      REAL, DIMENSION(24)   :: CC_GM_MI
      REAL, DIMENSION(24,5) :: RN_GM_MI
      REAL, DIMENSION(24,5) :: RN_GM_HR
      REAL, DIMENSION(24)   :: RFL_GM_MI
      REAL, DIMENSION(24)   :: RFL_GM_HR
      REAL, DIMENSION(24)   :: RFL_GM_GA
      REAL, DIMENSION(24)   :: RSTLS_GM_MI
      REAL, DIMENSION(24)   :: RSTLS_GM_HR
C
      DRNL_GM_HR =0.0
      DRNL_GM_MI =0.0
      HS_GM_TRP  =0.0
      HS_GM_MI   =0.0
      HS_GM_HR   =0.0
      CC_GM_HR   =0.0
      CC_GM_MI   =0.0
      RFL_GM_MI  =0.0
      RFL_GM_HR  =0.0
      RFL_GM_GA  =0.0
      RSTLS_GM_MI=0.0
      RSTLS_GM_HR=0.0
      RN_GM_MI   =0.0
      RN_GM_HR   =0.0
C
C     Initialize daily evaporative emissions arrays.
C
      EVAP_GM_MI =0.0
      EVAP_GM_HR =0.0
      EVAP_GM_TRP=0.0  !used for hot soak emissions only
      EVAP_GM_GA =0.0  !used for refueling emissions only
C
      IF(NO_EVAPCALC.EQ.1) GOTO 99  !Shuts off all Evap calculations and output.
C
C     Warn and exit evaporative calculations if the miles per day is zero.
C
      IF(BMYMPD(JDX,IV).LE.0.0) THEN
        CALL QUITER(FLOAT(MY),IVGAS,58,0)
        GOTO 99
      ENDIF
C
C     This sets the combined fraction of NGVs and ZEVs.
C     If both are present at the same time, the model first
C     sets aside the ZEV percentage from the total fleet.
C     The NGVFRAC is then applied to the remaining non ZEV fleet
C     fraction.
C
      ZEV = (1.0 - NGVFRAC) * (1.0 - PCTLEV(MY,IV,11))
C
C     Look up model year group pointers:
C     Carb
      IFDS=1
      ISTDC=IEVPTR(MY,IVGAS,1)
      FINJ=0.0   ! initialize all vehicles to carbureted vehicles
      IF(VVGAS(IV).EQ.1.AND.MY.GE.1981) THEN
        IVTEMP = IVGAS
        IF(IV.EQ.VTGASBUS) IVTEMP = 14
C
C       TBI
        IFDS=2
        ISTDT=IEVPTR(MY,IVGAS,1)
C
        IG5=IEVPTR(MY,IVGAS,2)
        FINJ(1)=TBI(IG5,IVTEMP)
C
C       PFI
        IFDS=3
        ISTDP=IEVPTR(MY,IVGAS,1)
C
        IG5=IEVPTR(MY,IVGAS,2)
        FINJ(2)=PFI(IG5,IVTEMP)
      ENDIF
C
C     DAILY AVERAGE CRANK CASE EMISSIONS:
C     Find the crankcase evaporative losses in units of grams/mile; crank
C     case losses are zero when miles per day are zero. (Only MOBILE5 daily
C     crank case emissions are calculated.)
C
      EVAP_GM_MI(ICRNCASE) = ZEV*CRANKC(MY,IV)
C
C     GV/T => tampering crankcase HC effects are to be added.
C     Note that IVGAS, used by the evap code, includes MC, while IVTAM, used
C     by the tampering code, does not. TCC does not contain data for MC.
C     (MC is not included in the VVGAS vehicle category.)
C                                        / IVGAS = 1,2,3,4,5,6,7,8,9,10,11,13,14,15
C    Mapping of IVGAS to IVTAM vehicles <          | | | | | | | | |  |  |  |  |  |
C                                        \ IVTAM = 1,2,3,4,5,6,7,8,9,10,11,13, 0,14
      IVTAM = 0
      IVTAM = MIN0(14,IVGAS*VVGAS(IV))
      IF(VVGAS(IV) .EQ. 1) THEN
        EVAP_GM_MI(ICRNCASE)= EVAP_GM_MI(ICRNCASE)+TCC(IDX,IVTAM)*ZEV
      ENDIF
C
C     DAILY AVERAGE DIURNAL EMISSIONS :
      CALL DIRNLDAY(MY,ACTUAL_MY,ICY,IVGAS)
C     DETERMINE THE MULTIDAY DIURNAL RATIOS
      CALL DIURATIO(IV, MY)
C
C     EVAPORATIVE EMISSIONS HOURLY LOOP
      DO IH = 1, MAXIH
C
C       RESTING LOSSES:
C       need hourly resting losses for hot soak, diurnals, and running losses
C       (units conversions, ZEV, and oxygenated fuel adjustments will be
C        applied to resting losses at the end of the hourly loop)
C
        CALL RSTLOS(IH,MY,ACTUAL_MY,ICY,IVGAS)   !units of grams/hour
C
C       DIURNAL LOSSES:
        IF ( IH .LE. MAX_DIURNL_HR ) THEN  !no diurnals past midnight
C
C         find the hourly diurnal emissions
          DRNL = 0.0
C         DIRNL_HR produces hourly emissions with units of grams/hour for a
C         given IH and MY/IDX.
C
          DRNL=ZEV*DIRNL_HR(IH,MY,IVGAS)
          CALL FUEL2(IDIURNAL,DRNL)
C
          DRNL_GM_HR(IH) = DRNL
          DRNL_GM_MI(IH) = DRNL*
     *              CNVRT_U(GPH,GPM,JDX,IV,IH,IDIURNAL,INERR)
C
C         average daily diurnal losses in units of grams/day:
          EVAP_GM_MI(IDIURNAL)=EVAP_GM_MI(IDIURNAL)+
     *                         DRNL_GM_MI(IH)*HVMT(IH)
          EVAP_GM_HR(IDIURNAL)=EVAP_GM_HR(IDIURNAL)+
     *                         DRNL_GM_HR(IH)/24.0
C
        END IF
C
C       HOT SOAK and RUNNING LOSS emissions have 24 daily intervals.
C       However, a weighted 14th interval includes the emissions from
C       the 14th through 24th intervals. To produce a total of 24
C       hourly hot soak and running loss values the (weighted) 14th
C       interval is subdivided equally into the 14th through 24th hours
C       in HOTSOK and in RNLOSS
C       (see the last calculation in HOTSOK and in RNLOSS).
C
C       HOT SOAK LOSSES:
        HSK=0.0
        HSK=ZEV*HOTSOK(IH,MY,ACTUAL_MY,IVGAS,JDX)
        CALL FUEL2(IHOTSOAK,HSK)
C
C       HOT SOAK hourly numbers
        HS_GM_TRP(IH) =  HSK  !hot soak units are in gram/trip
        HS_GM_MI(IH)  =  HSK*
     *           CNVRT_U(GPT,GPM,JDX,IV,IH,IHOTSOAK,INERR)
        HS_GM_HR(IH)  =  HSK*
     *           CNVRT_U(GPT,GPH,JDX,IV,IH,IHOTSOAK,INERR)
C
C       HOT SOAK average daily numbers
        EVAP_GM_MI(IHOTSOAK) =EVAP_GM_MI(IHOTSOAK)+HS_GM_MI(IH)*HVMT(IH)
        EVAP_GM_TRP          =EVAP_GM_TRP         +HS_GM_TRP(IH)/24.0
        EVAP_GM_HR(IHOTSOAK) =EVAP_GM_HR(IHOTSOAK)+HS_GM_HR(IH)/24.0
C
C
C       initialize RUNNING LOSSES:
C
        DO IROAD = IFREEWAY,IALLROADS
          RUNLOSS(IROAD)=0.0
        END DO
C
        IF( VVGAS(IV).EQ.1 ) THEN               !no motorcycles
C
C         RUNNING LOSSES:
          CALL RNLOSSSUM(IH,IV,MY,ACTUAL_MY,JDX,RUNLOSS)
C
          DO IROAD = IFREEWAY,IALLROADS
            RUNLOSS(IROAD)=ZEV*RUNLOSS(IROAD)   !units of grams/mile for a particular hour of the day
C
C           RUNNING LOSS hourly numbers
            RN_GM_MI(IH,IROAD) = RUNLOSS(IROAD) !running loss units are in grams/mile
            RN_GM_HR(IH,IROAD) = RUNLOSS(IROAD)*
     *              CNVRT_U(GPM,GPH,JDX,IV,IH,IRUNLOSS,INERR)
C
          END DO
C
C         RUNNING LOSS average daily numbers
          EVAP_GM_MI(IRUNLOSS) =EVAP_GM_MI(IRUNLOSS) +
     *                                  RN_GM_MI(IH,IALLROADS)*HVMT(IH)
          EVAP_GM_HR(IRUNLOSS) =EVAP_GM_HR(IRUNLOSS) +
     *                                  RN_GM_HR(IH,IALLROADS)/24.0
C
C         REFUELING LOSSES:
C         the hourly array is built in REFUEL called by EFCALX
          RFL_GM_MI(IH)=ZEV*RLRATE(IDX,IVTAM,IH)
          RFL_GM_HR(IH)=ZEV*RLRATE(IDX,IVTAM,IH)*
     *               CNVRT_U(GPM,GPH,JDX,IV,IH,IREFUEL,INERR)
          RFL_GM_GA(IH)=ZEV*RLRATE(IDX,IVTAM,IH)*
     *               CNVRT_U(GPM,GPG,JDX,IV,IH,IREFUEL,INERR)
C
C         average daily numbers
          EVAP_GM_MI(IREFUEL)=EVAP_GM_MI(IREFUEL) + RFL_GM_MI(IH)
     *                          *HVMT(IH)
          EVAP_GM_HR(IREFUEL)=EVAP_GM_HR(IREFUEL) + RFL_GM_HR(IH)/24.0
          EVAP_GM_GA         =EVAP_GM_GA          + RFL_GM_GA(IH)/24.0
C
C
        ENDIF
C
C       there are no RUNNING or REFUELING losses for motorcycles
        IF( IV .EQ. VTMC ) THEN
C
          RFL_GM_MI(IH)      = 0.0
          RFL_GM_HR(IH)      = 0.0
          RFL_GM_GA(IH)      = 0.0
          DO IROAD = IFREEWAY,IALLROADS
C
C         RUNNING LOSS hourly numbers
            RN_GM_MI(IH,IROAD) = 0.0
            RN_GM_HR(IH,IROAD) = 0.0
          END DO
C
C         RUNNING & REFUELING LOSS average daily numbers
          EVAP_GM_MI(IRUNLOSS) =0.0
          EVAP_GM_HR(IRUNLOSS) =0.0
          EVAP_GM_MI(IREFUEL)  =0.0
          EVAP_GM_HR(IREFUEL)  =0.0
          EVAP_GM_GA           =0.0
C
        ENDIF
C
C       CRANK CASE LOSSES:
C       hourly numbers need to be converted
        CC_GM_MI(IH) = EVAP_GM_MI(ICRNCASE)*
     *       CNVRT_U(GPM_D,GPM,JDX,IV,IH,ICRNCASE,INERR) ! CNVRT_U=1 for this case
        CC_GM_HR(IH) = CC_GM_MI(IH)*
     *       CNVRT_U(GPM,GPH,JDX,IV,IH,ICRNCASE,INERR)
C
C       average daily numbers
        EVAP_GM_HR(ICRNCASE)=EVAP_GM_HR(ICRNCASE)+CC_GM_HR(IH)/24.0
C
C
C       FINISH RESTING LOSSES:
C       adjust for ZEV and NGV
        RSTLSS_EF = ZEV*RSTLSS_EF      !units of grams/hour
C
C       units conversions :
        RSTLS_GM_HR(IH) = RSTLSS_EF   !resting losses are in grams/hour
        RSTLS_GM_MI(IH) = RSTLSS_EF*
     *             CNVRT_U(GPH,GPM,JDX,IV,IH,IRESTLOSS,INERR)
C
C       daily averages :
        EVAP_GM_MI(IRESTLOSS)=EVAP_GM_MI(IRESTLOSS)+
     *                        RSTLS_GM_MI(IH)*HVMT(IH)
        EVAP_GM_HR(IRESTLOSS)=EVAP_GM_HR(IRESTLOSS)+
     *                        RSTLS_GM_HR(IH)/24.0
C
C
      END DO
C     End of hourly loop.
C
      IROAD=5   ! Restingloss,Refuel,Hotsoak,Diurnal and Crankcase are
                ! road independent.
C
C    Build DATABASE evaporative emissions.
C
ccs      CALL DBEVAP(FRESTLOSS,IV,JDX,RSTLS_GM_MI,RSTLS_GM_HR,IROAD,
ccs     *            ACTUAL_MY,INERR)
ccs      CALL DBEVAP(FREFUEL,IV,JDX,RFL_GM_MI,RFL_GM_HR,IROAD,ACTUAL_MY,
ccs     *            INERR)
ccs      CALL DBEVAP(FHOTSOAK,IV,JDX,HS_GM_MI,HS_GM_HR,IROAD,ACTUAL_MY,
ccs     *            INERR)
ccs      CALL DBEVAP(FDIURNAL,IV,JDX,DRNL_GM_MI,DRNL_GM_HR,IROAD,ACTUAL_MY,
ccs     *            INERR)
ccs      CALL DBEVAP(FCRNCASE,IV,JDX,CC_GM_MI,CC_GM_HR,IROAD,ACTUAL_MY,
ccs     *            INERR)
ccs      DO IROAD=IFREEWAY,IRAMP
ccs        CALL DBEVAP(FRUNLOSS,IV,JDX,RN_GM_MI(1,IROAD),RN_GM_HR(1,IROAD),
ccs     &                              IROAD,ACTUAL_MY,INERR)
ccs      END DO

ccs   Call custom SMOKE functions to create aggregated database output
      CALL SMKEVPOTHER(FHOTSOAK,IV,JDX,HS_GM_MI,IROAD)
      CALL SMKEVPOTHER(FDIURNAL,IV,JDX,DRNL_GM_MI,IROAD)
      CALL SMKEVPOTHER(FRESTLOSS,IV,JDX,RSTLS_GM_MI,IROAD)
      CALL SMKEVPOTHER(FCRNCASE,IV,JDX,CC_GM_MI,IROAD)
      CALL SMKEVPOTHER(FREFUEL,IV,JDX,RFL_GM_MI,IROAD)
      
      DO IROAD = IFREEWAY, IRAMP
          CALL SMKEVPRUN(FRUNLOSS,IV,JDX,RN_GM_MI(1,IROAD),IROAD)
      END DO

   99 RETURN
      END SUBROUTINE EVP_HRLY
