      SUBROUTINE ATEVP_HRLY(ICY,IVGAS,IDX,IV,JDX,MY,ACTUAL_MY) 
C 
C  ATEVAP_HRLY calls the subroutines and functions to which determine hydrocarbon 
C  evaporative emissions for hot soak, diurnal, crankcase, resting losses, and 
C  running losses.  Hourly refueling losses were determined in a call to REFUEL.FOR 
C  in ATCALX. Hourly emissions are determined from the calls to the various individual  
C  evaporative routines. 
C 
C  Called by:  HCCALX 
C 
C  Calls:  CRANKC, DB_EVAP, DIRNL_HR, DIRNLDAY, DIURATIO, FUEL2, HOTSOK, IEVPTR, 
C          PCTLEV, QUITER, RNLOSSSUM, RSTLOS, CNVRT_U 
C 
C  Changes:  (Last change first) 
C 
C  26 Jul 01 AIR Task 20: New Subroutine. Based on EVP_HRLY.FOR. 
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
C    /ATOX4  /       HS_GM_MI,DRNL_GM_HR,RN_GM_HR,RFL_GM_MI,RSTLS_GM_MI 
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
C  CC_GM_MI      R    crankcase emissions in grams per mile 
C  DRNL          R    hourly diurnal losses 
C  DRNL_GM_MI    R    diurnal losses in grams per mile 
C  HS_GM_MI      R    hot soak emissions in grams per mile 
C  HSK           R    hourly hot soak losses; HOTSOK + tampering hs offset, if any 
C  IALLROADS     I    roadway index (5 : all roads weighted together) 
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
C  RFL_GM_MI     R    refueling emissions in grams per mile 
C  RN_GM_MI      R    running loss emissions in grams per mile 
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
      INCLUDE 'ATOX4.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'DAILY_EVAP.I' 
      INCLUDE 'EVAPAR.I' 
      INCLUDE 'INJECT.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'NGV.I' 
      INCLUDE 'RESTLS1.I' 
      INCLUDE 'RLCOM3.I' 
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
C 
C     Declare external functions. 
C 
      INTEGER IEVPTR 
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
      INTEGER, PARAMETER :: IRESTLOSS=5 
      INTEGER, PARAMETER :: IREFUEL=6 
C                       three of the five roadway type indices 
      INTEGER, PARAMETER :: IFREEWAY=1 
      INTEGER, PARAMETER :: IRAMP=4 
      INTEGER, PARAMETER :: IALLROADS=5 
C 
C     Initialize the hourly evaporative emissions array. 
C 
      ATHEV = 0.0 
C 
C     Warn and exit evaporative calculations if the miles per day is zero. 
C 
      IF(BMYMPD(JDX,IV).LE.0.0) THEN 
        CALL QUITER(FLOAT(MY),IVGAS,58,0) 
        GOTO 99 
      ENDIF 
C 
C     Set IVTAM 
C 
      IVTAM = MIN0(14,IVGAS*VVGAS(IV)) 
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
C         DIRNL_HR produces hourly emissions with units of grams/hour for a 
C         given IH and MY/IDX. 
C 
          DRNL=ZEV*DIRNL_HR(IH,MY,IVGAS) 
CDFK      IF(IV.EQ.1)WRITE(99,91) DRNL 
CDFK  91  FORMAT(' D1',F10.4) 
          CALL FUEL2(IDIURNAL,DRNL) 
CDFK      IF(IV.EQ.1)WRITE(99,92) DRNL 
CDFK  92  FORMAT(' D2',F10.4) 
C 
C         Convert to g/mi 
          DRNL = DRNL * CNVRT_U(GPH,GPM,JDX,ICY,IV,IVTAM,IH,IDIURNAL) 
CDFK      IF(IV.EQ.1)WRITE(99,93) DRNL 
CDFK  93  FORMAT(' D3',F10.4) 
C 
C         Diurnals are independent of road type, so use IROAD=5 
          ATHEV(IH,5,2)=DRNL 
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
        HSK=ZEV*HOTSOK(IH,MY,ACTUAL_MY,IVGAS,JDX) !hot soak units are in gram/trip 
        CALL FUEL2(IHOTSOAK,HSK) 
C 
C       Convert HOT SOAK to grams/mi 
        HSK =  HSK * CNVRT_U(GPT,GPM,JDX,ICY,IV,IVTAM,IH,IHOTSOAK) 
C 
C       Hot Soaks are independent of road type, so use IROAD=5 
        ATHEV(IH,5,1)=HSK 
C 
C 
C       RUNNING LOSSES: 
C 
        IF( VVGAS(IV).EQ.1 ) THEN ! there are no RUNNING or REFUELING losses for motorcycles 
C 
C         RUNNING LOSSES: 
          CALL RNLOSSSUM(IH,IV,MY,ACTUAL_MY,JDX,RUNLOSS) 
C 
          DO IROAD = IFREEWAY,IALLROADS  ! Running losses ARE road type dependent 
            ATHEV(IH,IROAD,3)=ZEV*RUNLOSS(IROAD)   !units of grams/mile for a particular hour of the day 
          END DO 
C 
C         REFUELING LOSSES: 
C         The hourly array is built in REFUEL called by ATCALX 
C         Refueling losses are independent of road type, so use IROAD=5 
          ATHEV(IH,5,5) =ZEV*RLRATE(IDX,IVTAM,IH) 
C 
        ENDIF 
C 
C 
C       FINISH RESTING LOSSES: 
C       adjust for ZEV and NGV 
        RSTLSS_EF = ZEV*RSTLSS_EF      !units of grams/hour 
C 
C       units conversions to grams/mi: 
        RSTLSS_EF = RSTLSS_EF* 
     *             CNVRT_U(GPH,GPM,JDX,ICY,IV,IVTAM,IH,IRESTLOSS) 
C 
C       Reseting losses are independent of road type, so use IROAD=5 
        ATHEV(IH,5,4) = RSTLSS_EF 
C 
C 
      END DO 
C     End of hourly loop. 
C 
   99 RETURN 
      END SUBROUTINE ATEVP_HRLY 
