      SUBROUTINE HCCALX(ICY,INERR)
C
C  HCCALX calculates HC emission factors.
C
C  Called by EFCALX.
C
C  Calls DBEVAP,DBEXHRS,DIRNLINI,EVAPCALC, EVP_HRLY, EFCVT, FAIL2,
C        HRLOOP, IEVPTR, PCTLEV, and QUITER.
C
C  Changes:
C
C  03 Jun 02 AIR Task Bug Fix 390 When spreadsheet output is requested,
C            the Expand Evap arrays (G_MI_EVP) must be calculated.
C  05 Apr 02 AIR Task Bug Fix 388, Modification of constants via QUITER.
C  12 Mar 02 AIR Task Bug384: Compute All Veh evap components so that
C            diesel evap and All Veh can be output properly(Note that
C            MOBILE6 does not calculate evaporative emissions for
C            diesels. The code associated with this note is for output
C            only and should not effect the results.)
C  15 Nov 01 AIR Task 27: added INERR to DBEVAP due to PMFLECON update.
C  21 Dec 00 @EPA-bag beta24, added code for database aggregate output
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed
C            MYCODE.I
C  19 Oct 00 @EPA-djb Task W2, Make the DBEXHRS and DBEVAP calls conditional.
C  26 Sep 00 AIR Task 08: Included trap to prevent use of NO EVP SAWTOOTH command
C            whenever Evap I/M programs exist.
C  25 Sep 00 @EPA-elg changed .GT. 0.0 condition to .GE. 0.0 condition throughout
C            this subroutine.
C  25 Sep 00 @EPA-elg bug 257 changed error condition that prevented ZEV vehicles.
C  15 Sep 00 AIR Task 08: Added call to EVAPCALC to computer the evap weighting
C            factors.
C  22 aug 00 @EPA-RAG bugs151&152; Added call to DBEVAP so that age 0 emissions
C            for motorcycles and heavy duty vehicles would be written to output;
C            also, motorcycles with ages less than 11 will have emissions sent
C            to database output
C  21 Jul 00 @EPA-RAG running loss; removed call to FAIL.FOR and
C            determination of pass/fail fractions, PFRATE's, for
C            purge/pressure evap. tests
C  07 Jul 00 AIR Task 05: Removed unused ICY from NGV_EF call.
C  16 May 00 @EPA-ddj BUG87, Removed affect of STARFG=GPS.
C  27 Apr 00 EPA-elg  Added call to NGV EF subroutine
C  24 Apr 00 AIR Task 02: Included INERR in HRLOOP parameter list
C  17 Apr 00 @EPA-djb Removed dead code.
C  20 MAR 00 @EPA-elg           changes for HRLOOP
C  30 Jan 00 @DynTel-bg  1-005  Calls to new diurnal subroutines and functions
C  26 Jan 00 @DynTel-MLA 1-011  Added var. HEF, call to DBEXHRS.
C  19 Jan 00 @DynTel-bag 1-006  renamed evap_hrly.for to evp_hrly.for
C  29 Nov 99 @DynTel-bag 1-006  Added grams/mile index, GPM, for
C            start emissions. Rearranged evaporative emissions variables and
C            emissions sums.
C   6 Sep 99 @DynTel-HXQ 2-684 The code to change unit for start
C            emssion from g/mi to g/start is moved to HRLOOP.FOR
C   5 Aug 99 @DynTel-ddj 2-699  Removed initialization of BYTAM.
C  22 Jun 99 @DynTel-MLA 2-694  Removed 'INCLUDE LEVBLK.' The LEV
C            variables are no longer used here.
C  2 Jun 99 @DynTel-BG 2-000 Put ICY in parameter list of HRLOOP
C  26 May 99 @DynTel-HXQ 2-462 Delete code to call EFNFTP
C  19 May 99 @DynTel-HXQ 2-650 replace the BMYTPD with STARTSPERDAY
C  18 May 99 @DynTel-BG 2-650 Added START(=ISR=1), RUN(=ISR=2), and
C            USE_SOAK parameters; USE_SOAK is used with the SOAK_FLAG to
C            determine if the start emissions soak time corrections
C            are to be applied
C  18 May 99 @DynTel-HQ 2-698 Delete the code to calculate SCF
C  6  may 99 @DynTel-mla 2-701 Check start, running, and total emissions sums.
C            TOTEVAP includes running, resting, diurnal, hotsoak, crankcase,
C            and refueling losses.
C  15 Apr 99 @DynTel-HQ 2-698 add the code to calculate SCF
C  19 Apr 99 @DynTel-ddj 2-000 Cleaned up initialization of variables.
C  07 Apr 99 @Dyntel-BAG 2-000 Merged 2-694, 2-696, and 2-701 changes.
C  01 Apr 99 @Dyntel-DDJ 2-696 Added hourly loop.
C  11 Mar 99 @Dyntel -RJD 2-694 Modified the BEFWGT parameter list.
C  18 Mar 99 @DynTel-MLA 2-701 Keep separate totals for start and running
C            emissions, and keep a separate total for evaporative emissions.
C  3  Mar 99 @Dyntel -bag 2-664 Merged hourly evaporative calculations.
C            This included deletion of calls to RNLOSS, RSTLOS,
C            restructuring of evaporative emissions code, and moving
C            some of the partial diurnal code moved to CCEVRT.
C                     
C  11 Feb 99 @Dyntel -RJD 2-663 Added STARFG enabling the user to select
C            starts emission units to be grams/mile (STARFG=GPM) or
C            grams/start (STARFG=GPS).
C   6 JAN 99 @DynTel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust
C             engine start and running
C  04 Nov 1998 @DynTel-ZK 2-000 Explicit Typing
C  14 Sep 98 @DynTel-RJD  2-663
C  10 Aug 98 @Dyntel -RJD 2-663 Replaced call to BEF with call to BEFWGT for 
C            light duty gasoline vehicles. Moved calculation of BYBEF4 and 
C            tampering adjustment from BEF to HCCALX.
C  09 Sep 98 @Dyntel -RJD 0-632 The basic emission factor for all light duty
C            gasoline vehicles and trucks is now obtained from BEFWGT instead
C            of BEF.
C
C  Input on call:
C
C    Constants: 
C    MAXVEH from include file IVTYPE.I
C
C    argument list: ICY,INERR
C    common blocks:
C    /ATOX1 / SSFLAG
C    /CUMCOM/ CUMMIL
C    /CVNAMS/ ILDGT   (IVTYPE.I)
C    /EVAPACT/SAWFLAG
C    /FLAGS2/ EVIMFLAG
C    /FLAGS3/ RLFLAG
C    /FLAGS4/ HCFLAG
C    /FLAGS5/ SOAK_FLAG,SPD_FLAG
C    /MAXIMA/ MAXYRS
C    /MYRCAL/ MEVMYR,TF
C    /REGION/ IREJN
C    /UNIT_PARAM/ GPM
C    /VMXCOM/ REGMIX,VMTMIX
C    /VVNAMS/ VVDSL, VVGAS, VVHEAVY, VVLDGAS, VVSTARTS (IVTYPE.I)
C    /VVTYPS/ VTGASBUS, VTMC  (IVTYPE.I)
C
C  Output on return:
C
C    argument list: INERR
C    common blocks:
C    /BYMYC2/ BYBEF4,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC5/ BYIMRE,REDSUM
C    /CVNAMS/ ILDGT   (IVTYPE.I)
C    /EVAPGR/ GREVP,IGREVP
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL,
C             IEFFTP,IEFEXH,IEFEVAP,IEFREFL,IEFRUNL,
C    /RESUL2/RLGGAL,RSTGPH,G_MI_EVP,IRLGGAL,IRSTGPH,IGMIEVP
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VDRNL,VHTSK,VCRNK
C    /RESUL4/ RUNEF,STARTEF,TOTEVAP,VTEVAP
C
C  Local array subscripts
C
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COMPCC     R    composite evaporative & crankcase emission factor for IDX
C  COMPEF     R    composite emission factor for IDX: summed across IDX to get
C                  EFFTP for the given IP & IV
C  EXHHC      R    exhaust HC component of total HC, for given IV x IDX case
C  GASCAP     R    refueling HC loss factor for given IDX and IV
C  ICY        I    calendar year
C  IEVP       I    evaporative emissions index : 1(=IHOTSOAK) -  hotsoak
C                                                2(=IDIURNAL) -  diurnal
C                                                3(=IRUNLOSS) -  running loss
C                                                4(=ICRNCASE) -  crankcase
C                                                5(=IRESTLOSS)-  resting loss
C                                                6(=IREFUEL)  -  refuelng loss
C  ICRNCASE   I    evaporative emissions index (see above)
C  IDIURNAL   I    evaporative emissions index (see above)
C  IHOTSOAK   I    evaporative emissions index (see above)
C  IREFUEL    I    evaporative emissions index (see above)
C  IRESTLOSS  I    evaporative emissions index (see above)
C  IRUNLOSS   I    evaporative emissions index (see above)
C  INERR      I    error count
C  IP         I    pollutant index (HC:1,CO:2,NOx;3)
C  IV         I    vehicle index
C  IVGAS      I    see (1,999)'s subscript dictionary
C  IVTAM      I    tampered vehicle index
C  MDX        I    model year index into road fuel economy arrays
C  RUNMYR     R    running emissions for one model year
C  STARTMYR   R    start emissions for one model year
C  USE_SOAK   I    (=1)used to check if hourly soak distributions are to
C                  be used
C
C  Notes:
C
C  HCCALX was modified for MOBILE4 In-House Version 02 to save by model year
C  computations for the BYMY Tables.
C  HCCALX was modified for MOBILE4.1 to compute grams/mile hot soak,
C  diurnal, and partial diurnal values.  Also, oxygenated and cert.
C  fuels support from MOBILE4 was removed in MOBILE4.1.
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor 
C  calculation function.
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4 
C  calculation to remove improperly making a SALHCF adjustment to by-model
C  Sep-19-1994 @ CSC-bsg Request 462 Made correction in HC BYMY FER
C  calculation, to properly account for tampering in the result printed.
C  Sep-21-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year
C  removed IY from the call to CCEVRT
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C  29 Oc. 1997 Dyntel@MMS 0-632a  deleted condition for calling PCTLEV
C  26 Feb 1998 Dyntel@RJD (after gjr) 2-204 Removed IY from parameter 
C              list and call to function BEF. 
C  03 Mar 98 @DynTel-RJD 2-204  Changed the name of the array storing the 
C             Model Year Registration fractions from JANMYR to MEVMYR to 
C             reflect the fact that January is no longer the only month of
C             evaluation that the Mobile model can handle.
C  11 Mar 98 @DynTel-RJD 2-622  Added calls to EFCVT so that emission
C             factors for composite LDGT vehicles are not calculated here.
C             VMLDGT has been removed from the parameter list.
C  23 Apr 98 @DynTel-RJD 2-654 Added ISR loops to distinguish between  
C            starts and running emissions. Added an extra dimension to
C            the BYBEF4 array so that starts and running emission can be
C            stored simultaneously.
C  27 Apr 98 @DynTel-BAG 2-661  Removed the travel fraction factor from BYFER
C  30 Jun 98 @DynTel-BAG 2-000  Removed I/M credit calculations, BYIMCR, BYIMRE,
C            and REDSUM             
C  8 July 98 Dyntel@MMS 2-666  used VVSTARTS(IV) in order to set Maxisr whenever
C            it was applicable and deleted  VTMC in certain positions , because 
C            VVSTARTS would include that.
C 23 Jul 98  @DynTel-MLA 2-000  Changed LDGT names to conform to
C            certification naming conventions.
C
      USE DATABASE, ONLY : DBFLAG, DBAGGR
C
      IMPLICIT NONE
      INCLUDE 'ATOX1.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC3.I'
      INCLUDE 'CUMCOM.I'
      INCLUDE 'DAILY_EVAP.I'
      INCLUDE 'EVAPACT.I'
      INCLUDE 'EVAPGR.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'FLAGS5.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'RESUL4.I'
      INCLUDE 'STDIST.I'
      INCLUDE 'UNIT_PARAM.I'
      INCLUDE 'VMXCOM.I'
      INCLUDE 'WEEKLY.I'
C
C     Declare external functions.
C
      REAL EFCVT
C
C     Declare parameter list.
C
      INTEGER, INTENT(IN)    :: ICY   ! calender year
      INTEGER, INTENT(INOUT) :: INERR ! error count
C
C     Declare local variables.
C
      REAL    COMPEF
      REAL    EXHHC
      REAL    HEF(24,5)
      REAL    NO_REFUEL
      REAL    RUNMYR
      REAL    STARTMYR
      REAL    VMTAGE
      INTEGER IEVP
      INTEGER IFAC
      INTEGER IHOUR
      INTEGER IP
      INTEGER IROAD
      INTEGER IV
      INTEGER IVGAS
      INTEGER IVTAM
      INTEGER IDX
      INTEGER JDX
      INTEGER ACTUAL_MY
      INTEGER MY
C
C
      INTEGER, SAVE :: EXPND_HC=2
      INTEGER, SAVE :: HC=1
      INTEGER, SAVE :: USE_SOAK=1
      INTEGER, SAVE :: IHOTSOAK=1,IDIURNAL=2,IRUNLOSS=3
      INTEGER, SAVE :: IRESTLOSS=5,ICRNCASE=4,IREFUEL=6
C                       two of the five roadway type indices
      INTEGER, PARAMETER :: IFREEWAY=1
      INTEGER, PARAMETER :: IRAMP=4
C                      evap. emission indices for database output
      INTEGER, SAVE :: FHOTSOAK=3
      INTEGER, SAVE :: FDIURNAL=4
      INTEGER, SAVE :: FRESTLOSS=5
      INTEGER, SAVE :: FRUNLOSS=6
      INTEGER, SAVE :: FCRNCASE=7
      INTEGER, SAVE :: FREFUEL=8
C                      for database output of vehicles, in particular motorcycles,
C                      in years without (mileage) data
      REAL,    SAVE :: NO_MILES_NO_EMISSIONS(24)=0.0
C
C     HEF is used to pass the database output values to DBEXHRS. If
C     DBEXHRS is called from here, the emission factors will all be
C     set to zero.
C
      DO IFAC = 1,5
        DO IHOUR = 1,24
          HEF(IHOUR,IFAC) = 0.0
        END DO
      END DO
C
C  Initialize
C
      NO_REFUEL   =0.0
      EFEXH(ILDGT)=0.0
      EFEVAP(ILDGT)=0.0
      REFL_GM_MI(ILDGT)=0.0
      RLGGAL(ILDGT)=0.0
      RSTGPH(ILDGT)=0.0
      EFRUNL(ILDGT)=0.0
      TOTEVAP(ILDGT)=0.0
      VTEVAP=0.0
      VEXH=0.0
      VEVAP=0.0
      VDRNL=0.0
      VHTSK=0.0
      VCRNK=0.0
      VLOSS=0.0
      VRUNLS=0.0
      VRSTLS=0.0
      DO IEVP=1,6
        G_MI_EVP(IEVP,ILDGT)=0.0 ! ILDGT = 29
        GREVP(IEVP,ILDGT)=0.0
      END DO
C
      IP=1
      EFFTP(IP,ILDGT)=0.0
      VFTP(IP)=0.0
      VRUN(IP)=0.0
      VSTART(IP)=0.0
C
C     Initialize gas vehicle counter, IVGAS, used for evaporative
C     emissions calculations
C
      IVGAS = 0
C
C     The following call to soak function, SOAKFCN, determines an array, SOAKFUN,
C     which are used to determine hourly soak time correction factors for use
C     in HRLOOP to calculate hourly emission factors. It is dependent only on
C     soak time, pollutant, and type of catalyst. So, calculate it for each
C     pollutant out of vehicle, start, hour, and year loops. (See parameters in
C     Table 4 of M6.STE.003)
C
      IF( SOAK_FLAG .EQ. USE_SOAK ) CALL SOAKFCN(HC)
C
C     Initialize the diurnal fractions
      CALL DIRNLINI()
C
C     To avoid problems with EVAPCALC call below, trap case where user has entered
C     Evap I/M programs but has also specified NO EVP SAWTOOTH. Reset SAWFLAG to 1
C     (sawtooth enabled) whenever this happens.
C
      IF(EVIMFLAG.GT.1 .AND. SAWFLAG.EQ.0) THEN
        CALL QUITER(0.,0,563,INERR)
        SAWFLAG = 1
      END IF
C
C     For each vehicle class, calculate exhaust, evap and partial diurnals,
C     and running loss values
C     Begin vehicle loop:
C
      DO 85 IV=1,MAXVEH
C
C       evaporative emissions vehicle index
C
        IVGAS = IVGAS + VVGASMC(IV)
C
        IVTAM = IV
        IF (IV.EQ.VTGASBUS) THEN
          IVTAM = 14
        ENDIF
C
        EFFTP(IP,IV)=0.0
        EFEXH(IV)=0.0
        EFEVAP(IV)=0.0
        REFL_GM_MI(IV)=0.0
        RLGGAL(IV)=0.0
        RSTGPH(IV)=0.0
        EFRUNL(IV)=0.0
        RUNEF(IP,IV) = 0.0
        STARTEF(IP,IV) = 0.0
        TOTEVAP(IV) = 0.0
        DO IEVP=1,6
          GREVP(IEVP,IV)=0.0
          G_MI_EVP(IEVP,IV)=0.0
        END DO
C
C       REGMIX is not affected by scenario input, only by
C       BLOCK DATA / one-time data values.  BYMY save areas
C       can be initialized after this check, since it would
C       be technically in error to have a VMTMIX(i) > 0.0
C       when REGMIX(i) = 0.0, so the situation described in
C       EFCALX for CO and NOX can not happen in HCCALX for HC.
C
        MY=0
C
C       calculate the speed correction factor for 2.5-65 m/h
C
        IF(VVLDGAS(IV).NE.1.AND.IV.NE.VTMC.AND.SPD_FLAG.EQ.1)
     &                       CALL SCFTAB(IV,IP)
C
C       Set up the Evap weighting factors, without and with I/M,
C       and without and with OBD, for all ages. Only call this
C       routine for all gasoline powered vehicles.
C
        IF(VVGASMC(IV).EQ.1) CALL EVAPCALC(ICY,IV)
C
C       Begin yearly loop:
C
        DO 70 IDX=1,MAXYRS
          JDX=(MAXYRS+1)-IDX
C
          BYFER(IP,JDX,IV)=0.0
          BYEVAP(JDX,IV)=0.0
          BYRSTL(JDX,IV)=0.0
C
C         No by model year tampering, running, or refueling losses
C         for diesel or motorcycles. These apply only to light duty
C         and heavy duty gas vehicles and gas buses. Initialize the
C         arrays which hold these numbers.
C
          IF(VVGAS(IV).EQ.1) THEN
            BYRUNL(JDX,IVTAM)=0.0
            BYREFL(JDX,IVTAM)=0.0
          END IF
C
          ACTUAL_MY=ICY+IDX-MAXYRS
          MY=ACTUAL_MY
C
          IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN
            MY = MIN(CAAMY,ACTUAL_MY)
          ENDIF
C
C         Skip the calculation for a model year in which a vehicle
C         has zero mileage, but produce database output if requested.
C
          IF((TF(IDX,IV).LE.0.0.AND.MEVMYR(JDX,IV).EQ.0.0) .OR.
     &       (VMTMIX(IV).EQ.0.0)) THEN
            IF(DBFLAG.EQ.2) THEN
ccs            CALL DBEXHRS(HEF,IP,IV,ACTUAL_MY,JDX,INERR)

ccs   Call custom SMOKE routine to aggregate exhaust emmisions
            CALL SMKEXOUT(HEF,IP,IV,JDX)

            IROAD=5   ! Restingloss,Refuel,Hotsoak,Diurnal and Crankcase are
                      ! road independent.
ccs            CALL DBEVAP(FRESTLOSS,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                  NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            CALL DBEVAP(FREFUEL,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                  NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            CALL DBEVAP(FHOTSOAK,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                  NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            CALL DBEVAP(FDIURNAL,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                  NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            CALL DBEVAP(FCRNCASE,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                  NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            DO IROAD=IFREEWAY,IRAMP
ccs              CALL DBEVAP(FRUNLOSS,IV,JDX,NO_MILES_NO_EMISSIONS,
ccs     *                    NO_MILES_NO_EMISSIONS,IROAD,ACTUAL_MY,INERR)
ccs            END DO

ccs   Call custom SMOKE routines to aggregate evaporative emissions
        CALL SMKEVPOTHER(FHOTSOAK,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        CALL SMKEVPOTHER(FDIURNAL,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        CALL SMKEVPOTHER(FRESTLOSS,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        CALL SMKEVPOTHER(FCRNCASE,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        CALL SMKEVPOTHER(FREFUEL,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        DO IROAD=IFREEWAY,IRAMP
          CALL SMKEVPRUN(FRUNLOSS,IV,JDX,NO_MILES_NO_EMISSIONS,IROAD)
        END DO
          
            ENDIF !DBFLAG conditional
            CYCLE
          ENDIF
C
          EXHHC  = 0.0
C
          COMPEF=0.0
          RUNMYR = 0.0
          STARTMYR = 0.0
C
          VMTAGE = CUMMIL(JDX,IV) / 10000.
          CALL NGV_EF(IP,IV,MY,IDX,VMTAGE)
C
          CALL HRLOOP(RUNMYR,STARTMYR,IP,IV,ICY,MY,ACTUAL_MY,IDX,INERR)
C
          EXHHC = RUNMYR + STARTMYR
C
          COMPEF = COMPEF + EXHHC
          EXHHC = EXHHC*TF(IDX,IV)
          EFEXH(IV) = EFEXH(IV) + EXHHC
          RUNEF(IP,IV) = RUNEF(IP,IV) + RUNMYR * TF(IDX,IV)
          STARTEF(IP,IV) = STARTEF(IP,IV) + STARTMYR * TF(IDX,IV)
C
C         Diesels (IV=5-7) do not have evaporative, on-board or RVP
C         components, so branch accordingly.
C
          IF(VVDSL(IV).EQ.1) GOTO 65
C
C         Now determine the composite evaporative and crankcase
C         factors for current MY, including any onboard effect.
C         Do the hourly computation for gas vehicles.
C
C         Call to EVP_HRLY computes hourly gram evap components.
C
C         Daily resting, running, crankcase, diurnal, hot soak, are determined by
C         adding respective hourly emissions in evp_hrly.for for gas fueled vehicles.
C
          CALL EVP_HRLY(ICY,IVGAS,IDX,IV,JDX,MY,ACTUAL_MY,INERR)
C
C         EXPANDED EVAPORATIVE OUTPUT(HCFLAG=2, 3, or 4)
C         For gas vehicles use the results from EVP_HRLY to determine composite
C         evaporative sums, GREVP(EVAP_TYPE,IV), over 25 years in grams/mile
          IF((HCFLAG.GE.EXPND_HC)
     *       .OR. DBAGGR .OR. SSFLAG) THEN  ! output in grams/mile and
                                            ! no refueling in composite sums.
                                            ! The database aggregate and
                                            ! spreadsheet output also
                                                   ! needs these arrays
C
            G_MI_EVP(IDIURNAL,IV)=G_MI_EVP(IDIURNAL,IV)
     *                        + EVAP_GM_MI(IDIURNAL)*TF(IDX,IV)
            G_MI_EVP(IHOTSOAK,IV)=G_MI_EVP(IHOTSOAK,IV)
     *                        + EVAP_GM_MI(IHOTSOAK)*TF(IDX,IV)
            G_MI_EVP(IRUNLOSS,IV)=G_MI_EVP(IRUNLOSS,IV)
     *                        + EVAP_GM_MI(IRUNLOSS)*TF(IDX,IV)
            G_MI_EVP(ICRNCASE,IV)=G_MI_EVP(ICRNCASE,IV)
     *                        + EVAP_GM_MI(ICRNCASE)*TF(IDX,IV)
            G_MI_EVP(IREFUEL,IV)= G_MI_EVP(IREFUEL,IV)
     *                        + EVAP_GM_MI(IREFUEL)*TF(IDX,IV)
            G_MI_EVP(IRESTLOSS,IV)=G_MI_EVP(IRESTLOSS,IV)
     *                        + EVAP_GM_MI(IRESTLOSS)*TF(IDX,IV)
C
          ENDIF
C
C         ALTERNATE EVAPORATIVE UNITS(HCFLAG=3)
C         composite evaporative sum over 25 years in grams/hour
C         hot soak is in grams/trip, crankcase is in grams/mile
C         running losses in grams/gal. are done below
          IF(HCFLAG .EQ. 3 ) THEN
C
            GREVP(IDIURNAL,IV)=GREVP(IDIURNAL,IV)
     *                   + EVAP_GM_HR(IDIURNAL)*MEVMYR(JDX,IV)
            GREVP(IHOTSOAK,IV)=GREVP(IHOTSOAK,IV)
     *                   + EVAP_GM_TRP*MEVMYR(JDX,IV)
            GREVP(IRUNLOSS,IV)=GREVP(IRUNLOSS,IV)
     *                   + EVAP_GM_HR(IRUNLOSS)*MEVMYR(JDX,IV)
            GREVP(ICRNCASE,IV)=GREVP(ICRNCASE,IV)
     *                   + EVAP_GM_HR(ICRNCASE)*MEVMYR(JDX,IV)
            GREVP(IREFUEL,IV)=GREVP(IREFUEL,IV)
     *                   + EVAP_GM_HR(IREFUEL)*MEVMYR(JDX,IV)
            RSTGPH(IV) = RSTGPH(IV)
     *                   + EVAP_GM_HR(IRESTLOSS)*MEVMYR(JDX,IV)
C
          ENDIF
C
C
C         by model year resting losses is in grams/mile
          BYRSTL(JDX,IV) = EVAP_GM_MI(IRESTLOSS)

C
C         by model year evaporative sum (diurnals,hot soak, and crank case) in grams/mile
          BYEVAP(JDX,IV)= EVAP_GM_MI(IDIURNAL)
     *                              +  EVAP_GM_MI(IHOTSOAK)
     *                              +  EVAP_GM_MI(ICRNCASE)
C
C         add diurnal, hot soak, crank case, and resting loss to
C         HC exhaust sum
          COMPEF=COMPEF + EVAP_GM_MI(IDIURNAL)
     *                  + EVAP_GM_MI(IHOTSOAK)
     *                  + EVAP_GM_MI(ICRNCASE)
     *                  + EVAP_GM_MI(IRESTLOSS)
C
C         composite by vehicle type sum over 25 years
          EFEVAP(IV)=EFEVAP(IV)+ ( EVAP_GM_MI(IDIURNAL)
     *                         +   EVAP_GM_MI(IHOTSOAK)
     *                         +   EVAP_GM_MI(ICRNCASE) )*TF(IDX,IV)
C
C         RUNNING LOSS additions (always output in grams/mile)
C         Motorcycles do not have running loss emissions.
          IF(IV.NE.VTMC) THEN
C
C           by model year output of running losses is in grams/mile
            BYRUNL(JDX,IVTAM) = EVAP_GM_MI(IRUNLOSS)
C
C           composite gram/mile output of running losses
            EFRUNL(IV)=EFRUNL(IV) + EVAP_GM_MI(IRUNLOSS)*TF(IDX,IV)
C
C           add in running loss emissions to evap. plus exhaust sum
C           for gas fueled vehicles.
            COMPEF=COMPEF + EVAP_GM_MI(IRUNLOSS)
C
          ENDIF
C
C
C         REFUELING LOSS additions
C         Refueling losses (GASCAP) only occur for gasoline vehicles (IV = 1-4)
C         under one of the nonzero refueling loss options (RLFLAG = 1-4).
C         For the HCFLAG=3 case, do not include refueling losses in composite
C         emission factor or on refueling loss g/mi line. They will be
C         output in g/gal section (see below).
C
          NO_REFUEL = 0.
          IF(VVGAS(IV).EQ.1.AND.RLFLAG.NE.5) THEN
C
C           REFUELING LOSSES are included in the total emissions sum when
C           HCFLAG is NOT equal to 4.
C           (7 Jan 93 - MLA)
            IF(HCFLAG.NE.4) THEN
              COMPEF = COMPEF + EVAP_GM_MI(IREFUEL)
              NO_REFUEL = 1.
            ENDIF
C
C           by model year output of refueling losses is in grams/mile
            BYREFL(JDX,IVTAM) = EVAP_GM_MI(IREFUEL)
C
C           SET EFRREL unconditionally now, so the numbers will appear in the
C           output. Jan 93 - MLA
C
            REFL_GM_MI(IV) = REFL_GM_MI(IV) +
     *                       EVAP_GM_MI(IREFUEL) * TF(IDX,IV)
C
C           Sum the refueling loss in g/gal for the alternate units output
            IF(HCFLAG.EQ.3) THEN
C
              RLGGAL(IV) = RLGGAL(IV) + EVAP_GM_GA * TF(IDX,IV)
C
            ENDIF
C
          ENDIF
C
C         evaporative sum (always in grams/mile)
C
          TOTEVAP(IV) = TOTEVAP(IV)
     *                   +(   EVAP_GM_MI(IDIURNAL)
     *                      + EVAP_GM_MI(IHOTSOAK)
     *                      + EVAP_GM_MI(IRUNLOSS)
     *                      + EVAP_GM_MI(IRESTLOSS)
     *                      + EVAP_GM_MI(ICRNCASE)
     *                      + EVAP_GM_MI(IREFUEL)*NO_REFUEL )*TF(IDX,IV)
C
C         Update FTP composite HC.
C         (No longer multiplied by the travel fraction. Apr 27 1998)
C
   65     BYFER(IP,JDX,IV) = COMPEF
C
C         FTP emissions in units of grams/mile use travel fraction
C         weighting for composite emissions sum output.
C
          EFFTP(IP,IV) = EFFTP(IP,IV) + (COMPEF*TF(IDX,IV))
C
C     Abort whenever any errors occur.
C
          IF(INERR.GT.0) GOTO 99
C
   70   CONTINUE   !End of yearly loop
C
        IF(EFFTP(IP,IV).LT.0.0.AND.VMTMIX(IV).GT.0.0)
     *  CALL QUITER(0.,IV,68,INERR)
C
C       Weight by VMTMIX for mileage components and REGMIX for gram components.
C       sum over the diferent vehicle classes weighted by mileage fractions to
C       get the total emissions for a particular emission type
        VFTP(IP) = VFTP(IP) + EFFTP(IP,IV)*VMTMIX(IV)
        VEXH     = VEXH     + EFEXH(IV)*VMTMIX(IV)
        VRUN(IP) = VRUN(IP) + RUNEF(IP,IV)*VMTMIX(IV)
        VSTART(IP) = VSTART(IP) + STARTEF(IP,IV)*VMTMIX(IV)
C
C       Evaporative Sums: Each vehicle type has the evaporative component
C       initialized to zero. Diesel vehicles do not contribute to these 
C       25 year vehicle emissions sums.
C
        VEVAP  = VEVAP  + EFEVAP(IV)*VMTMIX(IV)
        VDRNL  = VDRNL  + G_MI_EVP(IDIURNAL,IV)*VMTMIX(IV)
        VHTSK  = VHTSK  + G_MI_EVP(IHOTSOAK,IV)*VMTMIX(IV)
        VCRNK  = VCRNK  + G_MI_EVP(ICRNCASE,IV)*VMTMIX(IV)
        VLOSS  = VLOSS  + REFL_GM_MI(IV)*VMTMIX(IV)
        VRUNLS = VRUNLS + EFRUNL(IV)*VMTMIX(IV)
        VRSTLS = VRSTLS + G_MI_EVP(IRESTLOSS,IV)*VMTMIX(IV)
        VTEVAP = VTEVAP + TOTEVAP(IV)*VMTMIX(IV)
C
   85 CONTINUE
C     end of vehicle loop
C
C     Use registration fractions(grams/hour output) or
C     mileages fractions(grams/mile output) to get weighted
C     emission factor sums for composite vehicle types.
C
C     EXHAUST EMISSIONS:
C      FTP
      IF(EFFTP(IP,2).GE.0.0.AND.EFFTP(IP,3).GE.0.0.AND.
     * EFFTP(IP,4).GE.0.0.AND.EFFTP(IP,5).GE.0.0)
     * EFFTP(IP,ILDGT) = EFCVT(ILDGT,IEFFTP,IP,INERR)
C
C      START
      IF(STARTEF(IP,2).GE.0.0.AND.STARTEF(IP,3).GE.0.0.AND.
     * STARTEF(IP,4).GE.0.0.AND.STARTEF(IP,5).GE.0.0)
     * STARTEF(IP,ILDGT) = EFCVT(ILDGT,ISTARTEF,IP,INERR)
C
C      RUNNING
      IF(RUNEF(IP,2).GE.0.0.AND.RUNEF(IP,3).GE.0.0.AND.
     * RUNEF(IP,4).GE.0.0.AND.RUNEF(IP,5).GE.0.0)
     * RUNEF(IP,ILDGT) =  EFCVT(ILDGT,IRUNEF,IP,INERR)
C
C      TOTAL EXHAUST
      IF(EFEXH(2).GE.0.0.AND.EFEXH(3).GE.0.0.AND.EFEXH(4).GE.0.0.AND.
     * EFEXH(5).GE.0.0) EFEXH(ILDGT)= EFCVT(ILDGT,IEFEXH,IP,INERR)
C
C     diurnal, hot soak, resting loss, and crank case sum
      IF(EFEVAP(2).GE.0.0.AND.EFEVAP(3).GE.0.0.AND.EFEVAP(4).GE.0.0
     * .AND.EFEVAP(5).GE.0.0)
     * EFEVAP(ILDGT) = EFCVT(ILDGT,IEFEVAP,IP,INERR)
C
C     refueling (grams/mile)
      IF(REFL_GM_MI(2).GE.0.0.AND.REFL_GM_MI(3).GE.0.0
     * .AND.REFL_GM_MI(4).GE.0.0
     * .AND.REFL_GM_MI(5).GE.0.0)
     * REFL_GM_MI(ILDGT) = EFCVT(ILDGT,IEFREFL,IP,INERR)
C
C     refueling (grams/gal)
      IF(RLGGAL(2).GE.0.0.AND.RLGGAL(3).GE.0.0.AND.RLGGAL(4).GE.0.0
     * .AND.RLGGAL(5).GE.0.0)
     * RLGGAL(ILDGT) = EFCVT(ILDGT,IRLGGAL,IP,INERR)
C
C     hot soak(grams/mile)
      IF(G_MI_EVP(IHOTSOAK,2).GE.0.0.AND.G_MI_EVP(IHOTSOAK,3).GE.0.0
     *.AND.G_MI_EVP(IHOTSOAK,4).GE.0.0
     *.AND.G_MI_EVP(IHOTSOAK,5).GE.0.0)
     *G_MI_EVP(IHOTSOAK,ILDGT) = EFCVT(ILDGT,IGMIEVP,IHOTSOAK,INERR)
C
C     diurnal loss(grams/mile)
      IF(G_MI_EVP(IDIURNAL,2).GE.0.0.AND.G_MI_EVP(IDIURNAL,3).GE.0.0
     *.AND.G_MI_EVP(IDIURNAL,4).GE.0.0
     *.AND.G_MI_EVP(IDIURNAL,5).GE.0.0)
     *G_MI_EVP(IDIURNAL,ILDGT) = EFCVT(ILDGT,IGMIEVP,IDIURNAL,INERR)
C
C     resting loss(grams/mile)
      IF(G_MI_EVP(IRESTLOSS,2).GE.0.0.AND.G_MI_EVP(IRESTLOSS,3).GE.0.0
     *.AND.G_MI_EVP(IRESTLOSS,4).GE.0.0
     *.AND.G_MI_EVP(IRESTLOSS,5).GE.0.0)
     *G_MI_EVP(IRESTLOSS,ILDGT) = EFCVT(ILDGT,IGMIEVP,IRESTLOSS,INERR)
C
C     crank case(grams/mile)
      IF(G_MI_EVP(ICRNCASE,2).GE.0.0.AND.G_MI_EVP(ICRNCASE,3).GE.0.0
     *.AND.G_MI_EVP(ICRNCASE,4).GE.0.0
     *.AND.G_MI_EVP(ICRNCASE,5).GE.0.0)
     *G_MI_EVP(ICRNCASE,ILDGT) = EFCVT(ILDGT,IGMIEVP,ICRNCASE,INERR)
C
C     resting loss(grams/hr)
      IF(RSTGPH(2).GE.0.0.AND.RSTGPH(3).GE.0.0.AND.RSTGPH(4).GE.0.0
     * .AND.RSTGPH(5).GE.0.0)
     * RSTGPH(ILDGT) = EFCVT(ILDGT,IRSTGPH,IP,INERR)
C
C     running loss
      IF(EFRUNL(2).GE.0.0.AND.EFRUNL(3).GE.0.0.AND.EFRUNL(4).GE.0.0
     * .AND.EFRUNL(5).GE.0.0)
     * EFRUNL(ILDGT) = EFCVT(ILDGT,IEFRUNL,IP,INERR)
C
C     total evaporative emissions composite light duty trucks
      IF(TOTEVAP(2).GE.0.0.AND.TOTEVAP(3).GE.0.0.AND.
     *   TOTEVAP(4).GE.0.0.AND.TOTEVAP(5).GE.0.0)
     *   TOTEVAP(ILDGT) = (TOTEVAP(2)*REGMIX(2)+
     *                   TOTEVAP(3)*REGMIX(3)+
     *                   TOTEVAP(4)*REGMIX(4)+
     *                   TOTEVAP(5)*REGMIX(5)) /
     *       ( REGMIX(2) + REGMIX(3) + REGMIX(4) + REGMIX(5))
C
C     For alternate units output:
      IF(HCFLAG .EQ. 3) THEN
C     weight individual light duty truck type evaporative emissions
C     by registration fractions for composite light duty trucks
C     units depend on HCFLAG (see above)
C
        DO IEVP = 1, 3
          IF(GREVP(IEVP,2).GE.0.0.AND.GREVP(IEVP,3).GE.0.0.AND.
     *      GREVP(IEVP,4).GE.0.0.AND.GREVP(IEVP,5).GE.0.0)
     *      GREVP(IEVP,ILDGT) = (GREVP(IEVP,2)*REGMIX(2)+
     *                     GREVP(IEVP,3)*REGMIX(3)+
     *                     GREVP(IEVP,4)*REGMIX(4)+
     *                     GREVP(IEVP,5)*REGMIX(5)) /
     *         ( REGMIX(2) + REGMIX(3) + REGMIX(4) + REGMIX(5) )
        END DO
C
C
        IF(GREVP(4,2).GE.0.0.AND.GREVP(4,3).GE.0.0.AND.
     *    GREVP(4,4).GE.0.0.AND.GREVP(4,5).GE.0.0)
     *    GREVP(4,ILDGT) = EFCVT(ILDGT,IGREVP,4,INERR)
      ENDIF
C
  99  RETURN
      END
