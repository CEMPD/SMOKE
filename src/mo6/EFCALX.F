      SUBROUTINE EFCALX(ICY,INERR)
C
C  EFCALX calculates composite emission factors.
C
C  Called by PRCSCN.
C
C  Calls HRLOOP, BIGCFX, EFCVT, GETCUM, GSFIDX, HCCALX, QUITER, 
C        REFUEL, RVPEXH, OBD_HI_FRCN and TAMPER.
C
C  Changes: (Last change first)
C
C  05 Apr 02: AIR Task Bug Fix 388, Modification of constants via QUITER.
C  15 Nov 01 AIR Task 27: Added INERR support due to PMFLECON update.
C  27 Nov 00 @EPA-elg Added No Clean Air Act condition for tampering.
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed
C            MYCODE.I
C  06 Oct 00 @EPA-BG  Bug 155, Removed IDLFLG/FLAGS4.I 
C  25 Sep 00 @EPA-elg bug257 changed error condition that prevented ZEV vehicles. 
C  07 Jul 00 AIR Task 05: Removed unused ICY from NGV_EF call. 
C            Removed unused KKSTD from TEX1 call. 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  02 Jun 00 AIR Task 03: Removed call to BIGCFX. HRLOOP now calls BIGCFX. 
C   1 Jun 00 @EPA-ddj E11, Streamlined calculation of MY. 
C  16 May 00 @EPA-ddj BUG87, Removed STARFG conditionals. 
C  27 Apr 00 @EPA-elg  Added Call to Tampering routine TEX1. 
C  27 Apr 00 @EPA-elg  Calls NGV subroutine. 
C  24 Apr 00 AIR Task 02: Fixed call to REFUEL and HRLOOP to include INERR 
C  10 Apr 00 AIR Task 02: Replace PRTFLG with DBSELPOL 
C  17 Apr 00 @EPA-djb  Change GSFRAC to DSFRAC. Removed dead code. 
C  27 Apr 98 @DynTel-BAG 2-661  Removed the travel fraction factor from the 
C            final value of BYFER (= COMPEF only) 
C   9 Feb 00 @EPA-elg  Added IDX and MY parameters to HRLOOP call. 
C   2 Feb 00 @DynTel-ddj 1-032, Added CALL to T2CRT50K. 
C  26 Jan 00 @DynTel-MLA 1-011  Added var. HEF, call to DBEXHRS. 
C  29 Nov 99 @DynTel-bag 1-006  Added grams/mile index, GPM, for 
C            start emissions. 
C   4 Oct 99 @DynTel-BAG 2-694 Renamed H94OBD to OBD_HI_FRCN 
C   6 Sep 99 @DynTel-HXQ 2-684 The code to change unit for start 
C            emssion from g/mi to g/start is moved to HRLOOP.FOR 
C   5 Aug 99 @DynTel-ddj 2-699  Removed initialization of BYTAM. 
C  02 Jun 99 @DynTel-BG 2-000 Put ICY in parameter list of HRLOOP 
C  26 May 99 @DynTel-HXQ 2-462 Delete code to call EFNFTP 
C  19 May 99 @DynTel-HXQ 2-650 replace the BMYTPD with STARTSPERDAY 
C  18 May 99 @DynTel-BG 2-650 Added START(=ISR=1), RUN(=ISR=2), and 
C            USE_SOAK parameters; USE_SOAK is used with the SOAK_FLAG to 
C            determine if the start emissions soak time corrections 
C            are to be applied 
C  18 May 99 @DynTel-HQ 2-698 Delete the code to calculate SCF 
C  6  May 99 @DynTel-MLA 2-701 Initialized VRUN and VSTART to 0.0 
C  15 Apr 99 @DynTel-HQ 2-698 add the code to calculate SCF 
C  07 Apr 99 @Dyntel-BAG 2-000 Merged 2-694 changes with 2-701 changes. 
C  11 Mar 99 @Dyntel -RJD 2-694 Added call to H94OBD that corrects the 
C            high emitter fraction for OBD and OBD I/M programs. Made 
C            changes to the BEFWGT parameter list. 
C  18 Mar 99 @DynTel-MLA 2-701 Keep separate totals for start and running 
C            emissions. 
C  10 Mar 99 @DynTel-bag 2-664 IWEEK now passed through block data for 
C            weekday or weekend vehicle usage distributions 
C  11 Feb 99 @Dyntel -RJD 2-663 Added STARFG enabling the user to select 
C            starts emission units to be grams/mile (STARFG=GPM) or 
C            grams/start (STARFG=GPS). 
C   6 JAN 99 @DynTel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust 
C             engine start and running 
C  18 Sep 98 @DynTel -ZK  2-000 Explicit Typing 
C  14 Sep 98 @Dyntel -RJD 0-632 testing 
C  08 Sep 98 @Dyntel -RJD 0-632 The basic emission factor for all light duty 
C            gasoline vehicles and trucks is now obtained from BEFWGT instead 
C            of BEF. 
C  10 Aug 98 @Dyntel -RJD 2-663 Replaced call to BEF with call to BEFWGT for  
C            light duty gasoline vehicles. Moved calculation of BYBEF4 and  
C            tampering adjustment from BEF to EFCALX. 
C 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY,INERR 
C    common blocks: 
C    /CEVBMY/ BMYMPD,BMYTPD 
C    /CITRV2/ RVPX 
C    /CVNAMS/     ILDGT   (IVTYPE.I) 
C    /CUMCOM/     CUMMIL 
C    /FLAGS1/     TAMFLG 
C    /FLAGS3/     RLFLAG 
C    /FLAGS5/     SOAK_FLAG 
C    /GSFCOM/     DSFRAC 
C    /MAXIMA/     MAXIH,MAXPOL,MAXYRS 
C    /MYRCAL/     MEVMYR, TF 
C    /MYRSAV/     NEWCUM 
C    /REGION/     IREJN 
C    /RVPEX2/     RVPCF 
C    /UNIT_PARAM/ GPM 
C    /VMXCOM/     REGMIX,VMTMIX,VMLDGT 
C    /VVNAMS/     VVGAS, VVHEAVY,VVSTARTS  (IVTYPE.I) 
C    /VVTYPS/     VTGASBUS, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /BYMYC2/ BYBEF4,BYFER 
C    /RESUL1/ EFFTP 
C    /RESUL3/ VFTP 
C    /RESUL4/ RUNEF,STARTEF 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  COMPEF   R    composite emission factor for IDX: summed across IDX to get 
C                EFFTP for the given IP & IV 
C  HEF      R    Hourly emission factor 
C  RUNMYR   R    running emissions for one model year 
C  RVOPCF   R    exhaust HC & CO cf for RVP and open loop technology, IV=1,4 
C  STARTMYR R    start emissions for one model year 
C  VMTAGE   R    see parameter dictionary 
C 
C 
C  Notes: 
C 
C  For MOBILE4 EFCALX was modified to save by model year 
C  computations for the by model year tables. 
C  For MOBILE4.1 Oxygenated and certification fuels support was 
C  removed and in particularly the OPENLP correction. 
C  March 10, 1993 IY was added to the argument list to guarantee 
C  the 1990 inventory emission factors remain the same in newer 
C  versions of MOBILE as in MOBILE5. 
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor  
C  calculation  function. 
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4  
C  calculation to remove improperly making a SALHCF adjustment to by-model 
C  Sep-21-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar 
C  year removed IY from the call to TAMPER 
C  Oct-11-1994 @ CSC-bsg Request 467 Made correction to second call to  
C  EFNFTP function, the order of parameters was incorrect 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  26 Feb 98 Dyntel@RJD (after gjr) 2-204 Removed IY from parameter list and 
C            calls to subroutines refuel, hccalx and function bef.  This is  
C            due to a change in the way July emissions will be calculated. 
C            No interpolation will be used. 
C  11 Mar 98 @DynTel-RJD 2-622  Added calls to EFCVT so that emission 
C            factors for composite LDGT vehicles are not calculated here.  
C            Moved calculation of VMTMIX for LDGT vehicles to TFCALX. 
C  22 Apr 98 @DynTel-RJD 2-654 Added ISR loops to distinguish between   
C            starts and running emissions. Added an extra dimension to 
C            the BYBEF4 array so that starts and running emission can be 
C            stored simultaneously. 
C  27 Apr 98 @DynTel-BAG 2-661  Removed the travel fraction factor from the  
C            final value of BYFER (= COMPEF only) 
C  30 Jun 98 @DynTel-BAG 2-000  Removed the I/M credit calculations, 
C            REDSUM, BYIMCR, and BYIMRE 
C  09 July 98 @Dyntel -MMS 2-666 Added VVSTARTS in order to seperate 
C            computation of start and running modes 
C 
      USE DATABASE, ONLY: DBSELPOL 
      IMPLICIT NONE 
      INCLUDE 'BYMYC1.I' 
      INCLUDE 'BYMYC2.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'DSFCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RESUL3.I' 
      INCLUDE 'RESUL4.I' 
      INCLUDE 'RVPEX2.I' 
      INCLUDE 'SPEED6.I' 
      INCLUDE 'VMXCOM.I' 
      INCLUDE 'STDIST.I' 
      INCLUDE 'UNIT_PARAM.I' 
      INCLUDE 'WEEKLY.I' 
C 
C     Declare external functions. 
C 
      REAL    EFCVT 
C 
C     Declare parameter list. 
C 
      INTEGER, INTENT(IN)    :: ICY 
      INTEGER, INTENT(INOUT) :: INERR 
C 
C     Declare local variables. 
C 
      REAL    COMPEF 
      REAL    HEF(24,5) 
      REAL    RUNMYR, STARTMYR 
      REAL    VMTAGE 
      INTEGER IP 
      INTEGER ISR 
      INTEGER ITTT 
      INTEGER IV 
      INTEGER ACTUAL_MY 
      INTEGER MY 
      INTEGER IDX 
      INTEGER JDX 
      INTEGER, PARAMETER :: USE_SOAK=1 
      INTEGER, PARAMETER :: RUN=1 
      INTEGER, PARAMETER :: START=2 
C 
C     HEF is used to pass the database output values to DBEXHRS. If 
C     DBEXHRS is called from here, the emission factors will all be 
C     set to zero. 
C 
      HEF = 0.0 
C 
      IF(NEWCUM.NE.0) CALL GETCUM 
      NEWCUM=0 
C 
C  CUMMIL, a parameter array, is recast here for the BYMY tables as BYCUM, 
C  with zero replacing the potential CUMMIL, when TF=0.0. 
C 
      DO IV=1,MAXVEH 
        DO JDX=1,MAXYRS 
          IDX=(MAXYRS+1)-JDX 
          IF(TF(IDX,IV).EQ.0.0) BYCUM=0.0 
          IF(TF(IDX,IV).GT.0.0) BYCUM=CUMMIL 
        END DO 
      END DO 
C 
C  Setup for composite light duty gas truck. 
C 
C  Calculate this scenario's tampering and refueling offsets and RVP correction 
C  factors. 
C 
C  Tamper Subroutine is called twice. The first time is the normal call.  It calculates 
C  the tampering effects based on the user provided scenario parameters. 
C  The second call to TAMPER assumes No I/M and No ATP effects.  It is called so that 
C  the difference between the normal run and the No I/M or ATP can be calculated. 
C  This difference is subtracted from the 1981-1993 model year, and subsequent 
C  Tier1 vehicle emission rates to produce a no tampering emission rate. 
C  This is necessary because these emission rates include the effects of tampering. 
C 
C  Tampering effects for all three pollutants are called from EFCALX.  Tampering 
C  for HC is not called from HCCALX. 
C 
      ITTT = 1                               !Regular Tampering routine call. 
      IF(TAMFLG.GE.0) CALL TAMPER(ICY,ITTT) 
C 
      ITTT = 2                               !No I/M or ATP effects Tampering call. 
      IF(TAMFLG.GE.0) CALL TAMPER(ICY,ITTT) 
C 
      DO IDX = 1,25 
        MY=ICY+IDX-MAXYRS 
C 
        IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
          MY = MIN(CAAMY,MY) 
        ENDIF 
C 
        DO IV = 1,28 
          DO IP = 1,3 
            DO ISR = 1,2 
                CALL TEX1(IP,IV,IDX,ISR,MY) 
            END DO 
          END DO 
        END DO 
      END DO 
C 
      IF(TAMFLG.EQ.0) CALL QUITER(0.,0,119,INERR) 
C 
      IF(DBSELPOL(1).EQ.2) CALL REFUEL(ICY,INERR) 
      CALL RVPEXH(ICY) 
C 
C  Call OBD_HI_FRCN to calculate the OBD high emitter fraction arrays for 1994+ 
C  model year LDGV/Ts. 
C 
      IF (ICY.GE.1994) CALL OBD_HI_FRCN(ICY,INERR)
C 
C   Build Tier2 50K Certification Standard array CRT50KT2(IV,IP,IDX) 
C 
      CALL T2_CRT(ICY) 
C 
C     Finally ready for ef calculations - only for user selected pollutants. 
C     Now do FTP efs. 
C 
      IF(DBSELPOL(1).EQ.2) 
     *    CALL HCCALX(ICY,INERR) 
C 
      DO 90 IP=2,MAXPOL 
        IF(DBSELPOL(IP).EQ.1) CYCLE 
        VFTP(IP)=0.0 
        VRUN(IP)=0.0 
        VSTART(IP)=0.0 
        EFFTP(IP,ILDGT)=0.0 
C 
C       The following call to soak function, SOAKFCN, determines an array, SOAKFUN, 
C       which are used to determine hourly soak time correction factors for use 
C       in HRLOOP to determine hourly emission factors. It is dependent only on 
C       soak time, pollutant, and type of catalyst. So, calculate it for each 
C       pollutant out of vehicle, start, hour, and year loops. (See parameters in 
C       Table 4 of M6.STE.003) 
C 
        IF( SOAK_FLAG .EQ. USE_SOAK ) CALL SOAKFCN(IP) 
C 
        DO 80 IV=1,MAXVEH 
C 
          EFFTP(IP,IV)=0.0 
          RUNEF(IP,IV) = 0.0 
          STARTEF(IP,IV) = 0.0 
C 
C  Zero out BYMY save areas before VMTMIX check.  Since VMTMIX can be altered 
C  by scenario, a VMTMIX(i) > 0.0 case followed by a VMTMIX(i) = 0.0 case 
C  would lead to the latter scenario using the former scenario's BYMY data for 
C  IV=i. 
C 
          DO JDX=1,MAXYRS 
            BYBEF4(IP,JDX,IV,RUN)=0.0 
            BYBEF4(IP,JDX,IV,START)=0.0 
            BYFER(IP,JDX,IV)=0.0 
          END DO 
C 
          MY=0 
C 
C  calculate the speed correction factor for 2.5-65 m/h 
C 
          IF(VVLDGAS(IV).NE.1.AND.IV.NE.VTMC.AND.SPD_FLAG.EQ.1) 
     &                   CALL SCFTAB(IV,IP) 
C 
          DO 70 IDX=1, MAXYRS 
C 
            JDX = (MAXYRS + 1) - IDX 
            ACTUAL_MY=ICY+IDX-MAXYRS 
            MY=ACTUAL_MY 
C 
            IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
              MY = MIN(CAAMY,ACTUAL_MY) 
            ENDIF 
C 
C  If no mileage for this model year / vehicle class case, skip ef 
C  calculation, but database output should still be produced. 
C 
            IF((TF(IDX,IV).LE.0.0.AND.MEVMYR(JDX,IV).EQ.0.0) .OR. 
     &         (VMTMIX(IV).EQ.0.0)) THEN 
ccs              CALL DBEXHRS(HEF,IP,IV,ACTUAL_MY,JDX,INERR) 
ccs Add call to custom output routine 
              CALL SMKEXOUT(HEF,IP,IV,JDX) 
              ! IF (VMTMIX(IV).EQ.0.0) GOTO 80 
              CYCLE 
            ENDIF 
C 
C  Initialize emissions to zero. 
C 
            COMPEF=0.0 
            RUNMYR = 0.0 
            STARTMYR = 0.0 
C 
            VMTAGE = CUMMIL(JDX,IV) / 10000. 
            CALL NGV_EF(IP,IV,ACTUAL_MY,IDX,VMTAGE) 
C 
            CALL HRLOOP(RUNMYR,STARTMYR,IP,IV,ICY,MY,ACTUAL_MY, 
     *                  IDX,INERR) 
C 
            COMPEF = RUNMYR + STARTMYR 
C 
            BYFER(IP,JDX,IV)=COMPEF 
            EFFTP(IP,IV)=EFFTP(IP,IV)+COMPEF*TF(IDX,IV) 
            RUNEF(IP,IV)=RUNEF(IP,IV)+RUNMYR*TF(IDX,IV) 
C 
C     If start emissions are in g/mi, TF can be used for the model year 
C     weighting. If they are in g/start, use MEVMYR for the weighting. 
C 
            STARTEF(IP,IV)=STARTEF(IP,IV)+STARTMYR*TF(IDX,IV) 
C
C     Abort whenever any errors occur.
C
          IF(INERR.GT.0) GOTO 99
          
   70     CONTINUE 
C         END OF YEARLY LOOP, IDX 
C 
C  Sales fractions = 1.0 for all 4 HDV (HDG,HDD x L,H alt) groups  
C  and all MC's. 
C 
          IF(EFFTP(IP,IV).LT.0.0.AND.VMTMIX(IV).GT.0.0) 
     *    CALL QUITER(0.,IV,68,INERR) 
C 
          VFTP(IP)=VFTP(IP)+EFFTP(IP,IV)*VMTMIX(IV) 
          VRUN(IP)=VRUN(IP)+RUNEF(IP,IV)*VMTMIX(IV) 
          VSTART(IP)=VSTART(IP)+STARTEF(IP,IV)*VMTMIX(IV) 
C 
   80   CONTINUE 
C       END OF VEHICLE LOOP, IV 
C 
C  Assign composite light duty gas truck emission factor. 
C 
        IF(EFFTP(IP,2).GT.0.0.AND.EFFTP(IP,3).GT.0.0.AND. 
     *   EFFTP(IP,4).GT.0.0.AND.EFFTP(IP,5).GT.0.0) 
     *   EFFTP(IP,ILDGT)= EFCVT(ILDGT,IEFFTP,IP,INERR) 
        IF(STARTEF(IP,2).GT.0.0.AND.STARTEF(IP,3).GT.0.0.AND. 
     *   STARTEF(IP,4).GT.0.0.AND.STARTEF(IP,5).GT.0.0) 
     *   STARTEF(IP,ILDGT)= EFCVT(ILDGT,ISTARTEF,IP,INERR) 
        IF(RUNEF(IP,2).GT.0.0.AND.RUNEF(IP,3).GT.0.0.AND. 
     *   RUNEF(IP,4).GT.0.0.AND.RUNEF(IP,5).GT.0.0) 
     *   RUNEF(IP,ILDGT)= EFCVT(ILDGT,IRUNEF,IP,INERR) 
C 
   90 CONTINUE 
C     END OF POLLUTANT LOOP, IP 
C 
   99 RETURN
      END 
