      SUBROUTINE HRLOOP(RUNSUM,STARTSUM,IP,IV,ICY, 
     *                  MY,ACTUAL_MY,IDX,INERR) 
C 
C  HRLOOP calculates exhaust emission factors and 
C         transfers them to the Database Output and the 
C         Descriptive Output Routines. 
C 
C  Called by EFCALX, HCCALX 
C 
C  Calls: 
C 
C        AC_ADJ 
C        BIGCFX 
C        BIGIN2 
C        SPEED_CF 
C        DBEXHRS 
C        DEFEAT 
C        DEFEAT_DEVICE 
C        DEFEAT_STAT 
C        EF_CASE 
C        EF_HI_FRCN 
C        FUEL 
C        GBHP_GPM 
C        HDCPIC 
C        IEVPTR 
C        IVPTRT 
C        IMPROG 
C        METHANE 
C        PCLEFT 
C        P3WPOX 
C        SCFSLOW 
C        SCFTABMC 
C        SOAKCF 
C 
C  Changes:  (Last change first) 
C 
C  31 Dec 01 @EPA-BAG     Changed .GE. to >= in IMPROG conditionals 
C  15 Nov 01 AIR Task 27: DBEXHRS now includes INERR. 
C  17 Dec 01 @EPA-BAG        Remove I&M sawtooth 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  13 Feb 01 @EPA-elg increased the size of FVMT arrays to include vehicle class 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; added 
C            MY to parameter list of SOAKCF and P3WPOX 
C  19 Oct 00 @EPA-djb Task W2, Make the DBEXHRS call conditional. 
C  29 Sep 00 AIR Task 07: Bug Fix 271: Fixed logic with use of variables 
C            NGV_EF and NGV2 
C  30 Aug 00 @EPA-djb Bug 243, Removed CYCLE statements that prevented 
C            calculation of EMFAC array for *all* ISTDs. 
C  29 Aug 00 @EPA-elg Eliminated PCTREM storage.  I/M is now called in the 
C                     hourly loop and used as needed in IMPCT. 
C  16 Aug 00 @EPA-djb Bug 223, Apply the TCF after the SCF.  This requires 
C            sending the TCFRUN to the SPEED_CF subroutine. 
C  16 Aug 00 @EPA-djb Bug 227, Removed NGV start weighting and 
C            added NOXHH correction. 
C  15 Aug 00 @EPA-elg Bug 225, Altered the logic of assigning the values 
C                              to the array HIGHER. 
C  14 Aug 00 @EPA-djb Bug 226, removed the unitialized BASEF calculation. 
C  04 Aug 00 @EPA-rag Bug 203, put upper and lower limits on value of 
C            absolute humidity for NOx humidity correction 
C  07 Aug 00 @EPA-djb Bug 218, Added conditional to PCTREM initialization. 
C  27 Jul 00 @EPA-djb Bug 208, removed MAXSTD and EVASTD include. 
C                     Bug 178, moved I/M calculations outside IH loop. 
C  20 Jul 00 @EPA-ddj Bug200 Moved TAMP to inside major ISR loop, 
C            moved NGVFRAC outside of IH loop and indexed DDNOX with 
C            a dimension of 4 for IROAD. 
C  21 Jun 00 @EPA-DDJ bug162 redid bug 86 (see below; bug 86 fix was not integrated) 
C  14 Jun 00 @EPA-elg Changes for I&M programs; bug 146 
C  14 Jun 00 AIR Task 03: Modified code to properly handle temperature and 
C            RVP corrections and CO offsets. Added calls to BIGIN2 and BIGCFX. 
C  15 Jun 00 @EPA-BAG Removed INH from parameter list of METHANE subroutine; 
C                     INH is not used in the METHANE subroutine 
C  08 Jun 00 @EPA-djb Bug113 Altered ISTD/IFDS loop to calculate HIGHER 
C            and added ISTD to HIGHER array.  Made HIGHER a local array. 
C   2 Jun 00 @EPA-ddj E11, Changed varable NGVFRAC to an array (ISR,INH) 
C            and moved outside of IH loop. 
C   2 Jun 00 @EPS-ddj Bug 86, Changed conditional of FDS weighting, 
C            from VVLDGAS to VVGAS (now includes HDGVs). 
C  31 May 00 @EPA-djb Limited tampering cap to 1981+ vehicles. 
C  16 May 00 @EPA-ddj BUG27, Revised call to fuel correction routine. 
C  12 May 00 @EPA-djb Skip emission standards with no vehicles. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C   2 May 00 @EPA-ddj E-7, Changed sulfur call from SUL_COR to SULFUR. 
C  03 May 00 @EPA-djb Added code to skip start calculations for HDV. 
C  04 May 00 @EPA-bg    fix bug #72 "(ISTD.GE.1 .AND. MY.GE.1994)" rather than 
C                       "(ISTD.GT.1 .AND. MY.GE.1994)" ;Tier1 phase ins 
C  28 Apr 00 @EPA-elg   Changed Call to SPEED_CF to handle Offcycle and SFTP. 
C  27 Apr 00 @EPA-elg   Installed Tampering Effects. 
C  27 Apr 00 @EPA-elg   Fixed I/M for Start Emissions. 
C  27 Apr 00 @EPA-elg   Added NGV emission calculation code. 
C  27 Apr 00 @EPA-djb   Set BER checks back to zero. 
C  26 apr 00 @EPA -bg    removed IV from SOAKCF parameter list; IV is only used in 
C                        diagnostic write statements which are commented out 
C  24 Apr 00 AIR Task 02: Fixed SOAKCF call to include IV in parameter list. 
C            Included INERR in HRLOOP, IMPROG and EF_HI_FRCN parameter lists. 
C  17 Apr 00 @EPA-djb Set BER checks back to zero. 
C  11 Apr 00 @EPA-elg   task 31, 32 and hrloop task revisions 
C  21 MAR 00 @EPA-elg   Draft Master Version HRLOOP created. 
C  16 MAR 00 @EPA-elg   Added call to IMPROG from HRLOOP to calculate 
C                       High emitter fraction if no I/M. 
C  06 FEB 00 @EPA-elg   Remove calls to CF_BEF CF_BEF94 
C  28 JAN 00 @EPA-elg   BEGIN COMPLETE RESTRUCTURING OF HRLOOP. 
C  28 JAN 00 @EPA-elg   FTP Button code moved to FTPBUT.for. 
C  10 Mar 00 @DynTel-MLA 1-011 database code 
C  13 Jan 00 @DynTel-ddj Added call to AC_ADJ.  Limited NOXHH 
C            correction to 0 to accomodate revised celling for ABSHUM 
C            from 140 to 528 grains/lb. 
C  22 Dec 99 @Dyntel-nh removed diagnostic parameters in DEFEAT call. 
C  06 DEC 99 @DYNTEL-BAG Added grams/start index, GPS, for start emissions. 
C  07 Dec 99 @Dyntel-nh added code to invoke DEFEAT_STAT 
C  19 Nov 99 @DynTel-ddj 1-002 Changed conditional of STDWT,DEFEAT_STAT 
C            initialization from 1981 to 1994. 
C  02 Nov 99 @Dyntel-nh 1-004 Added code for DEFEAT DEVICE 
C  14 Oct 99 @DynTel-ddj 2-714  Removed addition of OFFCO.  OFFCO is now 
C            included in OMTCF. 
C  13 Oct 99 @DynTel-ddj 2-714  Maintained units for emission 
C            factors in grams except for HOUREF. 
C  28 Sep 99 @DynTel-ddj for EPA  Added ISTD index to OFFMTH. 
C  24 Aug 99 @DynTel-MLA 2-684 Added call to DBEXHRS. 
C   8 Sep 99 @DynTel-BAG 2-660 Renamed HIEMIT to EF_HI_FRCN 
C  25 Aug 99 @DynTel-ddj 2-710 Added conditional statement to bypass 
C            tampering calculation when TAMFLG > 0.  Also removed 
C            T_BEFSUM from 50 ISR do loop and added initialization of 
C            T_BEFSUM to 100 ISR do loop.   Added IDX back to SCFSLOW. 
C  10 Sep 99 @Dyntel-HXQ 2-705 Deleted the code for SUM which is used 
C             to modify LEVIMP(ISTD,LMY,IV) 
C   6 Sep 99 @DynTel-HXQ 2-684 hourly start emissions are summed 
C            over the day in g/start. Their unit is changed to g/mi 
C            if required 
C  17 Aug 99 @DynTel-ddj Removed IDX from SCFSLOW parameter list. 
C  27 Jul 99 @DynTel-ddj Moved ISR loop to front of IH loop and added 
C            ISTD loop after INH loop. 
C  13 Jul 99 @DynTel-ddj 2-699 Changed implementation logic of FTP_FLAG 
C            back so that FTP_FLAG is equal to 1 when invoked and 0 
C            otherwise. 
C  18 Jun 99 @DynTel-ddj 2-000 Added write statements to allow writing 
C            out of FTP emission #'s; these are currently commented out 
C   1 Jun 99 @DynTel-ddj 2-700 Changed implementation logic of FTP_FLAG 
C            to be consistent with other flags. ICY added to parameter list. 
C  28 May 99 @DynTel -bg 2-650 Added code to write out results for a check with 
C            EPA documents. Output code is commented out. 
C  19 May 99 @DynTel-ddj 2-700 Added coding for implementation of FTP 
C            exhaust output option. 
C  18 May 99 @Dyntel-HXQ 2-650 add the code of changing unit of start 
C                        emission from gram/start to gram/mile 
C  18 May 99 @Dyntel-BG  2-650 Rearranged correction factors in hourly loop. 
C            RVPCF now applied in first loop. Created constant parameters 
C            for INH(NORMAL=1 and HIGH=2) and ISR(RUN=1 and START=2) 
C  12 May 99 @DynTel-HXQ 2-698 Move SCF calculation out of IH Loop. 
C  10 May 99 @Dyntel-BG  2-000 Rearranged correction factors in hourly loop. 
C            RVPCF now applied in first loop. Created constant parameters 
C            for INH(NORMAL=1 and HIGH=2) and ISR(RUN=1 and START=2) 
C   6 May 99 @Dyntel-MLA 2-000 CHanged "IV.LT.6" to "VVLDGAS(IV).EQ.1". 
C            Also changed STARTDIST - it now contains data for 24 'hourly' 
C            intervals instead of the 14 intervals it held previously. 
C   3 May 99 @DynTel-ddj, Removed soak correction from hourly emission array. 
C            Soak correction already applied with temperature correction 
C            factor. Removed BEFSUM. It has been replaced by T_BEFSUM. 
C  10 Apr 99 @Dyntel-HXQ 2-698 Integrated this code with speed correction 
C  22 Apr 99 @DynTel-ddj Modified If conditional for when additive CO 
C            correction is applied.  Also changed the assignment of 
C            EF so same BASEF value is used for each IH. 
C  01 Apr 99 @DynTel-DDJ Integrated this code into MOBILE6, EFCALX and 
C            HCCALX. 
C  01 Dec 98 @DynTel-RJD Removed INERR from HIEMIT argument list. 
C  22 Oct 98 @DynTel-bag  Restructured soak time corrections 
C            calculations to reduce MOBILE6 run time. SOAKADJ, SOAKWGT, 
C            and CALSOAK have been replaced with SOAKCF and SOAKFCN. The 
C            call to P3WPOX to find calalyst type vehicle fractions has 
C            been moved to HRLOOP and the ratio (CAT=P3W+POX) is passed 
C            to SOAKCF as a parameter.  Also, the model year dependence 
C            is now passed as a common variable in MYCODE. 
C   4 Aug 98 @DynTel-MLA  New routine added for Mobile6. 
C 
C  Input on call: 
C 
C    Constants: 
C    MAXVEH from include file IVTYPE.I 
C    MAXEMIT, MAXIH from MAXIMA.I 
C 
C    argument list: IP,IV,ICY,MY,IDX,INERR 
C 
C  Output on return: 
C 
C    argument list: RUNSUM,STARTSUM 
C 
C    common blocks: 
C     /ALUHIN/  ABSHUM 
C     /BASEQ5/  EMFAC, HIGHER 
C     /BASEQ9/  MAXLYR 
C     /CEVBMY/  BMYMPD 
C     /CUMCOM/  CUMMIL 
C     /DDDATA/  FIRST_DEFEAT 
C     /FLAGS1/  SCFUSER 
C     /FLAGS2/  IMFLAG 
C     /FLAGS5/  AC_FLAG, SOAK_FLAG, SPD_FLAG 
C     /IMBLOCK1/ FUELCF 
C     /INJECT/  TBI, PFI 
C     /IVPCOM/  IVPTRS, IVPTRT 
C     /IVTYPE/  VTGASBUS, VVHDGV, VVLDGAS 
C     /MAXIMA/  MAXIH 
C     /NGV/     NGVFLG1,NATGAS,NGVFRAC 
C     /OMTCOM/  TCFSTART,TCFRUN,TCF12 
C     /OFFSET/  OFFSET 
C     /RVPEX2/  RVPCF 
C     /SPEED9/  FVMT, HVMT 
C     /STDIST/  STARTDIST, STARTSPERDAY 
C     /TAMOUT/  TEX 
C     /TEMPS/   TEMHRLY 
C     /WEEKLY/  IWEEK 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  AC_CORR  R    A/C Correction Factor 
C  ADDCF    R    Local variable that holds final additive correction factors 
C  AC_DUM   R    Emission factor by roadway type passed to A/C correction 
C                factor subroutines. 
C  BASEF    R    The basic EF used to calculate the hourly EF 
C  CAP      R    The limit on the fraction of the basic exhaust emission 
C                rate that can be tampering excess. 1=HC, 2=CO, 3=NOX 
C  CAT      R    Catalyst fractions 
C  CONV     R    Heavy-duty Diesel Vehicle Conversion Factor for Defeat 
C                Device Offset 
C  DDNOX    R    Defeat Device Offset 
C  EF       R    Temporary Variable holding the Basic Emission Factor 
C  NGV_EF   R    Temporary Variable holding the NGV Emission Factor 
C  FDSWT    R    Fuel delivery system weighting factor 
C  HEF(24)  R    Accumulates EFs for FDS/emitter groups for each MY 
C  HIGHER   R    Stored High Emitter Fractions 
C  IDX      I    Model year index (chronological order) 
C  IFDS     I    Fuel delivery system loop control variable 
C  IH       I    Hourly loop control variable 
C  INH      I    Emitter category loop control variable 
C  INHWT    R    Normal/High emitter weighting factor 
C  IP       I    Pollutant loop control variable. 
C  ISR      I    Start/running loop control variable 
C  IV       I    Vehicle type loop control variable 
C  JDX      I    Model year index (reverse chronological order) 
C  MPHR     R    Miles of VMT per Hour  (not speed) 
C  NGVFRAC  R    NGV Fleet Fraction 
C  NOXHH    R    NOx humidity correction factor 
C  OFFMTH   R    methane emission offset 
C  P3W      R    Percent 3-way catalysts 
C  POX      R    Percent oxidation catalysts 
C  RVPADJ   R    Fuel RVP adjustment factor 
C  SOAKCF   R    Final soak time correction factors. 
C  SKCF     R    Temperary variable to hold the soak correctin factor. 
C  STARTER  R    Variable holding the FTP start emission factor.  It is 
C                when calculating the HC oxyfuels effects for start emissions. 
C  STDWT    R    Certification standards weighting factor 
C  TEX      R    Tampering Offset 
C  TEF      R    Temporary Variable holding the Basic emission factor 
C  TPHR     R    Trips per hour 
C  VMTAGE   R    Vehicle miles travelled 
C 
C 
C  Notes: 
C 
      USE DDDATA 
      USE DATABASE, ONLY: DBFLAG 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ALUOUT.I' 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'BASEQ5.I' 
      INCLUDE 'BYMYC2.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IMBLOCK.I' 
      INCLUDE 'INJECT.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'NGV.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'OMTCOM.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RVPEX2.I' 
      INCLUDE 'SOAK.I' 
      INCLUDE 'SPEED2.I' 
      INCLUDE 'SPEED6.I' 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'STDIST.I' 
      INCLUDE 'TAMOUT.I' 
      INCLUDE 'TEMPC4.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'WEEKLY.I' 
C 
C     Declare parameters. 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(INOUT) :: INERR 
C 
      REAL,INTENT(OUT)    ::  RUNSUM 
      REAL,INTENT(OUT)    ::  STARTSUM 
C 
C     Declare external functions. 
C 
      INTEGER, EXTERNAL :: IEVPTR 
C 
      LOGICAL, EXTERNAL :: DEFEAT_DEVICE 
      LOGICAL                  DEFEATT 
C 
      REAL, EXTERNAL :: AC_ADJ 
      REAL, EXTERNAL :: SPEED_CF 
      REAL, EXTERNAL :: EF_CASE 
      REAL, EXTERNAL :: EF_HI_FRCN 
      REAL, EXTERNAL :: HDCPIC 
      REAL, EXTERNAL :: IMPROG 
      REAL, EXTERNAL :: PCLEFT 
      REAL, EXTERNAL :: SOAKCF 
C 
C     Declare local variables. 
C 
      REAL     AC_CORR 
      REAL     AC_DUM 
      REAL     AC_DUMM 
      REAL, DIMENSION(2)           :: BASEF 
      REAL, DIMENSION(3), SAVE     :: CAP=(/0.40,0.15,0.20/) 
      REAL     CAT 
      REAL     CONV 
      REAL, DIMENSION(4)           :: DDNOX 
      REAL     DDNOXX 
      REAL     EF 
      REAL, DIMENSION(3)           :: FDSWT 
      REAL, DIMENSION(24,5)        :: HEF 
      REAL, DIMENSION(2,3,2)       :: HIGHER 
      REAL     HUMIDITY 
      REAL, DIMENSION(2)           :: INHWT 
      REAL     MPHR 
      REAL, DIMENSION(2,2)         :: NGV1 
      REAL     NGV2 
      REAL     NGV_EF 
      REAL     NOXHH 
      REAL     P3W 
      REAL     POX 
      REAL     SKCF 
      REAL, DIMENSION(10)          :: STDWT 
      REAL     TAMP 
      REAL     TEF 
      REAL     TPHR 
      REAL     VMTAGE 
C 
      INTEGER  IFDS 
      INTEGER  IG5 
      INTEGER  IH 
      INTEGER  INH 
      INTEGER  IROAD 
      INTEGER  IROADD 
      INTEGER  ISR 
      INTEGER  ISTD 
      INTEGER  IVTEMP 
      INTEGER  JDX 
      INTEGER  JSTD 
      INTEGER  LMY 
      INTEGER  NGVPT 
      INTEGER  ST_INDX 
      INTEGER, DIMENSION(4)   :: RDMAP=(/3,2,1,4/) 
      INTEGER, PARAMETER      :: USE_SOAK=1 
      INTEGER, PARAMETER      :: RUN=1,START=2 
      INTEGER, PARAMETER      :: NORMAL=1,HIGH=2 
      INTEGER, PARAMETER      :: NOX=3 
C 
C  Call BIGIN2 to initialize the pointers, arrays and variables used 
C  by METHANE, BIGCFX and the IH loop below. 
C 
      JDX = (MAXYRS + 1) - IDX 
      CALL BIGIN2(IV,MY,JDX,IDX,IP) 
C 
C  Initialize DEFEAT values, when necessary  
C 
      IF (FIRST_DEFEAT) THEN 
        FIRST_DEFEAT=.FALSE. 
        CALL DEFEAT_STAT(IOUREP) 
      END IF 
C 
      DEFEATT = IP.EQ.NOX .AND. DEFEAT_DEVICE() .AND. VVHDDV(IV).EQ.1 
C 
C     Heavy-duty Diesel Vehicle NOx Defeat Device Calculations 
C 
      DO IROAD = 1, 4 
        IF (DEFEATT) THEN 
          IROADD=RDMAP(IROAD) 
          CALL DEFEAT(ACTUAL_MY,ICY,IV,IROADD,DDNOXX) 
          DDNOX(IROAD) = DDNOXX 
        ELSE 
          DDNOX(IROAD) = 0.0 
        END IF 
      END DO 
C 
      RUNSUM    = 0.0 
      STARTSUM  = 0.0 
      LMY = MY - 1993 
      IF(LMY.LE.1) LMY = 1 
      IF(LMY.GT.MAXLYR) LMY = MAXLYR 
C 
      IF((VVLDGAS(IV).EQ.0 .AND. VVLDDSL(IV).EQ.0) .OR. 
     &                      MY.LT.1994) THEN 
         STDWT(1) = 1.0 
         DO ISTD = 2, 10 
            STDWT(ISTD) = 0.0 
         END DO 
      ELSE 
         DO ISTD = 1, 10 
            STDWT(ISTD) = LEVIMP(ISTD,LMY,IVPTRC(IV)) 
         END DO 
      END IF 
C 
      VMTAGE = CUMMIL(JDX,IV) / 10000. 
C 
C     Calculate catalyst fractions for soak time correction factor 
C     calculation in hourly loop; gas vehicles only have catalysts. 
C 
      IF(VVGAS(IV) .EQ. 1) THEN 
         CALL P3WPOX(MY,IV,P3W,POX) 
         CAT = POX+P3W 
      ELSE 
         CAT = 0.0 
      ENDIF 
C 
C     Initialize the array used to accumulate the emission 
C     factors for the different emitter and FDS categories. 
C 
      IF(NGVFLG1.EQ.1) THEN 
         IF(ACTUAL_MY.LT.NGV_START) THEN 
            NGVFRAC = 0.0 
         ELSE 
            NGVPT = ACTUAL_MY - NGV_START + 1 
            NGVFRAC = NGVPHIN(NGVPT,IV) 
         END IF 
      ELSE 
         NGVFRAC = 0.0 
      END IF 
C 
      DO IH = 1,MAXIH 
        DO IROAD = 1,5 
          HEF(IH,IROAD) = 0.0 
        END DO 
      END DO 
C 
C     NGV Emission Calculation 
C 
      DO ISR = RUN, START 
         DO INH = NORMAL, HIGH 
               IF(NGVFLG1.EQ.1) THEN 
                  IF(NGVFRAC .GT.0.0) THEN 
                     NGV1(ISR,INH) = NATGAS(IV,IP,IDX,ISR,INH) 
                  ELSE 
                     NGV1(ISR,INH) = 0.0 
                  ENDIF 
               ELSE 
                  NGV1(ISR,INH) = 0.0 
               ENDIF 
         END DO 
      END DO 
C 
C     NOx humidity Correction carried over from MOBILE5. 
C 
      IF((VVLDGAS(IV).EQ.1. OR. IV.EQ.VTMC) .AND. IP.EQ.3) THEN 
          HUMIDITY=MAX(ABSHUM,21.)    !humidity test data lower limit 
          HUMIDITY=MIN(ABSHUM,124.)   !humidity test data upper limit 
          NOXHH = 1.0 - (HUMIDITY - 75.0)/263.2 
          NOXHH = MAX(NOXHH, 0.0) 
      ELSE 
          NOXHH = 1.0 
      END IF 
C 
C     Find the percentage of the current model year's vehicles 
C     using different fuel delivery system technologies. 
C 
      IF (VVGAS(IV).EQ.0.OR.MY.LT.1981) THEN 
C 
C       Only one FDS tech. group for Diesels, HDGV, MC, and 
C       all pre-1981 gasoline powered vehicles. 
C 
        FDSWT(1) = 1.00 
        FDSWT(2) = 0.00 
        FDSWT(3) = 0.00 
C 
      ELSEIF( MY>=1994 ) THEN 
C 
C       Post 1993 model years need only a single loop over the 
C       fuel delivery system; carbureted cars have no contribution 
C       and PFI numbers are equal to TBI numbers. So, use PFI 
C       numbers only. 
        FDSWT(1) = 0.00 
        FDSWT(2) = 0.00 
        FDSWT(3) = 1.00 
C 
      ELSE 
C 
C       1=Carb, 2=TBI, 3=PFI 
C 
        IVTEMP = IVPTRT(IV) 
        IF (IV.EQ.VTGASBUS) IVTEMP = 14 
        IG5 = IEVPTR(MY,IVTEMP,2) 
        FDSWT(2) = TBI(IG5,IVTEMP) 
        FDSWT(3) = PFI(IG5,IVTEMP) 
        FDSWT(1) = 1.00 - FDSWT(2) - FDSWT(3) 
C 
      ENDIF 
C 
C     Initialize the fraction of high emitters 
C 
      DO IFDS=1,3 
        DO ISTD=1,2 
          HIGHER(HIGH,IFDS,ISTD) = 0.0 
        END DO 
      END DO 
C 
C     Calculate and Store the Basic Emission Factors and High Emitter Fraction. 
C 
      DO IFDS = 1,3 
C 
        IF (FDSWT(IFDS).LE.0.00001) CYCLE 
C 
C     Only gasoline-powered vehicles are differentiated by 
C     fuel delivery system technologies. Diesels are not. 
C 
        IF(IFDS.GT.1.AND. 
     *     (VVGAS(IV).EQ.0.OR.MY.LT.1981)) CYCLE 
C 
C     1994 and newer model years are split into 11 emission standard 
C     groupings.  These groupings are: 
C 
C     ISTD = 1 : Tier 0 
C            2 : Intermediate Tier 1 
C            3 : Tier 1 
C            4 : Tier 2 
C            5 : Intermediate Transitional Low Emission Vehicle 
C            6 : Transitional Low Emission Vehicle 
C            7 : Intermediate Low Emission Vehicle 
C            8 : Low Emission Vehicle 
C            9 : Transitional Ultra Low Emission Vehicle 
C           10 : Ultra Low Emission Vehicle 
C           11 : Zero Emission Vehicle 
C 
C     The ISTD loop does *not* include the ZEVs, since by definition, 
C     the emissions from these vehicles are zero.  Leaving them out of 
C     the loop assures that the fraction of vehicles which are ZEVs 
C     are removed from the emission estimate. 
C 
        DO ISTD = 1,10 
C 
          JSTD = ISTD            !JSTD=1 : Tier 0 
          IF(JSTD.GT.2) JSTD=2   !JSTD=2 : Tier 1+ 
C 
C     Find the proportion of high and normal emitters in 
C     this MY/IP/IV/IFDS/ISTD grouping.  Always calculate 
C     the emitter fractions for the first two ISTD groups. 
C 
            IF (VVLDGAS(IV).EQ.0.OR.MY.LT.1981) THEN 
              HIGHER(HIGH,IFDS,JSTD)   = 0.0 
            ELSE 
              IF(VVLDGAS(IV) .EQ. 1 .AND. MY .LE. 1993) THEN 
                HIGHER(HIGH,IFDS,JSTD) = 
     *            EF_HI_FRCN(MY,IP,IV,VMTAGE,IFDS,1,INERR) 
              ELSE 
                IF(VVLDGAS(IV) .EQ. 1 .AND. MY.GE.1994) THEN 
C 
C     For cases of 1994+ MY's and NO I/M.  The ARRAY 
C     OBDHINIM returns the fraction of HIGH emitters in the 
C     fleet After OBD is applied. Tier 0 (ISTD=1) have no OBD. 
C 
                  IF(ISTD.EQ.1) THEN       !Tier 0 
                    HIGHER(HIGH,IFDS,1) = 
     *                EF_HI_FRCN(MY,IP,IV,VMTAGE,IFDS,ISTD,INERR) 
                  ELSE                     !Tier 1+ 
                    HIGHER(HIGH,IFDS,2) =  OBDHINIM(JDX,IP,IV) 
                  ENDIF 
                END IF 
              END IF 
            END IF 
C 
          HIGHER(HIGH,IFDS,JSTD) = 
     *         AMIN1(HIGHER(HIGH,IFDS,JSTD),1.0) 
          HIGHER(NORMAL,IFDS,JSTD) = 
     *         1.0 - HIGHER(HIGH,IFDS,JSTD) 
C 
          DO INH = 1,2 !1:Normal, 2:High 
C 
C           Skip Normal/High cases with no fraction. 
C 
            IF((INH.EQ.2) .AND. 
     &         (VVLDGAS(IV).NE.1 .OR. MY.LT.1981)) CYCLE 
C 
             DO ISR = 1,2 
C 
C     Heavy duty vehicles and buses do not have engine start emissions. 
C 
                BASEF(ISR) = 0.0 
                IF(VVHEAVY(IV).EQ.1 .AND. ISR.EQ.2) CYCLE 
C 
C     The basic emisssion factors are obtained from EF_CASE. 
C 
                IF(VVLDGAS(IV).EQ.1) THEN 
                  IF(ISTD.EQ.1 .AND. MY.LT.1994) THEN 
                    BASEF(ISR) = EF_CASE(MY,JDX,IP,IV,ISR, 
     *                                   VMTAGE,ISTD,INH,IFDS) 
                    IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0 
                  ELSE 
                    IF(ISTD.GE.1 .AND. MY.GE.1994) THEN 
                      BASEF(ISR) = EF_CASE(MY,JDX,IP,IV,ISR, 
     *                                     VMTAGE,ISTD,INH,IFDS) 
                      IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0 
                    ELSE 
                      CYCLE 
                    ENDIF 
                  ENDIF 
                ELSE 
                  BASEF(ISR) = EF_CASE(MY,JDX,IP,IV,ISR, 
     *                                 VMTAGE,ISTD,INH,IFDS) 
                  IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0 
                ENDIF 
C 
C     Add Tampering Offset to the basic emission factors. 
C 
C     First compute Tampering Offset. 
C 
                TAMP = TEX(IV,IP,ISR,IDX) 
C 
                IF(ABS(TAMP).GT.BASEF(ISR)*CAP(IP) .AND. 
     *                                 MY.GE.1981) THEN 
                  IF(TAMP.GT.0) THEN 
                     TAMP = BASEF(ISR)*CAP(IP) 
                  ELSEIF(TAMP.LE.0) THEN 
                     TAMP = -1.0*BASEF(ISR)*CAP(IP) 
                  ENDIF 
                ENDIF 
C 
                BASEF(ISR) = BASEF(ISR) + TAMP 
C 
C     METHANE subroutine returns the basic emission factor in the HC 
C             speciation units requested by the user. 
C 
                IF(IP.EQ.1) THEN 
                  IF(IV.LE.5) THEN 
                    CALL METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,BASEF) 
                  ELSEIF(IV.GE.6 .AND. ISTD.EQ.1) THEN 
                    CALL METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,BASEF) 
                  ENDIF 
                ENDIF 
C 
                BASEF(ISR) = BASEF(ISR) * NOXHH 
C 
C     EMFAC now contains the Base Emission Factor corrected for 
C           Tampering, Methane and Humidity. 
C 
                EMFAC(ISR,ISTD,IFDS,INH) = BASEF(ISR) 
C 
             END DO   !ISR LOOP 
C 
             CALL FUEL(MY,IP,IV,ICY,IFDS,INH,IDX,ISTD,BASEF,FUELCF) 
C 
          END DO      !ISTD 
        END DO        !INH LOOP 
      END DO          !IFDS 
C 
      CONV=HDCPIC(ACTUAL_MY,IV) 
C 
C     Begin Main Looping Structure of HRLOOP. 
C 
      DO IH = 1,MAXIH 
C 
C       Call BIGCFX to compute the Start and Running Temperature Corrections 
C       and Offsets applicatable for this hour. 
C 
        CALL BIGCFX(IDX,MY,IV,IP,IH) 
C 
        DO ISR = 1,2 
C 
C       Heavy duty vehicles and busses do not have engine start emissions. 
C 
          IF(VVHEAVY(IV).EQ.1 .AND. ISR.EQ.START) CYCLE 
C 
C     Calculation of trips per hour and miles driven in this hour 
C 
          IF(ISR.EQ.START) THEN 
             TPHR =  STARTSPERDAY(JDX,IV,IWEEK) * 
     *              STARTDIST(IH,IWEEK)               ! Starts per hour 
             MPHR = BMYMPD(JDX,IV) * HVMT(IH)         ! Miles per hour (VMT) 
          ENDIF 
          DO IROAD = 1,4 
C 
C     Standards, Fuel Delivery System and Normal/High loops 
C 
            DO ISTD = 1,10 
              JSTD = ISTD 
              IF(JSTD.GT.2) JSTD=2 
              IF (STDWT(ISTD).LE.0.00001) CYCLE 
              DO IFDS = 1,3       ! 1=carb,  2=tbi,  3=pfi 
                IF (FDSWT(IFDS).LE.0.00001) CYCLE 
                DO INH = 1,2      ! 1=normal,  2=high 
                  IF (HIGHER(INH,IFDS,JSTD).LE.0.00001) 
     *                CYCLE 
                  INHWT(INH) =  HIGHER(INH,IFDS,JSTD) 
C 
                  BASEF(ISR) = EMFAC(ISR,ISTD,IFDS,INH) 
C 
C  *******************   START EMISSIONS SECTION  ******************** 
C 
               IF(ISR.EQ.START .AND. IROAD.EQ.1) THEN !FOR START EMISSIONS 
C 
C     Apply start emissions temperature, offset and soak corrections, 
C     as required. If no saok is specified, then use a temperature  
C     correction corresponding to a 12-hour soak. Thereafter, apply 
C     the fuel and other corrections. 
C        
                   IF(SOAK_FLAG.EQ.USE_SOAK) THEN 
                     EF=0.0 
                     NGV_EF=0.0 
                     DO ST_INDX=1,69 
                        SKCF = SOAKCF(MY,CAT,IH,IWEEK,ST_INDX) 
                        EF = EF + (BASEF(START) * TCFSTART(ST_INDX) 
     *                     + OFFSET) * SKCF 
                        NGV_EF = NGV_EF + NGV1(ISR,INH) * SKCF 
                     END DO 
                    ELSE 
                     EF = BASEF(START) * TCF12 + OFFSET 
                     NGV_EF = NGV1(ISR,INH) 
                    ENDIF 
C 
                    EF = EF * FUELCF(START,INH,ISTD,IFDS,IH) 
C 
C     Apply the NOx humidity correction factor to NGVs. 
C 
                    NGV_EF = NGV_EF * NOXHH 
C 
                    EF = NGV_EF*NGVFRAC + EF*(1.0-NGVFRAC) 
C 
C     For PRE-1981 gas vehicles for all model years obtain % of emission factor 
C     left after I/M program effect is computed for the IH loop. Post-1980 vehicles 
C     are given a I&M credit in grams/mile to high emission vehicles only. Default 
C     (non-I/M, diesel, HD & MC) is no effect. 
C 
                    IF(IMFLAG>1) THEN 
                      IF(VVLDGAS(IV)==1.AND.MY>=1981.AND.INH==HIGH) THEN 
                        EF = EF-IMPROG(ICY,MY,ACTUAL_MY,JDX,IP,IV, 
     *                   VMTAGE,ISR,IFDS,ISTD,INERR,IDX,IH,IROAD,CAT,EF) 
                       ELSEIF((VVLDGAS(IV).EQ.1 .AND. MY.LT.1981) .OR. 
     *                     VVHDGV(IV)==1 ) THEN 
                         EF = EF * PCLEFT(MY,ICY,IP,IV) 
                       ENDIF 
                    ENDIF 
C 
C 
                    EF = EF * TPHR / MPHR                                ! Starts in g/mi units 
C 
                    EF = EF * INHWT(INH) * FDSWT(IFDS) * STDWT(ISTD)     ! Starts in g/mi units 
C 
                    HEF(IH,5) = HEF(IH,5) + EF                           ! Hourly Start emissions (g/mi) 
C 
                    STARTSUM = STARTSUM + EF * HVMT(IH)                  ! Daily start in g/mi units 
C 
C 
C ..  Array HEF(IH,IROAD,IU) contains the start emissions in cells where IROAD = 5. 
C     Subscript IU is the units subscript.  IU = 1 is for g/mi and IU = 2 is g/hour. 
C     IH is from 1-24.  1-24 correspond to the hour of the day. 
C     STARTSUM is the daily average start emission value in g/mi. 
C 
C 
C 
C  *****************   RUNNING EMISSIONS SECTION  ******************** 
C 
                  ELSE IF(ISR.EQ.RUN) THEN        !FOR RUNNING EMISSIONS 
C 
C 
C     Calculate and apply Speed Correction Factors, apply the 
C     Temperature Correction Factor and then calculate the full Air 
C     Conditioning Effect. 
C 
                    EF = BASEF(RUN) 
                    IF(SPD_FLAG.EQ.1)THEN 
                      IF(IV.EQ.VTMC) CALL SCFTABMC(IP,IV,MY,IDX,RUN) 
                      IF(VVLDGAS(IV).EQ.1.AND.SCFUSER.EQ.1) 
     *                CALL SCFSLOW(IP,IV,MY) 
                    END IF 
C 
                    IF(SPD_FLAG.EQ.1) THEN 
                      EF = SPEED_CF(IH,IV,INH,IP,MY,IDX,ISTD,IROAD, 
     &                      EF,AC_DUM,TCFRUN) 
                    ELSE 
                      AC_DUM = 0.0 
                      EF = EF * TCFRUN 
                    ENDIF 
                    TEF = MAX(EF,0.0) 
C 
C     Speed Correction for NGV.  A/C effects are not used. 
C 
                    IF(SPD_FLAG.EQ.1) THEN 
                      NGV2 = SPEED_CF(IH,IV,INH,IP,MY,IDX,ISTD, 
     *                       IROAD,NGV1(RUN,INH),AC_DUMM,TCFRUN) 
                    ELSE 
                      NGV2 = NGV1(RUN,INH) * TCFRUN 
                    ENDIF 
                    NGV2 = MAX(NGV2,0.0) 
C 
C     Calculates EFs for A/C correction.  These EF's must be 
C     corrected for Temperature and Speed Prior to A/C correction.  This 
C     is done in routine AC_ADJ(). 
C     A total of 14 speed combinations are needed for each roadway type. 
C 
C     AC Correction: 
C 
                    IF(VVLDGAS(IV).EQ.1.AND.AC_FLAG.EQ.1) THEN 
                      AC_CORR = AC_ADJ(IV,IH,ACTUAL_MY,ICY,AC_DUM) 
                    ELSE 
                      AC_CORR = 0.0 
                      AC_DUM = 0.0 
                    ENDIF 
C 
C     Apply the NOx humidity correction factor to NGVs. 
C 
                    NGV2 = NGV2 * NOXHH 
C 
                    TEF = TEF + AC_CORR + DDNOX(IROAD)*CONV 
C 
                    TEF = TEF * FUELCF(RUN,INH,ISTD,IFDS,IH) 
C 
                    TEF = NGV2*NGVFRAC + TEF*(1.0-NGVFRAC) 
C 
C     For PRE-1981 gas vehicles for all model years obtain % of emission factor 
C     left after I/M program effect is computed for the IH loop. Post-1980 vehicles 
C     are given a I&M credit in grams/mile to high emission vehicles only. Default 
C     (non-I/M, diesel, HD & MC) is no effect. 
C 
                    IF(IMFLAG>1) THEN 
                      IF(VVLDGAS(IV)==1.AND.MY>=1981.AND.INH==HIGH) THEN 
                        TEF = TEF - 
     *                     IMPROG(ICY,MY,ACTUAL_MY,JDX,IP,IV,VMTAGE,ISR, 
     *                        IFDS,ISTD,INERR,IDX,IH,IROAD,CAT,TEF) 
                      ELSEIF((VVLDGAS(IV)==1 .AND. MY<1981) .OR. 
     *                     VVHDGV(IV)==1 ) THEN 
                        TEF = TEF * PCLEFT(MY,ICY,IP,IV) 
                      ENDIF 
                    ENDIF 
C 
                    TEF = TEF * INHWT(INH) * 
     &                        FDSWT(IFDS) * STDWT(ISTD) 
C 
C 
                    HEF(IH,IROAD) = HEF(IH,IROAD) + TEF                ! Running in g/mi 
C 
                    RUNSUM = RUNSUM + TEF * HVMT(IH)*FVMT(IROAD,IH,IV)  ! Daily Running in g/mi 
C 
                  ENDIF                  !END RUNNING/START ENDIF 
C 
                END DO                   !END INH LOOP 
              END DO                  !END IFDS LOOP 
            END DO                 !END ISTD LOOP 
          END DO                !END IROAD LOOP 
        END DO               !END ISR LOOP 
      END DO              !END IH LOOP 
C 
ccs      IF(DBFLAG.EQ.2) CALL DBEXHRS(HEF,IP,IV,ACTUAL_MY,JDX,INERR) 
C
ccs   add call to new output function for special aggregation
      IF( DBFLAG == 2 ) CALL SMKEXOUT(HEF,IP,IV,JDX)

      RETURN 
      END SUBROUTINE HRLOOP 
