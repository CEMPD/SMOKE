      SUBROUTINE ATHRLOOP(IV,ICY,MY,ACTUAL_MY,IDX,CAT,INERR) 
C 
C  ATHRLOOP calculates Air Toxic VOC and TOG emission factors and 
C  stores the values by hour and road type. 
C 
C  Called by ATCALX 
C 
C  Calls: 
C 
C        AC_ADJ 
C        BIGCFX 
C        BIGIN2 
C        ATSPEED_CF 
C        EF_CASE 
C        EF_HI_FRCN 
C        FUEL 
C        GBHP_GPM 
C        IEVPTR 
C        IVPTRT 
C        IMPROG 
C        METHANE 
C        PCLEFT 
C        SCFSLOW 
C        SCFTABMC 
C        SOAKCF 
C 
C  Changes:  (Last change first) 
C 
C  31 Dec 01 @EPA-BAG bug368 Remove I&M sawtooth 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C  20 Sep 01 AIR Task 15: Comments 
C  09 Jul 01 AIR Task 15: New subroutine (Based on HRLOOP) 
C 
C  Input on call: 
C 
C    Constants: 
C    MAXVEH from include file IVTYPE.I 
C    MAXEMIT, MAXIH from MAXIMA.I 
C 
C    argument list: IV,ICY,MY,IDX,INERR 
C 
C  Output on return: 
C 
C    argument list: RUNSUM,STARTSUM 
C 
C    common blocks: 
C     /ALUHIN/  ABSHUM 
C     /BASEQ9/  MAXLYR 
C     /CEVBMY/  BMYMPD 
C     /CUMCOM/  CUMMIL 
C     /FLAGS1/  SCFUSER 
C     /FLAGS2/  IMFLAG 
C     /FLAGS5/  AC_FLAG, SOAK_FLAG, SPD_FLAG 
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
C                rate that can be tampering excess. 
C  CAT      R    Catalyst fractions 
C  EF       R    Temporary Variable holding the Basic Emission Factor 
C  NGV_EF   R    Temporary Variable holding the NGV Emission Factor 
C  FDSWT    R    Fuel delivery system weighting factor 
C  ATHEF    R    Accumulates EFs for FDS/emitter groups for each MY and VOC/TOG 
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
      USE DATABASE, ONLY: DBFLAG 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'ALUOUT.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'BYMYC2.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
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
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
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
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(INOUT) :: INERR 
      REAL   , INTENT(IN) :: CAT 
C 
C     Declare external functions. 
C 
      INTEGER, EXTERNAL :: IEVPTR 
C 
      REAL, EXTERNAL :: AC_ADJ 
      REAL, EXTERNAL :: ATSPEED_CF 
      REAL, EXTERNAL :: EF_CASE 
      REAL, EXTERNAL :: EF_HI_FRCN 
      REAL, EXTERNAL :: IMPROG 
      REAL, EXTERNAL :: PCLEFT 
      REAL, EXTERNAL :: SOAKCF 
C 
C     Declare local variables. 
C 
      REAL     AC_CORR 
      REAL     AC_DUM 
      REAL     AC_DUMM 
      REAL, DIMENSION(2)           :: ATVOC 
      REAL, DIMENSION(2)           :: ATTOG 
      REAL, DIMENSION(2)           :: BASEF 
      REAL, SAVE :: CAP=0.40 
      REAL     EF 
      REAL, DIMENSION(3)           :: FDSWT 
      REAL, DIMENSION(2,3,2)       :: HIGHER 
      REAL, DIMENSION(2)           :: INHWT 
      REAL     MPHR 
      REAL, DIMENSION(2,2)         :: NGV1 
      REAL     NGV2 
      REAL     NGV_EF 
      REAL     SKCF 
      REAL, DIMENSION(10)          :: STDWT 
      REAL     TAMP 
      REAL     TEF 
      REAL     TPHR 
      REAL     VMTAGE 
      REAL, DIMENSION(2,2,10,3,24,2) :: FUELCF 
C 
      INTEGER  IFDS 
      INTEGER  IG5 
      INTEGER  IH 
      INTEGER  INH 
      INTEGER  IROAD 
      INTEGER  ISR 
      INTEGER  ISTD 
      INTEGER  IVTEMP 
      INTEGER  IVT 
      INTEGER  JDX 
      INTEGER  JSTD 
      INTEGER  LMY 
      INTEGER  NGVPT 
      INTEGER  ST_INDX 
      INTEGER, PARAMETER      :: USE_SOAK=1 
      INTEGER, PARAMETER      :: RUN=1,START=2 
      INTEGER, PARAMETER      :: NORMAL=1,HIGH=2 
      INTEGER, PARAMETER      :: IP=1   !HC 
C 
C  Call BIGIN2 to initialize the pointers, arrays and variables used 
C  by METHANE, BIGCFX and the IH loop below. 
C 
      JDX = (MAXYRS + 1) - IDX 
      CALL BIGIN2(IV,MY,JDX,IDX,IP) 
C 
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
          DO IVT = 1,2 
            ATHEF(IH,IROAD,IVT) = 0.0 
          END DO 
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
C 
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
C     Initialize the composite Air Toxic normal/high emitter weight array 
C 
      ATNHWT = 0.0 
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
                IF(ABS(TAMP).GT.BASEF(ISR)*CAP .AND. 
     *                                 MY.GE.1981) THEN 
                  IF(TAMP.GT.0) THEN 
                     TAMP = BASEF(ISR)*CAP 
                  ELSEIF(TAMP.LE.0) THEN 
                     TAMP = -1.0*BASEF(ISR)*CAP 
                  ENDIF 
                ENDIF 
C 
                BASEF(ISR) = BASEF(ISR) + TAMP 
C 
C     Compute VOC and TOG using Air Toxic version of METHANE routine. 
C 
                ATVOC(ISR) = BASEF(ISR) 
                ATTOG(ISR) = BASEF(ISR) 
C 
                IF(IV.LE.5 .OR. (IV.GE.6.AND.ISTD.EQ.1))  
     *            CALL ATMETHANE(IV,MY,ISTD,IFDS,IDX,ISR,ATVOC,ATTOG) 
C 
C     ATEMFAC now contains the Base VOC and TOG Emission Factors corrected for 
C           Tampering, and Methane. 
C 
                ATEMFAC(ISR,ISTD,IFDS,INH,1) = ATVOC(ISR) 
                ATEMFAC(ISR,ISTD,IFDS,INH,2) = ATTOG(ISR) 
C 
             END DO   !ISR LOOP 
C 
C     Compute the fuel correction factors for VOC and TOG, using AT version 
C     of FUEL. Must call twice to fill FUELCF array because fuel correction 
C     depends on magnitude of input EF (ATVOC or ATTOG) 
C 
             CALL ATFUEL(MY,IV,ICY,IFDS,INH,IDX,ISTD,1,ATVOC,FUELCF) 
             CALL ATFUEL(MY,IV,ICY,IFDS,INH,IDX,ISTD,2,ATTOG,FUELCF) 
C 
C     Compute the composite AT normal/high emitter array 
C 
             ATNHWT(INH) = ATNHWT(INH)  
     *                   + HIGHER(INH,IFDS,JSTD) 
     *                   * FDSWT(IFDS)  
     *                   * STDWT(ISTD) 
C 
          END DO      !ISTD 
        END DO        !INH LOOP 
      END DO          !IFDS 
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
        DO ISR = 1,2  ! 1=Run, 2=Start 
C 
C       Heavy duty vehicles and busses do not have engine start emissions. 
C 
          IF(VVHEAVY(IV).EQ.1 .AND. ISR.EQ.START) CYCLE 
C 
C     Calculation of trips per hour and miles driven in this hour 
C 
          IF(ISR.EQ.START) THEN 
             TPHR = STARTSPERDAY(JDX,IV,IWEEK) * 
     *              STARTDIST(IH,IWEEK)               ! Starts per hour 
             MPHR = BMYMPD(JDX,IV) * HVMT(IH)         ! Miles per hour 
          ENDIF 
C 
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
C     Loop for VOC (IVT=1) and TOG (IVT=2) 
C     (Programming Note: Must put loop here because various correction 
C      factors, (speed, oxy, etc) depend on the magnitude of EF.) 
C 
                DO IVT=1,2 
C 
                  BASEF(ISR) = ATEMFAC(ISR,ISTD,IFDS,INH,IVT) 
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
                    EF = EF * FUELCF(START,INH,ISTD,IFDS,IH,IVT) 
C 
C     Weight the EFs. 
C 
                    EF = NGV_EF*NGVFRAC + EF*(1.0-NGVFRAC) 
C 
C     For gas vehicles for all model years obtain % of emission factor 
C     left after I/M program effect is computed for the IH loop. 
C     Default (non-I/M, diesel, HD & MC) I/M is no effect or no credit. 
C 
                    IF(IMFLAG>1) THEN 
                      IF(VVLDGAS(IV)==1.AND.MY>=1981.AND.INH==HIGH) THEN 
                         EF = EF - 
     *                     IMPROG(ICY,MY,ACTUAL_MY,JDX,IP,IV,VMTAGE,ISR, 
     *                     IFDS,ISTD,INERR,IDX,IH,IROAD,CAT,EF) 
                      ELSEIF((VVLDGAS(IV)==1 .AND. MY<1981) .OR. 
     *                     VVHDGV(IV)==1 ) THEN 
                           EF = EF*PCLEFT(MY,ICY,IP,IV) 
                       ENDIF 
                    ENDIF 
C 
                    EF = EF * TPHR / MPHR                                ! Starts in g/mi units 
C 
                    EF = EF * INHWT(INH) * FDSWT(IFDS) * STDWT(ISTD)     ! Starts in g/mi units 
C 
                    ATHEF(IH,5,IVT) = ATHEF(IH,5,IVT) + EF               ! Hourly Start emissions (g/mi) 
C 
C ..  Array ATHEF(IH,IROAD,IVT) contains the start emissions, g/mi,  
C     in cells where IROAD = 5. IH is from 1-24.  1-24 correspond  
C     to the hour of the day. IVT is 1 for VOC, 2 for TOG 
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
                      EF = ATSPEED_CF(IH,IV,INH,IP,MY,ISTD,IROAD, 
     &                      IVT,EF,AC_DUM,TCFRUN) 
                    ELSE 
                      AC_DUM = 0.0 
                      EF = EF * TCFRUN 
                    ENDIF 
                    TEF = MAX(EF,0.0) 
C 
C     Speed Correction for NGV.  A/C effects are not used. 
C 
                    IF(SPD_FLAG.EQ.1) THEN 
                      NGV2 = ATSPEED_CF(IH,IV,INH,IP,MY,ISTD, 
     *                       IROAD,IVT,NGV1(RUN,INH),AC_DUMM,TCFRUN) 
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
C     Apply the air conditioning and fuel corrections. 
C 
                    TEF = TEF + AC_CORR 
C 
                    TEF = TEF * FUELCF(RUN,INH,ISTD,IFDS,IH,IVT) 
C 
C     Weight the EFs. 
C 
                    TEF = NGV2*NGVFRAC + TEF*(1.0-NGVFRAC) 
C 
C     For gas vehicles for all model years obtain % of emission factor 
C     left after I/M program effect is computed for the IH loop. 
C     Default (non-I/M, diesel, HD & MC) I/M is no effect or no credit. 
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
                    ATHEF(IH,IROAD,IVT) = ATHEF(IH,IROAD,IVT) + TEF     ! Running in g/mi 
C 
                  ENDIF                  !END RUNNING/START ENDIF 
C 
                  END DO                   !END IVT (VOC and TOG) LOOP 
                END DO                   !END INH LOOP 
              END DO                  !END IFDS LOOP 
            END DO                 !END ISTD LOOP 
          END DO                !END IROAD LOOP 
        END DO               !END ISR LOOP 
      END DO              !END IH LOOP 
C 
      RETURN 
      END SUBROUTINE ATHRLOOP 
