      SUBROUTINE CO2HRLOOP(IV,MY,IDX,INERR)
C
C  CO2LOOP calculates the THC and CO exhaust emission factors that are
C  used to compute CO2. These values are computed only for FTP conditions.
C  No adjustemts are made for I/M, tampering, fuels, temperature, speed,
C  air conditioning, offsets, etc. Soak corrections are applied, but only
C  for the standard 68-84F condition.
C
C  Called by CO2CALX
C
C  Calls:
C
C        BIGIN2
C        EF_CASE
C        EF_HI_FRCN
C        IEVPTR
C        IVPTRT
C        CO2METHANE
C        P3WPOX
C        SOAKCF
C
C  Changes:  (Last change first)
C  11 Jul 02  AIR Task 36: NOx humidity correction made EMFAC array hourly
C             dependent, modified code accordingly.
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  20 Mar 02: Added support for NO HIGH EMITTERS.
C  04 Jan 01: Removed unused code, updated comments
C  13 Dec 01: Expanded NGV trap to include NGVFLG1 and NGVFLG2.
C  06 Dec 01 Applied 17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C  13 Nov 01 AIR Task 27: Special version of HRLOOP (10Jul01) for CO2.
C
C  Input on call:
C
C    Constants:
C    MAXVEH from include file IVTYPE.I
C    MAXEMIT, MAXIH from MAXIMA.I
C
C    argument list: IV,MY,IDX,INERR
C
C  Output on return:
C
C    argument list: RUNSUM,STARTSUM
C
C    common blocks:
C     /BASEQ5/  EMFAC, HIGHER
C     /BASEQ9/  MAXLYR
C     /CEVBMY/  BMYMPD
C     /CUMCOM/  CUMMIL
C     /FLAGS5/  AC_FLAG, SOAK_FLAG, SPD_FLAG, NOHIGH
C     /INJECT/  TBI, PFI
C     /IVPCOM/  IVPTRS, IVPTRT
C     /IVTYPE/  VTGASBUS, VVHDGV, VVLDGAS
C     /MAXIMA/  MAXIH
C     /NGV/     NGVFLG1,NGVFLG2
C     /OMTCOM/  TCFSTART,TCFRUN,TCF12
C     /SPEED9/  FVMT, HVMT
C     /STDIST/  STARTDIST, STARTSPERDAY
C     /TEMPS/   TEMHRLY
C     /WEEKLY/  IWEEK
C
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BASEF    R    The basic EF used to calculate the hourly EF
C  CAT      R    Catalyst fractions
C  EF       R    Temporary Variable holding the Basic Emission Factor
C  FDSWT    R    Fuel delivery system weighting factor
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
C  P3W      R    Percent 3-way catalysts
C  POX      R    Percent oxidation catalysts
C  SOAKCF   R    Final soak time correction factors.
C  SKCF     R    Temperary variable to hold the soak correctin factor.
C  STDWT    R    Certification standards weighting factor
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
      INCLUDE 'BASEQ5.I'
      INCLUDE 'BASEQ9.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'CEVBMY.I'
      INCLUDE 'CO2COM.I'
      INCLUDE 'CUMCOM.I'
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
      INCLUDE 'OMTCOM.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'RFORM2.I'
      INCLUDE 'SOAK.I'
      INCLUDE 'SPEED9.I'
      INCLUDE 'STDIST.I'
      INCLUDE 'TEMPC4.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'WEEKLY.I'
C
C     Declare parameters.
C
      INTEGER, INTENT(IN) :: IDX
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(INOUT) :: INERR
C
C     Declare external functions.
C
      INTEGER, EXTERNAL :: IEVPTR
C
      REAL, EXTERNAL :: EF_CASE
      REAL, EXTERNAL :: EF_HI_FRCN
      REAL, EXTERNAL :: HDCPIC
      REAL, EXTERNAL :: SOAKCF
C
C     Declare local variables.
C
      REAL, DIMENSION(2)           :: BASEF
      REAL     CAT
      REAL     EF
      REAL, DIMENSION(3)           :: FDSWT
      REAL, DIMENSION(2,3,2)       :: HIGHER
      REAL, DIMENSION(2)           :: INHWT
      REAL     MPHR
      REAL     P3W
      REAL     POX
      REAL     SKCF
      REAL, DIMENSION(10)          :: STDWT
      REAL     TPHR
      REAL     VMTAGE
C
      INTEGER  IP
      INTEGER  IFDS
      INTEGER  IG5
      INTEGER  IH
      INTEGER  INH
      INTEGER  ISR
      INTEGER  ISTD
      INTEGER  IVTEMP
      INTEGER  JDX
      INTEGER  JSTD
      INTEGER  LMY
      INTEGER  ST_INDX
      INTEGER, PARAMETER      :: USE_SOAK=1
      INTEGER, PARAMETER      :: RUN=1,START=2
      INTEGER, PARAMETER      :: NORMAL=1,HIGH=2
C
C  Since CO2 data for NGV is sparse, do not allow CO2 with NGV
C  to run.
C
      IF(NGVFLG1.EQ.1 .OR. NGVFLG2.EQ.1) THEN
      CALL QUITER(0.,0,83,INERR)
      RETURN
      END IF
C
C  Initialize the THC and CO holding array for CO2.
C
      CO2HEF=0.0
C
C  Loop for HC (THC) and CO only.
C
      DO IP=1,2
C
C  Call the soak function, SOAKFCN, to set up the array, SOAKFUN, for each IP.
C  These values are used to determine hourly soak time correction factors.
C  They are dependent only on soak time, pollutant, and type of catalyst.
C
      IF(SOAK_FLAG .EQ. 1) CALL SOAKFCN(IP)
C
C  Call BIGIN2 to initialize the pointers, arrays and variables used
C  by CO2METHANE and the IH loop below.
C
      JDX = (MAXYRS + 1) - IDX
      CALL BIGIN2(IV,MY,JDX,IDX,IP,INERR)
      IF(INERR.GT.0) GOTO 99
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
        IG5 = IEVPTR(MY,IVTEMP,2,INERR)
        IF(INERR.GT.0) GOTO 99
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
            IF (VVLDGAS(IV).EQ.0.OR.MY.LT.1981.OR.NOHIGH.EQ.1) THEN
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
     &         (VVLDGAS(IV).NE.1 .OR.
     &                   MY.LT.1981 .OR. NOHIGH.EQ.1)) CYCLE
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
     *                                   VMTAGE,ISTD,INH,IFDS,INERR)
                    IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0
                  ELSE
                    IF(ISTD.GE.1 .AND. MY.GE.1994) THEN
                      BASEF(ISR) = EF_CASE(MY,JDX,IP,IV,ISR,
     *                                     VMTAGE,ISTD,INH,IFDS,INERR)
                      IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0
                    ELSE
                      CYCLE
                    ENDIF
                  ENDIF
                ELSE
                  BASEF(ISR) = EF_CASE(MY,JDX,IP,IV,ISR,
     *                                 VMTAGE,ISTD,INH,IFDS,INERR)
                  IF(BASEF(ISR).LT.0.00001) BASEF(ISR) = 0.0
                ENDIF
                IF(INERR.GT.0) GOTO 99
C
C     CO2METHANE subroutine returns the basic emission factor in THC
C
                IF(IP.EQ.1) THEN
                  IF(IV.LE.5) THEN
                    CALL CO2METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,BASEF)
                  ELSEIF(IV.GE.6 .AND. ISTD.EQ.1) THEN
                    CALL CO2METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,BASEF)
                  ENDIF
                ENDIF
C
C     EMFAC now contains the Base Emission Factor corrected for Methane.
C     (EMFAC array is used via common, so store hourly values)
C
                DO IH=1,MAXIH
                  EMFAC(ISR,ISTD,IFDS,INH,IH) = BASEF(ISR)
                END DO
             END DO   !ISR LOOP
C
          END DO      !ISTD
        END DO        !INH LOOP
      END DO          !IFDS
C
C     Begin Main Looping Structure of HRLOOP.
C
      DO IH = 1,MAXIH
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
                  BASEF(ISR) = EMFAC(ISR,ISTD,IFDS,INH,IH) ! Hourly dependent
C
C  *******************   START EMISSIONS SECTION  ********************
C
               IF(ISR.EQ.START) THEN !FOR START EMISSIONS
C
C     Apply start emissions soak corrections, as required. If no saok is
C     specified, then use a temperature correction corresponding to a 12-hour soak. 
C     NOTE: Soak corrections are standard 68-84F only.
C       
                   IF(SOAK_FLAG.EQ.USE_SOAK) THEN
                     EF=0.0
                     DO ST_INDX=1,69
                        SKCF = SOAKCF(MY,CAT,IH,IWEEK,ST_INDX)
                        EF = EF + BASEF(START) * SKCF
                     END DO
                    ELSE
                     EF = BASEF(START) * TCF12
                    ENDIF
C
                    EF = EF * TPHR / MPHR                                ! Starts in g/mi units
C
                    EF = EF * INHWT(INH) * FDSWT(IFDS) * STDWT(ISTD)     ! Starts in g/mi units
C
C     Add start emissions to CO2 holding array. Note that the values
C     are weighted over the entire 24-hours.
C
                    CO2HEF(IP) = CO2HEF(IP) + EF * HVMT(IH)              ! Hourly Start emissions (g/mi)
C
C  *****************   RUNNING EMISSIONS SECTION  ********************
C
                  ELSE IF(ISR.EQ.RUN) THEN        !FOR RUNNING EMISSIONS
C
                    EF = BASEF(RUN)
C
                    EF = EF * INHWT(INH) * FDSWT(IFDS) * STDWT(ISTD)
C
C     Add start emissions to CO2 holding array.  Note that the values
C     are weighted over the entire 24-hours.
C
                    CO2HEF(IP) = CO2HEF(IP) + EF * HVMT(IH)    ! Running in g/mi
C
                  ENDIF                  !END RUNNING/START ENDIF
C
                END DO                   !END INH LOOP
              END DO                  !END IFDS LOOP
            END DO                 !END ISTD LOOP
        END DO               !END ISR LOOP
       END DO              !END IH LOOP
      END DO             !END IP LOOP
C
   99 RETURN
      END SUBROUTINE CO2HRLOOP
