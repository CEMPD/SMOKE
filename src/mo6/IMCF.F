      SUBROUTINE IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,INH,
     *                IH,IROAD,CAT,EF,INERR)
C
C  IMCF calculates and applies all of the correction factors found
C  in HRLOOP to the intermediate I/M Emission levels.
C
C  Called by IMPROG
C
C  Calls
C
C  Changes: (Last change first)
C
C  08 Jul 02 AIR Task 36: added hourly absolute humidity processing to
C            NOx humidity corrections.
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  04 Dec 00 @EPA-bag Task X6, added model year, MY, to SOAKCF's parameter list;
C            removed MY, JDX, and IDX from MYCODE.I; added ACTUAL_MY
C  29 Sep 00 AIR Task 07: Bug Fix 271: Fixed logic with use of variables
C            NGV_EF and NGV2
C  29 Aug 00 @EPA-elg created new subroutine.
C
C  Input on call:
C
C    argument list: ICY,IV,IP,ISR,IDX,MY,ISTD,IFDS,INH,IH,IROAD,CAT,EF
C
C    common blocks:
C    /ALUHIN/  ABSHUM
C    /FLAGS1/   SCFUSER
C    /FLAGS5/   AC_FLAG, SOAK_FLAG, SPD_FLAG
C    /IMBLOCK1/ FUELCF
C    /IVTYPE/   VTGASBUS, VVHDGV, VVLDGAS
C    /NGV/      NGVFLG1, NATGAS, NGVFRAC
C    /OFFSET/   OFFSET
C    /OMTCOM/   TCFSTART, TCFRUN, TCF12
C    /TAMOUT/   TEX
C    /WEEKLY/   IWEEK
C
C  Output on return:
C
C
C
C  Notes:
C
      IMPLICIT NONE
C
      INCLUDE 'ALUHIN.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS5.I'
      INCLUDE 'IMBLOCK.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'NGV.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'OMTCOM.I'
      INCLUDE 'TAMOUT.I'
      INCLUDE 'WEEKLY.I'
C
C     Declare external functions.
C
      REAL, EXTERNAL :: AC_ADJ
      REAL, EXTERNAL :: SPEED_CF
      REAL, EXTERNAL :: SOAKCF
      REAL, EXTERNAL :: CABSHUM
C
      REAL  EF,CAT
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: ACTUAL_MY
      INTEGER, INTENT(IN) :: ICY
      INTEGER, INTENT(IN) :: IDX
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: IROAD
      INTEGER, INTENT(IN) :: ISR
      INTEGER, INTENT(IN) :: IFDS
      INTEGER, INTENT(IN) :: ISTD
      INTEGER, INTENT(IN) :: IH
      INTEGER, INTENT(INOUT) :: INERR
C
      INTEGER INH
C
C     Declare local variables.
C
      REAL  AC_DUM, AC_DUMM, AC_CORR
      REAL, DIMENSION(2)           :: EF_A
      REAL, DIMENSION(3), SAVE     :: CAP=(/0.40,0.15,0.20/)
      REAL  EFSUM, NGV_EF, SKCF
      REAL  HUMIDITY
      REAL  NGV2
      REAL  NOXHH
      REAL  TAMP
      REAL  HRABSHUM
C
      INTEGER ST_INDX
C
C
C     Add Tampering Offset to the basic emission factors.
C
C     First compute Tampering Offset.
C
      TAMP = TEX(IV,IP,ISR,IDX)
C
      IF(ABS(TAMP).GT.EF*CAP(IP) .AND. MY.GE.1981) THEN
        IF(TAMP.GT.0) THEN
          TAMP = EF*CAP(IP)
        ELSEIF(TAMP.LE.0) THEN
          TAMP = -1.0*EF*CAP(IP)
        ENDIF
      ENDIF
C
      EF = EF + TAMP
C
C
C     METHANE subroutine returns the basic emission factor in the HC
C             speciation units requested by the user.
C
      IF(IP.EQ.1) THEN
        EF_A(1) = EF
        EF_A(2) = EF
        IF(IV.LE.5) THEN
          CALL METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,EF_A)
        ELSEIF(IV.GE.6 .AND. ISTD.EQ.1) THEN
          CALL METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,EF_A)
        ENDIF
        EF = EF_A(ISR)
      ENDIF
C
C
C     NOx humidity Correction carried over from MOBILE5.
C     Use actual abshum whenever user hourly rel hum is available.
C     (User-supplied Rel Hum data supercedes user-supplied abshum entry.)
C
      IF((VVLDGAS(IV).EQ.1. OR. IV.EQ.VTMC) .AND. IP.EQ.3) THEN
            IF(RH_DIST) THEN
              HRABSHUM=CABSHUM(IH) ! Use user-supplied data, when available
            ELSE
              HRABSHUM=ABSHUM
            END IF
            HUMIDITY=MAX(HRABSHUM,21.)    !humidity test data lower limit
            HUMIDITY=MIN(HRABSHUM,124.)   !humidity test data upper limit
          NOXHH = 1.0 - (HUMIDITY - 75.0)/263.2
          NOXHH = MAX(NOXHH, 0.0)
      ELSE
          NOXHH = 1.0
      END IF
C
      EF = EF * NOXHH
C
C     NGV Emission Calculation
C
      IF(NGVFLG1.EQ.1) THEN
        IF(NGVFRAC .GT.0.0) THEN
           NGV2 = NATGAS(IV,IP,IDX,ISR,INH)
        ELSE
           NGV2 = 0.0
        ENDIF
      ELSE
        NGV2 = 0.0
      ENDIF
C
C
      IF(ISR.EQ.2) THEN !FOR START EMISSIONS
C
C     Apply start emissions temperature, offset and soak corrections,
C     as required. If no saok is specified, then use a temperature 
C     correction corresponding to a 12-hour soak. Thereafter, apply
C     the fuel and other corrections.
C       
        IF(SOAK_FLAG.EQ.1) THEN
          EFSUM=0.0
          NGV_EF=0.0
          DO ST_INDX=1,69
             SKCF = SOAKCF(MY,CAT,IH,IWEEK,ST_INDX)
             EFSUM = EFSUM + (EF * TCFSTART(ST_INDX)
     *          + OFFSET) * SKCF
             NGV_EF = NGV_EF + NGV2 * SKCF
          END DO
        ELSE
          EFSUM = EF * TCF12 + OFFSET
          NGV_EF = NGV2
        ENDIF
        EF = EFSUM
C
      ELSE      !RUNNING ONLY CALCULATIONS
C
C     Calculate and apply Speed Correction Factors, apply the
C     Temperature Correction Factor and then calculate the full Air
C     Conditioning Effect.
C
        EFSUM = EF
        IF(SPD_FLAG.EQ.1) THEN
          IF(IV.EQ.VTMC) CALL SCFTABMC(IP,IV,MY,IDX,ISR,INERR)
          IF(VVLDGAS(IV).EQ.1.AND.SCFUSER.EQ.1) 
     &       CALL SCFSLOW(IP,IV,MY,INERR)
C
          EFSUM = SPEED_CF(IH,IV,INH,IP,MY,IDX,ISTD,IROAD,
     &                     EFSUM,AC_DUM,TCFRUN)
          NGV2 = SPEED_CF(IH,IV,INH,IP,MY,IDX,ISTD,
     *                     IROAD,NGV2,AC_DUMM,TCFRUN)
        ELSE
          AC_DUM = 0.0
          EFSUM = EFSUM * TCFRUN
          NGV2  = NGV2  * TCFRUN
        ENDIF
        EF = MAX(EFSUM,0.0)
        NGV_EF = MAX(NGV2,0.0)
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
        EF = EF + AC_CORR
C
      ENDIF
C
C     Apply all fuel correction factors.
C
      EF = EF * FUELCF(ISR,INH,ISTD,IFDS,IH)
C
C     Apply the NOx humidity correction factor to NGVs.
C
      NGV_EF = NGV_EF * NOXHH
C
      EF = NGV_EF*NGVFRAC + EF*(1.0-NGVFRAC)
C
C
      RETURN
      END
