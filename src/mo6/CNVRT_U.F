      REAL FUNCTION 
     *CNVRT_U(IN_UNIT,OUT_UNIT,JDX,ICY,IV,IVTAM,IH,EMISSION_TYPE) 
C 
C  CNVRT_U does unit conversions for the evaporative emissions. 
C 
C  Called by EVP_HRLY 
C 
C  Calls 
C 
C  Changes:  (Last change first) 
C 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  10 Mar 00 @DynTel-mla  1-011 moved HSRED fractions to block data 
C  28 fEB 00 @DynTel-bg  1-006 Changed conversion factor to 1.0 for daily 
C            gram/mile average to hourly gram/mile  
C  12 Jan 00 @DynTel-bg  1-006 Renamed to cnvrt_u.for 
C 
C    common blocks: 
C    /CEVBMY/  BMYMPD,BMTPDF 
C    /EVAHS1/  HS_HRLY_TRP_FRCN, HSRED 
C    /SPEED9/  HVMT 
C    /STDIST/  STARTSPERDAY(28,2) 
C    /RLCOM1/  ROADFE 
C    /WEEKLY/  IWEEK 
C 
      IMPLICIT NONE 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'EVAHS1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'RLCOM1.I' 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'STDIST.I' 
      INCLUDE 'WEEKLY.I' 
      INCLUDE 'MYRCAL.I' 
CCC      INCLUDE 'UNIT_PARAM.I' 
C 
C     declare parameter list 
C 
      INTEGER, INTENT(IN) :: EMISSION_TYPE 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IVTAM 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: OUT_UNIT 
      INTEGER, INTENT(IN) :: IN_UNIT 
C 
C     declare local variables 
C 
      INTEGER  MDX 
      INTEGER  MY 
      REAL     TPD 
      INTEGER, PARAMETER :: IHOTSOAK=1 
      INTEGER, PARAMETER :: GPT=1,GPH=2,GPG=3,GPM_D=4,GPM=5 
      INTEGER, PARAMETER :: CAR=1,TRUCK=2 
C 
C     Set miles per day and trips per day 
C     according to vehicle class, age, week or weekend. 
C 
      TPD = STARTSPERDAY(JDX,IV,IWEEK) 
C 
C     Find units conversion factors for hourly evaporative emissions 
C 
      SELECT CASE (IN_UNIT) 
C 
        CASE (GPT) 
C 
C         The number of hot soaks per day are reduced from the number of 
C         trips per day by vehicle class (see p.5  of M6.FLT.004) 
          IF (EMISSION_TYPE.EQ.IHOTSOAK .AND. IV.EQ.VTLDGV) 
     *                                            TPD=TPD*HSRED(CAR) 
          IF (EMISSION_TYPE.EQ.IHOTSOAK .AND. IV.NE.VTLDGV) 
     *                                            TPD=TPD*HSRED(TRUCK) 
C 
          IF(OUT_UNIT .EQ. GPM) CNVRT_U= 
     *       (HS_HRLY_TRP_FRCN(IH,IWEEK)*TPD)/(HVMT(IH)*BMYMPD(JDX,IV)) 
C 
          IF(OUT_UNIT .EQ. GPH) CNVRT_U= 
     *                               HS_HRLY_TRP_FRCN(IH,IWEEK)*TPD 
C 
        CASE (GPH) 
C 
          IF(OUT_UNIT .EQ. GPM ) CNVRT_U = 
     *                               1.0/(BMYMPD(JDX,IV)*HVMT(IH)) 
C 
        CASE (GPM) 
C 
          IF(OUT_UNIT .EQ. GPH) CNVRT_U = HVMT(IH)*BMYMPD(JDX,IV) 
C 
          IF(OUT_UNIT .EQ. GPG) THEN 
C 
C           Set the road fuel economy MY pointer. 
            MY=ICY-JDX+1 
            IF(MY.LE.1969)MDX=1 
            IF(MY.GT.1969.AND.MY.LT.2000)MDX=MY-1969+1 
            IF(MY.GE.2000)MDX=32 
C 
            CNVRT_U = ROADFE(MDX,IVTAM) 
C 
          ENDIF 
C 
        CASE (GPM_D) !convert from a daily average (over 24 hours) 
C                     to the individual hourly numbers 
C 
          IF(OUT_UNIT .EQ. GPM) CNVRT_U = 1.0 
C 
        CASE DEFAULT 
          CNVRT_U =0.0 
C 
      END SELECT 
C 
      RETURN 
      END FUNCTION CNVRT_U 
