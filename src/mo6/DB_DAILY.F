      SUBROUTINE DB_DAILY(EF, HSRFAC, DBHEADER, RFE, IROAD, IV, JDX, MY) 
C 
C     This subroutine outputs the DAILY evaporative, exhaust and 
C     refueling emissions to the Database. 
C 
C 
C     Called by:  DBEVAP and DBEXHRS 
C 
C     Calls:  None 
C 
C     CHANGES  :  (LAST CHANGE FIRST) 
C 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  13 Feb 01 @EPA-elg increased the size of FVMT arrays to include vehicle class 
C  20 Oct 00 AIR Task 03: Added MYR to output 
C  17 Jul 00 @EPA-ddj Bug187 determination of miles travelled in an hour 
C  10 Jul 00 @EPA-ddj Bug187 Removed 2nd of 8 parameters from DB_DAILY 
C            parameter list. 
C  18 MAY 00 @EPA -DJ; BUG #96, changed output format for ZGMDAY from f9.4 to f10.3 
C  05 MAY 00 @EPA -HMC; CHANGED OUTPUT FORMATS, AND ADDED DAILY OUTPUT 
C 
C    common blocks: 
C 
C    /SPEED9/   FVMT,HVMT 
C    /CEVBMY/   BMYMPD 
C    /MYRCAL/   MEVMYR 
C    /STDIST/   STARTDIST, STARTSPERDAY 
C    /WEEKLY/   IWEEK 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'SPEED9.I'          ! FVMT,HVMT 
      INCLUDE  'CEVBMY.I'          ! BMYMPD 
      INCLUDE  'EVAHS1.I'          ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'MYRCAL.I'          ! MEVMYR 
      INCLUDE  'STDIST.I'          ! STARTDIST, STARTSPERDAY 
      INCLUDE  'VMXCOM.I'          ! DBVCOUNT 
      INCLUDE  'WEEKLY.I'          ! IWEEK 
C 
C Input arguments 
C 
      REAL, INTENT(IN)    :: EF(24) 
      REAL, INTENT(IN)    :: HSRFAC 
      REAL, INTENT(IN)    :: RFE 
C 
      INTEGER, INTENT(IN) :: IROAD 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: MY 
C 
      CHARACTER(LEN=25), INTENT(IN) :: DBHEADER 
C 
C Local variables 
C 
      INTEGER             :: IH       ! HOUR INDEX 
C 
      REAL                :: FACVMT 
      REAL, DIMENSION(24) :: MILES    ! Stores hourly miles 
      REAL                :: ZGMDAY   ! SUM of DAILY GRAMS/DAY 
      REAL                :: ZSTARTS  ! SUM of DAILY STARTS 
      REAL                :: ZENDS    ! SUM of DAILY ENDS 
      REAL                :: ZMILES   ! SUM of DAILY MILES 
      REAL                :: ZGMMILE  ! Stores DAILY GRAMS/MILE 
      REAL                :: FMILES   ! SUM of MILES/DAY FOR AN FTYPE 
C 
C     SUM HOURLY NUMBERS TO DAILY 
C 
      ZGMDAY = 0.0 
      ZSTARTS = 0.0 
      ZENDS  = 0.0 
      ZMILES  = 0.0 
      DO IH = 1,24 
          MILES(IH) = BMYMPD(JDX,IV)*HVMT(IH) 
          ZGMDAY = ZGMDAY + EF(IH)*MILES(IH) 
          ZSTARTS= ZSTARTS + 
     &              STARTSPERDAY(JDX,IV,IWEEK)*STARTDIST(IH,IWEEK) 
          ZENDS =  ZENDS +  STARTSPERDAY(JDX,IV,IWEEK)* 
     &             HSRFAC*HS_HRLY_TRP_FRCN(IH,IWEEK) 
          ZMILES = ZMILES + MILES(IH) 
      END DO 
      IF(ZMILES .LT. 0.000001) THEN 
          ZGMMILE = 0.0 
          FACVMT = 1.0 
      ELSE 
          ZGMMILE = ZGMDAY/ZMILES 
C 
C           CALCULATE DAILY FACVMT 
C 
          IF (IROAD .EQ. 5) THEN 
             FACVMT = 1.0 
          ELSE 
             FMILES = 0.0 
             DO IH = 1,24 
                FMILES = FMILES + FVMT(IROAD,IH,IV) * MILES(IH) 
             END DO 
             FACVMT = FMILES/ZMILES 
          END IF 
      END IF 
C 
C           WRITE DAILY DATABASE OUTPUT 
C 
      WRITE (DBUNIT(1),220) DBHEADER, 
     &       ZGMMILE, TABCHAR, ZGMDAY, TABCHAR, ZSTARTS, TABCHAR, 
     &       ZENDS, TABCHAR, ZMILES, TABCHAR, RFE, TABCHAR, 
     &       FACVMT, TABCHAR, MEVMYR(JDX,IV), TABCHAR, DBVCOUNT(IV), 
     &       TABCHAR, MY 
  220 FORMAT(A25,F9.4,A1,F10.3,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.4,A1,F8.4,A1,F10.4,A1,I4) 
C 
      RETURN 
      END SUBROUTINE DB_DAILY 
