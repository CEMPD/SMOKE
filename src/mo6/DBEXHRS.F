      SUBROUTINE DBEXHRS(HEF,IP,IV,MY,JDX,INERR) 
C 
C  DBEXHRS writes Exhaust Start and Exhaust Running emissions data to 
C  the EMISSIONS database table. 
C 
C  Called by:  HRLOOP 
C 
C  CallS:  DB_DAILY 
C 
C  Changes: (Last change first) 
C 
C  15 Nov 01 AIR Task 27: Added INERR support due to PMFLECON update. 
C  06 Nov 01 AIR Task 27: Updated fuel economy calculations via PMFLECON. 
C  19 Oct 01 AIR Task: Put in code (again) so that all vehicles have MPG output. 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  13 Feb 01 @EPA-elg increased the size of FVMT arrays to include vehicle class 
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection 
C 19 Oct 00 @EPA-djb Task W2, Added DBAGGR condition on output. 
C 20 Oct 00 AIR Task 03: Added MYR to daily and hourly output. 
C  06 Oct 00 @EPA-bag BUG155 Remove FTP_FLAG 
C 10 Jul 00 @EPA-ddj Bug187 Removed 2nd of 8 parameters (GRAMS_HR) 
C           from calls to DB_DAILY. 
C 30 May 00 @EPA-elg Replaced BDSPEED1.I with SPEED9.I 
C  05 MAY 00 @EPA        E22    ADDED DAILY DATABASE OUTPUT 
C  10 Mar 00 @DynTel-MLA 1-013  database code 
C 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: HEF, IP, IV, JDX, MY 
C 
C    common blocks: 
C 
C    /SPEED9/   FVMT,HVMT 
C    /CEVBMY/   BMYMPD 
C    /FLAGS5/   FTP_INPUT_FLAG 
C    /IVTYPE/   VVSTARTS 
C    /MYRCAL/   MEVMYR 
C    /STDIST/   STARTDIST, STARTSPERDAY 
C    /TEMPS/    TEMHRLY 
C    /VMXCOM/   DBVCOUNT 
C    /WEEKLY/   IWEEK 
C 
C  Output on return: 
C 
C    common blocks: 
C 
C  Local array subscripts: 
C 
C  Local variable dictionary: 
C 
C  Name     Type              Description 
C  ------   ----  ------------------------------------------------------- 
C  DBHEADER  C    Non-changing part of database output record. 
C  IFAC      I    Facility type index. 
C  IH        I    Hourly interval. 
C  IVC       I    Model year group index for vehicle counts. 
C  IW        I    Weekday/weekend flag. 
C  RFE       R    Fuel economy value. 
C  ZGMDAY    R    USED TO SUM DAILY GRAMS/DAY 
C  ZSTARTS   R    USED TO SUM DAILY STARTS 
C  ZENDS     R    USED TO SUM DAILY ENDS 
C  ZMILES    R    USED TO SUM DAILY MILES 
C  ZGMMILE   R    USED TO CALCULATE DAILY GRAMS/MILE 
C  FMILES    R    USED TO SUM MILES/DAY FOR AN FTYPE 
C 
C  Notes: 
C 
C 
      USE DATABASE, ONLY : DBFLAG,DBSELPOL,DBSELVEH,DBSELEFT,DBSELFAC, 
     &                     DBNFILE,DBNRUN,DBNSCEN,DBDAILY,DBAGGR, 
     &                     DBUNIT,TABCHAR,DBSELAGE,DBSELHR 
C 
      IMPLICIT NONE 
      INCLUDE  'SPEED9.I'     ! FVMT,HVMT 
      INCLUDE  'CEVBMY.I'     ! BMYMPD 
      INCLUDE  'EVAHS1.I'     ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'FLAGS5.I'     ! FTP_INPUT_FLAG 
      INCLUDE  'IVTYPE.I'     ! VVSTARTS 
      INCLUDE  'MAXIMA.I'     ! MAXYRS 
      INCLUDE  'MYRCAL.I'     ! MEVMYR,TF 
      INCLUDE  'STDIST.I'     ! STARTDIST, STARTSPERDAY 
      INCLUDE  'TEMPS.I'      ! TEMHRLY 
      INCLUDE  'VMXCOM.I'     ! DBVCOUNT 
      INCLUDE  'WEEKLY.I'     ! IWEEK 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
      REAL,DIMENSION(24,5),INTENT(IN)  ::  HEF 
      INTEGER,INTENT(IN)               ::  IP 
      INTEGER,INTENT(IN)               ::  IV 
      INTEGER,INTENT(IN)               ::  JDX 
      INTEGER,INTENT(IN)               ::  MY 
      INTEGER,INTENT(INOUT)            ::  INERR 
C 
      CHARACTER(25)                ::  DBHEADER 
C 
      INTEGER                      ::  EFTYPE 
      INTEGER                      ::  IFAC 
      INTEGER                      ::  IH 
      INTEGER                      ::  IW 
      INTEGER, SAVE                ::  UPPER_LMT=2 
      INTEGER, SAVE                ::  LOWER_LMT=1 
      LOGICAL                      ::  AGE_OK   ! FLAG FOR AGES/MODEL YEARS TO PROCESS 
C 
      REAL                         ::  HSRFAC 
      REAL                         ::  RFE 
C 
C     Write the output only if database output is selected for the current 
C     run and the database selection criteria include the current pollutant, 
C     the current vehicle type, and the Exhaust Start or Running emission 
C     type. 
C 
      AGE_OK=(DBSELAGE(LOWER_LMT) <= (JDX-1) .AND. 
     &           (JDX-1) <= DBSELAGE(UPPER_LMT)) 
      IF (DBFLAG /= 2 .OR. DBSELPOL(IP) /= 2 .OR. 
     &    DBSELVEH(IV) /= 2 .OR. .NOT.AGE_OK) THEN 
         GOTO 99 
      ENDIF 
C 
C     Don't write the database output if either FTP output 
C     option was selected. 
C 
      IF (FTP_INPUT_FLG /= 0) THEN 
         GOTO 99 
      ENDIF 
C 
      IW = IWEEK 
C 
      RFE = PMFLECON(MY,IV,INERR) ! Compute fuel economy for all vehicle classes 
C 
      IF (IV == VTLDGV) THEN 
        HSRFAC = HSRED(1) 
      ELSE 
        HSRFAC = HSRED(2) 
      ENDIF 
C 
C     Exhaust Running emissions. 
C 
      IF (DBSELEFT(1) == 2) THEN 
C 
         EFTYPE = 1 
C 
         DO IFAC = 1,4 
C 
            IF (DBSELFAC(IFAC) /= 2) CYCLE 
C 
C           WRITE THE HEADER 
C 
            WRITE(DBHEADER,110) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &        DBNSCEN,TABCHAR,IP,TABCHAR,IV,TABCHAR,EFTYPE,TABCHAR, 
     &        IFAC,TABCHAR,JDX-1,TABCHAR 
  110 FORMAT(I3,A1,I3,A1,I3,A1,I1,A1,I2,A1,I2,A1,I1,A1,I2,A1) 
C 
            IF (DBDAILY) THEN 
C 
                CALL DB_DAILY(HEF(1:24,IFAC),HSRFAC,DBHEADER, 
     &                        RFE,IFAC,IV,JDX,MY) 
C 
            ELSEIF (.NOT.DBAGGR) THEN 
C 
C               Write the remaining database fields for hourly output. (The 
C               remaining fields are the hour, gm/mi, gm/hr, Starts, Ends, 
C               Miles, MPG, HrVMT, FacVMT, RegDist, VCount, AmbTemp, and 
C               DiurTemp, and MY. 
C 
                WRITE(DBUNIT(1),115) (DBHEADER(1:25),IH,TABCHAR, 
     &            HEF(IH,IFAC),TABCHAR, 
     &            HEF(IH,IFAC)*BMYMPD(JDX,IV)*HVMT(IH),TABCHAR, 
     &            STARTSPERDAY(JDX,IV,IW)*STARTDIST(IH,IW),TABCHAR, 
     &            STARTSPERDAY(JDX,IV,IW)* 
     &            HSRFAC*HS_HRLY_TRP_FRCN(IH,IW), 
     &            TABCHAR,BMYMPD(JDX,IV)*HVMT(IH),TABCHAR,RFE,TABCHAR, 
     &            HVMT(IH),TABCHAR,FVMT(IFAC,IH,IV),TABCHAR, 
     &            MEVMYR(JDX,IV), 
     &            TABCHAR,DBVCOUNT(IV),TABCHAR,TEMHRLY(IH),TABCHAR, 
     &            TEMHRLY(IH),TABCHAR,MY,IH=DBSELHR(LOWER_LMT), 
     &            DBSELHR(UPPER_LMT)) 
C 
  115 FORMAT(A25,I2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F10.4,A1,F6.1,A1,F6.1, 
     &       A1,I4) 
C 
            ENDIF 
C 
         END DO 
C 
      END IF 
C 
C     Exhaust Start emissions. 
C 
      IF (DBSELEFT(2) == 2 .AND. VVSTARTS(IV) == 1 .AND. 
     &    DBSELFAC(5) == 2) THEN 
C 
         EFTYPE = 2 
         IFAC = 5 
C 
C        PREPARE HEADER INFORMATION FOR EITHER DATABASE FORMAT 
C 
         WRITE(DBHEADER,110) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &     DBNSCEN,TABCHAR,IP,TABCHAR,IV,TABCHAR,EFTYPE,TABCHAR, 
     &     IFAC,TABCHAR,JDX-1,TABCHAR 
C 
         IF (DBDAILY) THEN 
C 
            CALL DB_DAILY(HEF(1:24,5),HSRFAC,DBHEADER,RFE, 
     &                    IFAC,IV,JDX,MY) 
C 
         ELSEIF (.NOT.DBAGGR) THEN 
C 
C           Write the remaining database fields in hourly format. (The 
C           remaining fields are the hour, gm/mi, gm/hr, Starts, Ends, 
C           Miles, MPG, HrVMT, FacVMT, RegDist, VCount, AmbTemp, and 
C           DiurTemp, and MY. Ends is null for all ETypes except Hot Soak. 
C 
            WRITE(DBUNIT(1),115) (DBHEADER(1:25),IH,TABCHAR, 
     &        HEF(IH,5),TABCHAR,HEF(IH,5)*BMYMPD(JDX,IV)*HVMT(IH), 
     &        TABCHAR,STARTSPERDAY(JDX,IV,IW)*STARTDIST(IH,IW),TABCHAR, 
     &        STARTSPERDAY(JDX,IV,IW)*HSRFAC*HS_HRLY_TRP_FRCN(IH,IW), 
     &        TABCHAR,BMYMPD(JDX,IV)*HVMT(IH),TABCHAR,RFE,TABCHAR, 
     &        HVMT(IH),TABCHAR,1.0,TABCHAR,MEVMYR(JDX,IV),TABCHAR, 
     &        DBVCOUNT(IV),TABCHAR,TEMHRLY(IH),TABCHAR,TEMHRLY(IH), 
     &        TABCHAR,MY,IH=DBSELHR(LOWER_LMT),DBSELHR(UPPER_LMT)) 
         END IF 
C 
      END IF 
C 
   99 RETURN 
      END SUBROUTINE DBEXHRS 
