      SUBROUTINE PMOUTDB(ICY,INERR) 
C 
C  PMOUTDB writes the PM emissions data to 
C  the EMISSIONS database table. 
C 
C  Called by:  PMCALX 
C 
C  CallS:  DB_DAILY 
C 
C  Changes: (Last change first) 
C 
C 15 Nov 01: Added INERR processing. 
C 19 Apr 01: Added code to prevent output of PM IDLE. 
C 17 Apr 01: Added code to stop output of GASPM with diesels 
C            and ECARBON/OCARBON with gas. 
C 06 Mar 01: New MOBILE6 routine. 
C 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: ICY 
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
      USE DATABASE, ONLY : DBFLAG,DBSELPART,DBSELVEH,DBSELEFT,DBSELFAC, 
     &                     DBNFILE,DBNRUN,DBNSCEN,DBDAILY,DBAGGR, 
     &                     DBUNIT,TABCHAR,DBSELAGE,DBSELHR 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'CEVBMY.I'     ! BMYMPD 
      INCLUDE  'EVAHS1.I'     ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'FLAGS5.I'     ! FTP_INPUT_FLAG 
      INCLUDE  'IVTYPE.I'     ! VVSTARTS 
      INCLUDE  'MAXIMA.I'     ! MAXYRS 
      INCLUDE  'MYRCAL.I'     ! MEVMYR,TF 
      INCLUDE  'PART1.I'      ! MAXIPPM 
      INCLUDE  'PART3.I'      ! DBPM 
      INCLUDE  'SPEED9.I'     ! FVMT,HVMT 
      INCLUDE  'STDIST.I'     ! STARTDIST, STARTSPERDAY 
      INCLUDE  'TEMPS.I'      ! TEMHRLY 
      INCLUDE  'VMXCOM.I'     ! DBVCOUNT 
      INCLUDE  'WEEKLY.I'     ! IWEEK 
C 
      INTEGER,INTENT(IN)               ::  ICY 
      INTEGER,INTENT(INOUT)            ::  INERR 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
      CHARACTER(26)                ::  DBHEADER 
C 
      INTEGER                      ::  EFTYPE 
      INTEGER                      ::  IFAC 
      INTEGER                      ::  IH 
      INTEGER                      ::  IW 
      INTEGER, SAVE                ::  UPPER_LMT=2 
      INTEGER, SAVE                ::  LOWER_LMT=1 
      INTEGER                      ::  IPPM 
      INTEGER                      ::  IP 
      INTEGER                      ::  MY 
      INTEGER                      ::  JDX 
      INTEGER                      ::  IV 
C 
      REAL                         ::  HSRFAC 
      REAL                         ::  RFE 
      REAL,DIMENSION(24,5)         ::  HEF 
C 
C 
C     Don't write the database output if not selected or 
C     if either FTP output option was selected. 
C 
      IF(DBFLAG.EQ.1) RETURN 
C 
      IF (FTP_INPUT_FLG /= 0) RETURN 
C 
C     Initialize 
C 
      IW = IWEEK 
C 
C     Write the output only if the database selection criteria include 
C     the PM pollutant, vehicle type, and age. 
C 
      DO IPPM=1,MAXIPPM 
C 
       IF(DBSELPART(IPPM).EQ.1) CYCLE 
C 
       IF(IPPM.EQ.11) CYCLE ! Prevent DB output for PM Idle 
C 
       IP=IPPM+6 
C 
       DO IV=1,MAXVEH 
C 
        IF(DBSELVEH(IV).EQ.1) CYCLE 
C 
C     Skip if GASPM for DIESELS, OCARBON/ECARBON for GAS 
C 
        IF(IPPM.EQ.4.AND.VVDSL(IV).EQ.1 .OR. 
     *     IPPM.GE.2.AND.IPPM.LE.3.AND.VVGASMC(IV).EQ.1) CYCLE 
C 
C     Set up 
C 
        IF (IV == VTLDGV) THEN 
          HSRFAC = HSRED(1) 
        ELSE 
          HSRFAC = HSRED(2) 
        ENDIF 
C 
        DO JDX=1,MAXYRS 
C 
        IF(DBSELAGE(LOWER_LMT) .GT. (JDX-1) .OR. 
     &     DBSELAGE(UPPER_LMT) .LT. (JDX-1)) CYCLE 
C 
      MY=ICY-JDX+1 
C 
      RFE = PMFLECON(MY,IV,INERR) 
C 
C     Load the 24-hour array HEF, for use with PMDB_DAILY and PMDB_AGGR below 
C 
      DO IFAC=1,4 
       DO IH=1,24 
        HEF(IH,IFAC)=DBPM(IPPM,IFAC,IH,JDX,IV) 
       END DO 
      END DO 
C 
C     All PM emissions are Running Exhaust ETYPEs (1), except for Brake and Tire, 
C     which have their own ETYPE categories (9 and 10). 
C 
      IF (DBSELEFT(1) .EQ. 2 .AND. IPPM .LE. 7) THEN 
C 
         EFTYPE = 1 
C 
         DO IFAC = 1,4 
C 
            IF (DBSELFAC(IFAC) .NE. 2) CYCLE 
C 
C           WRITE THE HEADER 
C 
            WRITE(DBHEADER,110) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &        DBNSCEN,TABCHAR,IP,TABCHAR,IV,TABCHAR,EFTYPE,TABCHAR, 
     &        IFAC,TABCHAR,JDX-1,TABCHAR 
  110 FORMAT(I3,A1,I3,A1,I3,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1) 
C 
            IF (DBDAILY) THEN 
C 
                CALL PMDB_DAILY(HEF(1:24,IFAC),HSRFAC,DBHEADER, 
     &                        RFE,IFAC,IV,JDX,MY) 
C 
            ELSEIF (.NOT.DBAGGR) THEN 
C 
C               Write the remaining database fields for hourly output. (The 
C               remaining fields are the hour, gm/mi, gm/hr, Starts, Ends, 
C               Miles, MPG, HrVMT, FacVMT, RegDist, VCount, AmbTemp, and 
C               DiurTemp, and MY. 
C 
                WRITE(DBUNIT(1),115) (DBHEADER(1:26),IH,TABCHAR, 
     &           HEF(IH,IFAC),TABCHAR, 
     &           HEF(IH,IFAC)*BMYMPD(JDX,IV)*HVMT(IH),TABCHAR, 
     &           STARTSPERDAY(JDX,IV,IW)*STARTDIST(IH,IW),TABCHAR, 
     &           STARTSPERDAY(JDX,IV,IW)*HSRFAC*HS_HRLY_TRP_FRCN(IH,IW), 
     &           TABCHAR,BMYMPD(JDX,IV)*HVMT(IH),TABCHAR,RFE,TABCHAR, 
     &           HVMT(IH),TABCHAR,FVMT(IFAC,IH,IV),TABCHAR, 
     &           MEVMYR(JDX,IV),TABCHAR,DBVCOUNT(IV),TABCHAR, 
     &           TEMHRLY(IH),TABCHAR,TEMHRLY(IH),TABCHAR,MY, 
     &           IH=DBSELHR(LOWER_LMT),DBSELHR(UPPER_LMT)) 
C 
  115 FORMAT(A26,I2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F10.4,A1,F6.1,A1,F6.1, 
     &       A1,I4) 
C 
            ENDIF 
C 
         END DO !IFAC 
C 
      END IF 
C 
C     Brake or Tire Emissions. Since they are independent of 
C     road (facility) type, use the arterial value in calls to 
C     PMDB_DAILY and general DB output (eg HEF(1:24,1) and HEF(IH,1)) 
C 
      IF (DBSELEFT( 9) .EQ. 2 .AND. IPPM .EQ. 8 .OR. 
     &    DBSELEFT(10) .EQ. 2. AND. IPPM .EQ. 9 ) THEN 
C 
         EFTYPE = IPPM + 1 
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
            CALL PMDB_DAILY(HEF(1:24,1),HSRFAC,DBHEADER,RFE, 
     &                    IFAC,IV,JDX,MY) 
C 
         ELSEIF (.NOT.DBAGGR) THEN 
C 
C           Write the remaining database fields in hourly format. (The 
C           remaining fields are the hour, gm/mi, gm/hr, Starts, Ends, 
C           Miles, MPG, HrVMT, FacVMT, RegDist, VCount, AmbTemp, and 
C           DiurTemp, and MY. Ends is null for all ETypes except Hot Soak. 
C 
            WRITE(DBUNIT(1),115) (DBHEADER(1:26),IH,TABCHAR, 
     &        HEF(IH,1),TABCHAR,HEF(IH,1)*BMYMPD(JDX,IV)*HVMT(IH), 
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
      END DO !JDX 
      END DO !IV 
      END DO !IPPM 
C 
      IF(DBAGGR) CALL PMDB_AGGR(ICY,INERR) 
C 
   99 RETURN 
C 
      END SUBROUTINE PMOUTDB 
