      SUBROUTINE CO2OUTDB(ICY,INERR) 
C 
C  CO2OUTDB writes the CO2 emissions data to 
C  the EMISSIONS database table. The output 
C  units are g/mi. 
C 
C  Called by:  ATCALX 
C 
C  Calls:  CO2DB_DAILY, CO2DB_AGGR 
C 
C  Changes: (Last change first) 
C 
C 14 Nov 01: AIR Task 27: New MOBILE6 routine. 
C 
C    Argument list: ICY 
C 
C  Notes: 
C 
C 
      USE DATABASE, ONLY : DBFLAG,DBSELPOL,DBSELVEH,DBSELEFT,DBSELFAC, 
     &                     DBNFILE,DBNRUN,DBNSCEN,DBDAILY,DBAGGR, 
     &                     DBUNIT,TABCHAR,DBSELAGE,DBSELHR 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'CEVBMY.I'     ! BMYMPD 
      INCLUDE  'CO2COM.I'     ! DBCO2 
      INCLUDE  'EVAHS1.I'     ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'FLAGS5.I'     ! FTP_INPUT_FLAG 
      INCLUDE  'IVTYPE.I'     ! VVSTARTS,VVDSL,VTLDGV,VVHEAVY 
      INCLUDE  'MAXIMA.I'     ! MAXYRS 
      INCLUDE  'MYRCAL.I'     ! MEVMYR,TF 
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
      CHARACTER(26)                ::  DBHEADER    ! Non-changing part of database output record. 
C 
      INTEGER                      ::  ETYPE       ! Emission type index 
      INTEGER                      ::  IDX         ! Yearly loop 
      INTEGER                      ::  IFAC        ! Facility type index. 
      INTEGER                      ::  IH          ! Hourly interval. 
      INTEGER                      ::  IP          ! Pollutant number 
      INTEGER                      ::  IV          ! Vehicle class 
      INTEGER                      ::  IW          ! Weekday/weekend flag. 
      INTEGER                      ::  JDX         ! Age 
      INTEGER                      ::  MY          ! Model year 
      INTEGER, SAVE                ::  LOWER_LMT=1 ! Lower Index into DBSELAGE and DBSELHR 
      INTEGER, SAVE                ::  UPPER_LMT=2 ! Upper index into DBSELAGE and DBSELHR 
C 
      REAL                         ::  HSRFAC      ! Hot soak fraction 
      REAL                         ::  RFE         ! Fuel economy value. 
      REAL,DIMENSION(24,5)         ::  HEF         ! Hourly emission factor 
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
C     Write the CO2 output only if the database selection criteria include 
C     the pollutant, vehicle type, and age. 
C 
      IP = 4 ! CO2 
C 
       IF(DBSELPOL(IP).EQ.1) RETURN 
C 
       DO IV=1,MAXVEH 
C 
        IF(DBSELVEH(IV).EQ.1) CYCLE 
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
      IDX=(MAXYRS+1)-JDX 
C 
      MY=ICY-JDX+1 
C 
      RFE = PMFLECON(MY,IV,INERR) 
C 
C     CO2 only has Running Exhaust Etype 
C 
      ETYPE=1 
C 
      IF(DBSELEFT(ETYPE).EQ.1) RETURN 
C 
C     Load the 24-hour array HEF. 
C 
      DO IFAC=1,4 
       DO IH=1,24 
        HEF(IH,IFAC)=DBCO2(IFAC,IH,IDX,IV)  ! DBCO2 units are g/mi 
       END DO 
      END DO 
C 
         DO IFAC = 1,4 
C 
            IF (DBSELFAC(IFAC) .NE. 2) CYCLE 
C 
C           WRITE THE HEADER 
C 
            WRITE(DBHEADER,110) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &        DBNSCEN,TABCHAR,IP,TABCHAR,IV,TABCHAR,ETYPE,TABCHAR, 
     &        IFAC,TABCHAR,JDX-1,TABCHAR 
  110 FORMAT(I3,A1,I3,A1,I3,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1) 
C 
            IF (DBDAILY) THEN 
C 
                CALL CO2DB_DAILY(HEF(1:24,IFAC),HSRFAC,DBHEADER, 
     &                        RFE,IFAC,IV,JDX,MY) 
C 
            ELSEIF (.NOT.DBAGGR) THEN 
C 
C  Loop over desired hours 
C 
           DO IH=DBSELHR(LOWER_LMT),DBSELHR(UPPER_LMT) 
C 
C  Write the remaining database fields for hourly output. (The 
C  remaining fields are the hour, gm/mi, gm/hr, Starts, Ends, 
C  Miles, MPG, HrVMT, FacVMT, RegDist, VCount, AmbTemp, and 
C  DiurTemp, and MY. 
C 
        WRITE(DBUNIT(1),115) DBHEADER(1:26),IH, 
     &  TABCHAR, HEF(IH,IFAC), 
     &  TABCHAR, HEF(IH,IFAC)*BMYMPD(JDX,IV)*HVMT(IH), 
     &  TABCHAR, STARTSPERDAY(JDX,IV,IW)*STARTDIST(IH,IW), 
     &  TABCHAR, STARTSPERDAY(JDX,IV,IW)*HSRFAC*HS_HRLY_TRP_FRCN(IH,IW), 
     &  TABCHAR, BMYMPD(JDX,IV)*HVMT(IH), 
     &  TABCHAR, RFE, 
     &  TABCHAR, HVMT(IH), 
     &  TABCHAR, FVMT(IFAC,IH,IV), 
     &  TABCHAR, MEVMYR(JDX,IV), 
     &  TABCHAR, DBVCOUNT(IV), 
     &  TABCHAR, TEMHRLY(IH), 
     &  TABCHAR, TEMHRLY(IH), 
     &  TABCHAR, MY 
C 
  115 FORMAT(A26,I2,A1,F9.1,A1,F9.1,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F10.4,A1,F6.1,A1,F6.1, 
     &       A1,I4) 
C 
              END DO !IH 
C 
            ENDIF 
C 
         END DO !IFAC 
C 
      END DO !ETYPE 
      END DO !JDX 
C 
      IF(DBAGGR) CALL CO2DB_AGGR(ICY,INERR) 
C 
   99 RETURN 
C 
      END SUBROUTINE CO2OUTDB 
