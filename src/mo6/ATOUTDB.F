      SUBROUTINE ATOUTDB(ICY,INERR) 
C 
C  ATOUTDB writes the AT emissions data to 
C  the EMISSIONS database table. The output 
C  units are g/mi. 
C 
C  Called by:  ATCALX 
C 
C  Calls:  ATDB_DAILY, ATDB_AGGR 
C 
C  Changes: (Last change first) 
C 
C 28 Nov 01: Air TASK 28: Added ATHEAP for AT dynamic array support. 
C 15 Nov 01 AIR Task 27: Added INERR support due to PMFLECON update. 
C 05 OCT 01: AIr TASK 15-22: Added filter to prevent output of evap from 
C            start/exhaust only species. 
C 20 Sep 01: AIR Task 15-22: Comments 
C 10 Jul 01: AIR Task 15-22: New MOBILE6 routine. 
C 
C    Argument list: ICY 
C 
C  Notes: 
C 
C 
      USE DATABASE, ONLY : DBFLAG,DBSELTOX,DBSELVEH,DBSELEFT,DBSELFAC, 
     &                     DBNFILE,DBNRUN,DBNSCEN,DBDAILY,DBAGGR, 
     &                     DBUNIT,TABCHAR,DBSELAGE,DBSELHR 
      USE ATHEAP 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'ATOX1.I'      ! MAXIAT 
      INCLUDE  'ATOX4.I'      ! DBAT 
      INCLUDE  'CEVBMY.I'     ! BMYMPD 
      INCLUDE  'EVAHS1.I'     ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'FLAGS5.I'     ! FTP_INPUT_FLAG 
      INCLUDE  'IVTYPE.I'     ! VVSTARTS,VVDSL,VTLDGV,VVHEAVY 
      INCLUDE  'MAXIMA.I'     ! MAXYRS 
      INCLUDE  'MYRCAL.I'     ! MEVMYR,TF 
      INCLUDE  'PART1.I'      ! MAXIPPM 
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
      INTEGER                      ::  IAT         ! Air Toxic pollutant index 
      INTEGER                      ::  IB          ! Beginning road type index 
      INTEGER                      ::  IDX         ! Yearly loop 
      INTEGER                      ::  IE          ! Ending road type index 
      INTEGER                      ::  IFAC        ! Facility type index. 
      INTEGER                      ::  IH          ! Hourly interval. 
      INTEGER                      ::  IP          ! Pollutant number 
      INTEGER                      ::  IV          ! Vehicle class 
      INTEGER                      ::  IW          ! Weekday/weekend flag. 
      INTEGER                      ::  IXV         ! eXhaust/eVap index 
      INTEGER                      ::  JDX         ! Age 
      INTEGER                      ::  MY          ! Model year 
      INTEGER, SAVE                ::  LOWER_LMT=1 ! Lower Index into DBSELAGE and DBSELHR 
      INTEGER, SAVE                ::  UPPER_LMT=2 ! Upper index into DBSELAGE and DBSELHR 
C 
      REAL                         ::  HSRFAC      ! Hot soak fraction 
      REAL                         ::  RFE         ! Fuel economy value. 
      REAL,DIMENSION(24,5)         ::  HEF         ! Hourly emission factor 
      REAL                         ::  FACVMT      ! Road type VMT fraction 
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
C     the AT pollutant, vehicle type, and age. 
C 
      DO IAT=1,MAXIAT 
C 
       IF(DBSELTOX(IAT).EQ.1) CYCLE 
C 
       IP=IAT+MAXIPPM+4  ! AT emission numbers start after classic and PM pollutants 
C 
       DO IV=1,MAXVEH 
C 
        IF(DBSELVEH(IV).EQ.1) CYCLE 
C 
C     Exhaust MTBE is output for gasoline vehicles only 
C 
        IF(VVDSL(IV).EQ.1 .AND. IAT.EQ.2) CYCLE 
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
C     Loop over the ETYPES 
C 
          DO ETYPE=1,8 
C 
            IF(DBSELEFT(ETYPE).EQ.1) CYCLE 
C 
C     Heavy duty vehicles and busses do not have engine start emissions. 
C 
            IF(VVHEAVY(IV).EQ.1 .AND. ETYPE.EQ.2) CYCLE 
C 
C     Prevent output of evap for all diesels 
C 
            IF(VVDSL(IV).EQ.1 .AND. ETYPE.GE.3) EXIT 
C 
C     Prevent output of evap emissions from anything but Benzene and MTBE. 
C 
            IF(IAT.GE.3 .AND. ETYPE.GE.3) EXIT 
C 
C     Map AT emissions into ETYPE space 
C 
            IF (ETYPE.LE.2) IXV=1  !Exhaust Running and Starts 
            IF (ETYPE.EQ.3) IXV=2  !Evap Hot Soak 
            IF (ETYPE.EQ.4) IXV=3  !Evap Diurnal 
            IF (ETYPE.EQ.5) IXV=5  !Evap Resting 
            IF (ETYPE.EQ.6) IXV=4  !Evap Running 
            IF (ETYPE.EQ.7) CYCLE  !Evap Crankcase (Not Available) 
            IF (ETYPE.EQ.8) IXV=6  !Evap Refueling 
C 
C     Set the facility (road) types to output. Exhaust Running and 
C     Evap Running Loss are road dependent (Roads 1-4). Starts and 
C     the remaining Evaps are road independent (Road 5) 
C 
            IF(ETYPE.EQ.1 .OR. ETYPE.EQ.6) THEN 
              IB=1 
              IE=4 
            ELSE 
              IB=5 
              IE=5 
            END IF 
C 
C     Load the 24-hour array HEF, for use with ATDB_DAILY below. 
C 
      DO IFAC=IB,IE 
       DO IH=1,24 
        HEF(IH,IFAC)=DBAT(IDX,IH,IFAC,IXV,IAT,IV)  ! DBAT units are g/mi 
       END DO 
      END DO 
C 
         DO IFAC = IB,IE 
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
                CALL ATDB_DAILY(HEF(1:24,IFAC),HSRFAC,DBHEADER, 
     &                        RFE,IFAC,IV,JDX,MY) 
C 
            ELSEIF (.NOT.DBAGGR) THEN 
C 
C  Loop over desired hours 
C 
           DO IH=DBSELHR(LOWER_LMT),DBSELHR(UPPER_LMT) 
C 
C  Select the road fraction 
C 
        IF(IFAC.LE.4) THEN 
          FACVMT=FVMT(IFAC,IH,IV) 
        ELSE 
          FACVMT=1.0 
        END IF 
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
     &  TABCHAR, FACVMT, 
     &  TABCHAR, MEVMYR(JDX,IV), 
     &  TABCHAR, DBVCOUNT(IV), 
     &  TABCHAR, TEMHRLY(IH), 
     &  TABCHAR, TEMHRLY(IH), 
     &  TABCHAR, MY 
C 
  115 FORMAT(A26,I2,A1,F9.6,A1,F9.6,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
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
      END DO !IV 
      END DO !IAT 
C 
      IF(DBAGGR) CALL ATDB_AGGR(ICY,INERR) 
C 
   99 RETURN 
C 
      END SUBROUTINE ATOUTDB 
