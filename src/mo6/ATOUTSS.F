      SUBROUTINE ATOUTSS(ICY) 
C 
C  ATOUTSS writes the AT emissions data to 
C  the SPREADSHEET database table. The output 
C  units are g/mi. 
C 
C  Called by:  ATCALX 
C 
C  Changes: (Last change first) 
C 
C  39 Aug 01 AIR Task 22: New routine. 
C 
C  Called by ATCALX 
C 
C  Calls ATOUTLDT, ATOUTHDV, ATOUTBUS 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    modules : None 
C    common blocks: 
C    see below 
C 
C  Local array subscripts: 
C 
C  NAMPOL(8)  -  NAMPOL ( IAT ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  ITER      I   write statement index 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  NAMPOL   C*4  pollutant names 
C  NAMUNT   C*12 component units, HCFLAG=3,5: mixed, HCFLAG=4: g/mi 
C  NAMVEH   C*8  vehicle class names 
C  SUM       R   summation variable 
C 
C 
      USE DATABASE, ONLY : DBNFILE,DBNRUN,DBNSCEN,TABCHAR,DBVEHNAMES, 
     *                     DBSCTITLE 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'  ! MAXIAT,ATSSFLAG,AROM,ATOXY,BENZ,E200,E300,OLEF      
      INCLUDE 'ATOX4.I'  ! EFAT 
      INCLUDE 'CITRV1.I' ! RVPBAS 
      INCLUDE 'EVAL.I'   ! MEVAL 
      INCLUDE 'FLAGS2.I' ! IMFLAG 
      INCLUDE 'IOUCOM.I' ! IOUATSS 
      INCLUDE 'IVTYPE.I' ! MAXVEH 
      INCLUDE 'NGV.I'    ! NGVFLG1 
      INCLUDE 'PART1.I'  ! D_PPM,MAXIPPM 
      INCLUDE 'REGION.I' ! IREJN 
      INCLUDE 'SPEED0.I' ! USRSPD 
      INCLUDE 'STRING.I' ! NOYES 
      INCLUDE 'SULFUR.I' ! S_PPM 
      INCLUDE 'TEMPS.I'  ! TEMMAX,TEMMIN 
      INCLUDE 'VMXCOM.I' ! VMTMIX 
      INCLUDE 'WEEKLY.I' ! IWEEK 
C 
C     Arguments: 
C 
      INTEGER, INTENT(IN) :: ICY 
C 
C     Local variables:                          
C 
      INTEGER I 
      INTEGER IAT 
      INTEGER IMF 
      INTEGER INGV 
      INTEGER IOXY 
      INTEGER IP 
      INTEGER IV 
      INTEGER IXV 
      INTEGER,DIMENSION(6)::ETYPE=(/1,3,4,6,5,8/) 
C 
      CHARACTER(20)  ::  SSHEADER 
      CHARACTER(50)  ::  SSTRAIL1 
      CHARACTER(90)  ::  SSTRAIL2 
C 
      CHARACTER(LEN=17),DIMENSION(16)::NAMPOL=(/ 
     *'Benzene          ', 
     *'MTBE             ', 
     *'1,3 Butadiene    ', 
     *'Formaldehyde     ', 
     *'Acetaldehyde     ', 
     *'Acrolein         ', 
     *'Benzene Hot Soak ', 
     *'Benzene Diurnal  ', 
     *'Benzene Running  ', 
     *'Benzene Resting  ', 
     *'Benzene Refueling', 
     *'MTBE Hot Soak    ', 
     *'MTBE Diurnal     ', 
     *'MTBE Running     ', 
     *'MTBE Resting     ', 
     *'MTBE Refueling   '/) 
C 
C 
C  If no spreadsheet output selected, return. 
C 
      IF(.NOT.ATSSFLAG) RETURN 
C 
C  Write the spreadsheet column headings, only once. 
C 
      IF(DBNFILE.EQ.1 .AND. 
     *   DBNRUN .EQ.1 .AND. 
     *   DBNSCEN.EQ.1)  
     *WRITE(IOUATSS,101)(TABCHAR,I=1,6), 
     *                 ((TRIM(DBVEHNAMES(IV)),TABCHAR),IV=1,MAXVEH), 
     *                  (TABCHAR,I=1,23) 
 101  FORMAT('File',A1,'Run',A1,'Scenario',A1, 
     *       'CY',A1,'EType',A1,'Pol Name',A1, 
     *       28(A,A1), 
     *       'Month',A1,'Altitude',A1,'TMin',A1,'Tmax',A1, 
     *       'Nom RVP',A1,'Gas Sulfur',A1,'Dsl Sulfur',A1, 
     *       'I/M?',A1,'Avg Spd',A1,'NGV?',A1,'E200',A1, 
     *       'E300',A1,'Arom',A1,'Olef',A1,'Benz',A1, 
     *       'MTBE vol%',A1,'MTBE MktFrac',A1, 
     *       'ETBE vol%',A1,'ETBE MktFrac',A1, 
     *       'ETOH vol%',A1,'ETOH MktFrac',A1, 
     *       'TAME vol%',A1,'TAME MktFrac',A1, 
     *       'Description') 
C 
C  Set the I/M index 
C 
      IF (IMFLAG.EQ.1) THEN 
        IMF = 1 
      ELSE 
        IMF = 2 
      ENDIF 
C 
C  Set the NGV index 
C 
      IF(NGVFLG1.EQ.0) THEN 
        INGV=1 
      ELSE 
        INGV=2 
      END IF 
C 
C  Create header section, Columns 1-4 
C 
      WRITE(SSHEADER,110) DBNFILE,TABCHAR, 
     &                    DBNRUN,TABCHAR, 
     &                    DBNSCEN,TABCHAR, 
     &                    ICY 
  110 FORMAT(3(I3,A1),I4) 
C 
C  Create the descriptive trailer, Columns 36-46 
C 
      WRITE(SSTRAIL1,111) TABCHAR,MEVAL, 
     &                    TABCHAR,IREJN, 
     &                    TABCHAR,TEMMIN, 
     &                    TABCHAR,TEMMAX, 
     &                    TABCHAR,RVPBAS, 
     &                    TABCHAR,S_PPM, 
     &                    TABCHAR,D_PPM, 
     &                    TABCHAR,NOYES(IMF), 
     &                    TABCHAR,USRSPD, 
     &                    TABCHAR,NOYES(INGV) 
  111 FORMAT(2(A1,I1),3(A1,F5.1),2(A1,F5.0),A1,A,A1,F5.1,A1,A) 
C 
C  Create the air toxic trailer, Columns 47-59 
C 
      WRITE(SSTRAIL2,112) TABCHAR,E200, 
     &                    TABCHAR,E300, 
     &                    TABCHAR,AROM, 
     &                    TABCHAR,OLEF, 
     &                    TABCHAR,BENZ, 
     &                  ((TABCHAR,ATOXY(I,IOXY),I=3,2,-1),IOXY=1,4) 
  112 FORMAT(5(A1,F5.2),8(A1,F5.3)) 
C 
C  Write the VMTMIX line 
C 
      WRITE(IOUATSS,230) TRIM(SSHEADER),TABCHAR,TABCHAR, 
     *                    (TABCHAR,VMTMIX(IV),IV=1,MAXVEH), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  230 FORMAT(A,A1,'VMT',A1,'0',28(A1,F10.4),2A,A1,A) 
C 
C     Print the exhaust emission factor group. 
C 
      IXV = 1 ! Exhaust index for AT 
C 
      DO IAT=1,MAXIAT 
C 
       IP=IAT+MAXIPPM+4  ! AT exhaust pollutant numbers start after classic and PM pollutants 
C 
      WRITE(IOUATSS,231) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                    (TABCHAR,EFAT(IXV,IAT,IV),IV=1,MAXVEH), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  231 FORMAT(A,A1,I1,A1,I2,28(A1,F10.2),2A,A1,A) 
C 
      END DO 
C 
C     Print the evap emission factor group. 
C 
      DO IAT=1,2 
        DO IXV=2,6 
C 
        IP = MAXIPPM+4+MAXIAT + (IAT-1)*5 + IXV-1  ! AT evap pollutant numbers start after AT exhaust 
C 
      WRITE(IOUATSS,231) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                    (TABCHAR,EFAT(IXV,IAT,IV),IV=1,MAXVEH), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
C 
      END DO 
      END DO 
C 
      RETURN 
      END 
