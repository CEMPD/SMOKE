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
C  08 Mar 02 AIR Task 15,18: Changed formats to ES to accomodate the 
C            wide range of ADDITIONAL HAPS values. Added USE ATHEAP 
C            due to ADD HAPS dynamic array expansion. Updated printing 
C            logic. 
C  04 Dec 01 AIR Task 27-28: Include ALL VEH. 
C  28 Sep 01 AIR Task 25: Changed ATSSFLAG to SSFLAG to accomodate PM output. 
C            Added Particle Size to output. 
C  10 Sep 01 Air Task 18: Added HAPS 
C  30 Aug 01 AIR Task 22: New routine. 
C 
C  Called by ATCALX 
C 
C  Calls NONE 
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
C 
C 
      USE DATABASE, ONLY : DBNFILE,DBNRUN,DBNSCEN,TABCHAR,DBVEHNAMES, 
     *                     DBSCTITLE 
      USE ATHEAP, ONLY: EFAT, VAT, ATHSNAME, ATHIV, ATHPMAP 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'  ! MAXIAT,SSFLAG,AROM,ATOXY,BENZ,E200,E300,OLEF      
      INCLUDE 'ATOX5.I'  ! ATHAPS,NMAP 
      INCLUDE 'CITRV1.I' ! RVPBAS 
      INCLUDE 'EVAL.I'   ! MEVAL 
      INCLUDE 'FLAGS2.I' ! IMFLAG 
      INCLUDE 'IOUCOM.I' ! IOUSS 
      INCLUDE 'IVTYPE.I' ! MAXVEH,VVDSL 
      INCLUDE 'NGV.I'    ! NGVFLG1 
      INCLUDE 'PART1.I'  ! D_PPM,MAXIPPM,PSZCUT 
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
      INTEGER IMAP 
      INTEGER IMF 
      INTEGER INGV 
      INTEGER IOXY 
      INTEGER IP 
      INTEGER IRB 
      INTEGER IV 
      INTEGER IVC 
      INTEGER IXV 
      INTEGER JXV 
      INTEGER,DIMENSION(6)::ETYPE=(/1,3,4,6,5,8/) ! 1=Exh 3=Hot Soak 4=Diurnal 5=Resting 6=Running 8=Refueling 
C 
      CHARACTER(10)  ::  CEF(29) 
      CHARACTER(30)  ::  SSHEADER 
      CHARACTER(50)  ::  SSTRAIL1 
      CHARACTER(95)  ::  SSTRAIL2 
C 
      LOGICAL        ::  PRNT 
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
      IF(.NOT.SSFLAG) RETURN 
C 
C  Write the spreadsheet column headings, only once. Include ALL VEH 
C 
      IF(DBNFILE.EQ.1 .AND. 
     *   DBNRUN .EQ.1 .AND. 
     *   DBNSCEN.EQ.1)  
     *WRITE(IOUSS,101)(TABCHAR,I=1,7), 
     *                 (TRIM(DBVEHNAMES(IV)),TABCHAR,IV=1,MAXVEH+1), 
     *                  (TABCHAR,I=1,24) 
 101  FORMAT('File',A1,'Run',A1,'Scenario',A1, 
     *       'CY',A1,'EType',A1,'Pol #',A1,'Pol Name',A1, 
     *       29(A,A1), 
     *       'Month',A1,'Altitude',A1,'TMin',A1,'Tmax',A1, 
     *       'Nom RVP',A1,'Gas Sulfur',A1,'Dsl Sulfur',A1, 
     *       'I/M?',A1,'Avg Spd',A1,'NGV?',A1,'E200',A1, 
     *       'E300',A1,'Arom',A1,'Olef',A1,'Benz',A1, 
     *       'MTBE vol%',A1,'MTBE MktFrac',A1, 
     *       'ETBE vol%',A1,'ETBE MktFrac',A1, 
     *       'ETOH vol%',A1,'ETOH MktFrac',A1, 
     *       'TAME vol%',A1,'TAME MktFrac',A1, 
     *       'Part Size',A1,'Description') 
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
     &                  ((TABCHAR,ATOXY(I,IOXY),I=3,2,-1),IOXY=1,4), 
     &                    TABCHAR,PSZCUT 
  112 FORMAT(5(A1,F5.2),8(A1,F5.3),A1,F5.2) 
C 
C  Write the VMTMIX line, include ALL VEH 
C 
      WRITE(IOUSS,230) TRIM(SSHEADER),TABCHAR,TABCHAR,TABCHAR, 
     *                    (TABCHAR,VMTMIX(IV),IV=1,MAXVEH),TABCHAR, 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  230 FORMAT(A,A1,'VMT',A1,'0',A1,'VMT',28(A1,F10.4),A1,'1',2A,A1,A) 
C 
C     Print the exhaust emission factor group. 
C 
      IXV = 1 ! Exhaust index for AT 
C 
      DO IAT=1,MAXIAT 
C 
       IP=IAT+MAXIPPM+4  ! AT exhaust pollutant numbers start after classic and PM pollutants 
C 
      WRITE(IOUSS,231) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                 TABCHAR,TRIM(NAMPOL(IAT)), 
     *                 (TABCHAR,EFAT(IXV,IAT,IV),IV=1,MAXVEH), 
     *                  TABCHAR,VAT(IXV,IAT), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
  231 FORMAT(A,A1,I1,A1,I2,A1,A,29(A1,F10.2),2A,A1,A) 
C 
      END DO 
C 
C     Print the HAPS exhaust emission factor group. 
C 
      IF(ATHAPS) THEN 
C 
      JXV = 1 ! Exhaust index for HAPS 
C 
      DO IMAP=1,NMAP 
C 
      IP=ATHPMAP(IMAP)  ! HAPS exhaust pollutant numbers are between 50 and 99 
      IAT=MAXIAT+IMAP 
C 
      PRNT = .FALSE. 
      IVC=0 
      CEF='NA' 
c 
      DO IV=1,MAXVEH  ! Convert data to characters so that "NA" can signify missing data 
      DO IRB=1,4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
      WRITE(CEF(IV),10) EFAT(IXV,IAT,IV) 
   10 FORMAT(ES10.3) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      EXIT 
      END IF 
      END DO 
      END DO 
      IF(IVC.EQ.MAXVEH) WRITE(CEF(29),10) VAT(IXV,IAT) ! ALL VEH 
C 
      IF(PRNT) 
     *WRITE(IOUSS,232) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                   TABCHAR,TRIM(ATHSNAME(IMAP)), 
     *                    (TABCHAR,TRIM(CEF(IV)),IV=1,MAXVEH+1), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  232 FORMAT(A,A1,I1,A1,I2,A1,A,29(A1,A),2A,A1,A) 
C 
      END DO 
      END IF !HAPS 
C 
C     Print the evap emission factor group. 
C 
      DO IAT=1,2 
C 
       IP=IAT+MAXIPPM+4  ! AT exhaust pollutant numbers start after classic and PM pollutants 
C 
        DO IXV=2,6 
C 
      WRITE(IOUSS,231) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                   TABCHAR,TRIM(NAMPOL(MAXIAT+(IAT-1)*5+IXV-1)), 
     *                    (TABCHAR,EFAT(IXV,IAT,IV),IV=1,MAXVEH), 
     *                     TABCHAR,VAT(IXV,IAT), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
C 
      END DO 
      END DO 
C 
C     Print the HAPS evap emission factor group. 
C 
      IF(ATHAPS) THEN 
C 
      JXV = 2 ! Evap index for HAPS 
C 
      DO IMAP=1,NMAP 
C 
      IP=ATHPMAP(IMAP)  ! HAPS exhaust pollutant numbers are between 50 and 99 
      IAT=MAXIAT+IMAP 
C 
      DO IXV=2,6 
C 
      PRNT = .FALSE. 
      IVC=0 
      CEF='NA' 
C 
      DO IV=1,MAXVEH  ! Convert data to characters so that "NA" can signify missing data 
C 
      IF(VVDSL(IV).EQ.1) THEN   ! Diesel EVAP ADD HAPS are always zero. They alone do not trigger printing. 
        WRITE(CEF(IV),10) EFAT(IXV,IAT,IV) 
        IVC=IVC+1 
        CYCLE 
      END IF 
      DO IRB=1,4 
C 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
      WRITE(CEF(IV),10) EFAT(IXV,IAT,IV) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      EXIT 
      END IF 
      END DO 
      END DO 
      IF(IVC.EQ.MAXVEH) WRITE(CEF(29),10) VAT(IXV,IAT) ! ALL VEH 
C 
      IF(PRNT) 
     *WRITE(IOUSS,232) TRIM(SSHEADER),TABCHAR,ETYPE(IXV),TABCHAR,IP, 
     *                   TABCHAR,TRIM(ATHSNAME(IMAP)), 
     *                    (TABCHAR,TRIM(CEF(IV)),IV=1,MAXVEH+1), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
C 
      END DO !IXV 
      END DO !IMAP 
      END IF !HAPS 
C 
      RETURN 
      END 
