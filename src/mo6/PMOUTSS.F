      SUBROUTINE PMOUTSS(ICY) 
C 
C  PMOUTSS writes the PM emissions data to 
C  the SPREADSHEET database table. The output 
C  units are g/mi. 
C 
C  Called by:  PMCALX 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 01 AIR Task 27-28: Include ALL VEH. 
C  28 Sep 01: AIR Task 25: Changed ATSSFLAG to SSFLAG to accomodate PM output. 
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
     *                     DBSCTITLE,DBSELPART 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'  ! SSFLAG,AROM,ATOXY,BENZ,E200,E300,OLEF      
      INCLUDE 'CITRV1.I' ! RVPBAS 
      INCLUDE 'EVAL.I'   ! MEVAL 
      INCLUDE 'FLAGS2.I' ! IMFLAG 
      INCLUDE 'IDLE.I'   ! IDLEPMFLG 
      INCLUDE 'IOUCOM.I' ! IOUSS 
      INCLUDE 'IVTYPE.I' ! MAXVEH,VVGAS,VVDSL 
      INCLUDE 'NGV.I'    ! NGVFLG1 
      INCLUDE 'PART1.I'  ! D_PPM,MAXIPPM,PSZCUT 
      INCLUDE 'PART3.I'  ! EFPM, VPM 
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
      INTEGER ETYPE 
      INTEGER I 
      INTEGER IMAP 
      INTEGER IMF 
      INTEGER INGV 
      INTEGER IOXY 
      INTEGER IP 
      INTEGER IPPM 
      INTEGER IV 
      INTEGER, SAVE, DIMENSION(11) :: MAPOUT=(/5,4,3,2,1,8,9,6,7,10,11/) 
C 
      CHARACTER(10)  ::  CEF(29) 
      CHARACTER(30)  ::  SSHEADER 
      CHARACTER(50)  ::  SSTRAIL1 
      CHARACTER(95)  ::  SSTRAIL2 
C 
      CHARACTER(LEN=8),DIMENSION(11)::NAMPOL= 
     *(/'     SO4', 
     *  ' OCARBON', 
     *  ' ECARBON', 
     *  '   GASPM', 
     *  '    Lead', 
     *  '     SO2', 
     *  '     NH3', 
     *  '   Brake', 
     *  '    Tire', 
     *  'EMFAC_PM', 
     *  ' PM Idle'/) 
C 
C  If no spreadsheet output selected, return. 
C 
      IF(.NOT.SSFLAG) RETURN 
C 
C  Write the spreadsheet column headings, only once. 
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
C     Print the PM exhaust emission factors 
C 
      DO IMAP=1,MAXIPPM 
C 
        IPPM=MAPOUT(IMAP) 
C 
C     Skip output if: EMFAC_PM, not selected, or idle not computed 
C 
      IF(IPPM.EQ.10) CYCLE 
C 
       IF(DBSELPART(IPPM).EQ.1) CYCLE 
C 
       IF(IPPM.EQ.11.AND..NOT.IDLEPMFLG) CYCLE 
C 
       IP=IPPM+6 
C 
C     All PM emissions are Running Exhaust ETYPEs (1), except for Brake and Tire, 
C     which have their own ETYPE categories (9 and 10). 
C 
      IF (IPPM .LE. 7 .OR.  
     &    IPPM .EQ. 11) THEN 
        ETYPE = 1 
      ELSE IF(IPPM.EQ.8) THEN 
        ETYPE = 9 
      ELSE IF(IPPM.EQ.9) THEN 
        ETYPE = 10 
      END IF 
C 
C  Convert data to characters so that "NA" can signify not available data 
C 
      DO IV=1,MAXVEH   
      CEF(IV)='NA' 
      IF((IPPM.EQ.4.OR.IPPM.EQ.5).AND.VVDSL(IV).EQ.1) CYCLE  ! No lead data or GASPM for diesels 
      IF((IPPM.EQ.2.OR.IPPM.EQ.3).AND.VVGASMC(IV).EQ.1) CYCLE  ! No ECARBON or OCARBON data for gas vehicles 
      WRITE(CEF(IV),10) EFPM(IPPM,IV) 
   10 FORMAT(F10.4) 
      END DO 
      WRITE(CEF(29),10) VPM(IPPM) ! All Veh 
C 
      WRITE(IOUSS,232) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                   TABCHAR,TRIM(NAMPOL(IPPM)), 
     *                    (TABCHAR,TRIM(CEF(IV)),IV=1,MAXVEH+1), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  232 FORMAT(A,A1,I2,A1,I2,A1,A,29(A1,A),2A,A1,A) 
C 
      END DO 
C 
      RETURN 
      END 
