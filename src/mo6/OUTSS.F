      SUBROUTINE OUTSS(ICY) 
C 
C  OUTSS writes the HC, CO, NOx, CO2 emissions data to 
C  the SPREADSHEET database table. The output units are g/mi. 
C 
C  Changes: (Last change first) 
C 
C  11 Mar 02: AIR Task Bug384: Change logic so that diesel evaps are 
C             reported as zero, instead of NA (per ADD HAPS change.) 
C  04 Dec 01: AIR Task 29: Added ALL VEH, MPG, and EVAP output. 
C  27 Nov 01: AIR Task 28: New routine. 
C 
C  Called by:  PRCSCN 
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
C  NAMPOL(8)  -  NAMPOL ( IP ) 
C  NAMNMH(10) -  NAMNMH ( NMHFLG ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C 
      USE DATABASE, ONLY : DBNFILE,DBNRUN,DBNSCEN,TABCHAR,DBVEHNAMES, 
     *                     DBSCTITLE,DBSELPOL 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'  ! SSFLAG,AROM,ATOXY,BENZ,E200,E300,OLEF   
      INCLUDE 'CITRV1.I' ! RVPBAS 
      INCLUDE 'CO2COM.I' ! MPG, FMPG 
      INCLUDE 'EVAL.I'   ! MEVAL 
      INCLUDE 'FLAGS2.I' ! IMFLAG 
      INCLUDE 'FLAGS4.I' ! NMHFLG 
      INCLUDE 'IOUCOM.I' ! IOUSS 
      INCLUDE 'IVTYPE.I' ! MAXVEH,VVGAS,VVDSL 
      INCLUDE 'MAXIMA.I' ! MAXPOL 
      INCLUDE 'NGV.I'    ! NGVFLG1 
      INCLUDE 'PART1.I'  ! D_PPM,PSZCUT 
      INCLUDE 'REGION.I' ! IREJN 
      INCLUDE 'RESUL1.I' ! EFFTP    
      INCLUDE 'RESUL2.I' ! G_MI_EVP    
      INCLUDE 'RESUL3.I' ! VFTP    
      INCLUDE 'RESUL4.I' ! TOTEVAP, VTEVAP, VDRNL, VHTSK 
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
      INTEGER IMF 
      INTEGER INGV 
      INTEGER IOXY 
      INTEGER IP 
      INTEGER IV 
C 
      INTEGER, SAVE :: IDIURNAL=2,IHOTSOAK=1 
      INTEGER, SAVE :: ICRNCASE=4,IRESTLOSS=5 
C 
      CHARACTER( 5)  ::  EMISS 
      CHARACTER(10)  ::  EVAP 
      CHARACTER(10)  ::  CEF(29) 
      CHARACTER(30)  ::  SSHEADER 
      CHARACTER(50)  ::  SSTRAIL1 
      CHARACTER(95)  ::  SSTRAIL2 
C 
      CHARACTER(LEN=5), DIMENSION(5) :: NAMNMH=(/ 
     *                   ' THC ',' NMHC',' VOC ',' TOG ',' NMOG'/) 
C 
      CHARACTER(LEN=5),DIMENSION(4)::NAMPOL=(/ 
     *                   ' HC  ',' CO  ',' NOX ',' CO2 '/) 
C 
C  If no spreadsheet output selected, return. 
C 
      IF(.NOT.SSFLAG) RETURN 
C 
C  Write the spreadsheet column headings, only once. Include ALL VEH. 
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
C  Write the MPG line if CO2 calculated, include ALL VEH 
C 
      IF(DBSELPOL(4).EQ.2) 
     *WRITE(IOUSS,220) TRIM(SSHEADER),TABCHAR,TABCHAR,TABCHAR, 
     *                    (TABCHAR,MPG(IV),IV=1,MAXVEH),TABCHAR,FMPG, 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  220 FORMAT(A,A1,'MPG',A1,'0',A1,'MPG',29(A1,F10.1),2A,A1,A) 
C 
C  Write the VMTMIX line, include ALL VEH 
C 
      WRITE(IOUSS,230) TRIM(SSHEADER),TABCHAR,TABCHAR,TABCHAR, 
     *                    (TABCHAR,VMTMIX(IV),IV=1,MAXVEH),TABCHAR, 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  230 FORMAT(A,A1,'VMT',A1,'0',A1,'VMT',28(A1,F10.4),A1,'1',2A,A1,A) 
C 
C  HC, CO, NOX and CO2 emissions are Running Exhaust ETYPEs (1). 
C 
      ETYPE = 1 
C 
C  Print the exhaust emission factors, including CO2 
C 
      DO IP=1,MAXPOL+1 
       IF(DBSELPOL(IP).EQ.1) CYCLE 
C 
        IF (IP.EQ.1) THEN 
          EMISS = NAMNMH(NMHFLG) 
        ELSE 
          EMISS = NAMPOL(IP) 
        ENDIF 
C 
C  Convert data to characters to allow for formatting changes 
C 
      DO IV=1,MAXVEH   
        IF(IP.LE.3) THEN 
        WRITE(CEF(IV),10) EFFTP(IP,IV) 
   10   FORMAT(F10.3) 
        ELSE 
        WRITE(CEF(IV),20) EFFTP(IP,IV) 
   20   FORMAT(F10.1) 
        END IF 
      END DO 
C 
        IF(IP.LE.3) THEN 
        WRITE(CEF(29),10) VFTP(IP) ! ALL VEH 
        ELSE 
        WRITE(CEF(29),20) VFTP(IP) ! ALL VEH 
        END IF 
C 
      WRITE(IOUSS,232) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                   TABCHAR,TRIM(EMISS), 
     *                    (TABCHAR,TRIM(CEF(IV)),IV=1,MAXVEH+1), 
     *                     TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                     TABCHAR,TRIM(DBSCTITLE) 
  232 FORMAT(A,A1,I2,A1,I2,A1,A,29(A1,A),2A,A1,A) 
C 
      END DO 
C 
C  Print out the evap constituents, if available.  
C 
      IP=1 
      IF(DBSELPOL(IP).EQ.2) THEN 
C 
C  Hot Soak 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) G_MI_EVP(IHOTSOAK,IV) 
      END DO 
      WRITE(CEF(29),10) VHTSK 
C 
      EVAP='Hot Soak' 
      ETYPE=3 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
  290 FORMAT(A,A1,I2,A1,I2,A1,A,29(A1,A),2A,A1,A) 
C 
C  Diurnal 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) G_MI_EVP(IDIURNAL,IV) 
      END DO 
      WRITE(CEF(29),10) VDRNL 
C 
      EVAP='Diurnal' 
      ETYPE=4 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
C  Resting 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) G_MI_EVP(IRESTLOSS,IV) 
      END DO 
      WRITE(CEF(29),10) VRSTLS 
C 
      EVAP='Resting' 
      ETYPE=5 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
C  Running 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) EFRUNL(IV) 
      END DO 
      WRITE(CEF(29),10) VRUNLS 
C 
      EVAP='Running' 
      ETYPE=6 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
C  Crankcase 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) G_MI_EVP(ICRNCASE,IV) 
      END DO 
      WRITE(CEF(29),10) VCRNK 
C 
      EVAP='Crankcase' 
      ETYPE=7 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
C  Refueling 
C 
      CEF='NA' 
      DO IV=1,MAXVEH 
      IF(VVGASMC(IV).EQ.1) WRITE(CEF(IV),10) REFL_GM_MI(IV) 
      END DO 
C 
      EVAP='Refueling' 
      ETYPE=8 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
C  Total Evap 
C 
      DO IV=1,MAXVEH 
      WRITE(CEF(IV),10) TOTEVAP(IV) 
      END DO 
      WRITE(CEF(29),10) VTEVAP 
C 
      EVAP='Total Evap' 
      ETYPE=9 
C 
      WRITE(IOUSS,290) TRIM(SSHEADER),TABCHAR,ETYPE,TABCHAR,IP, 
     *                 TABCHAR,TRIM(EVAP), 
     *                (TABCHAR,CEF(IV),IV=1,MAXVEH+1), 
     *                 TRIM(SSTRAIL1),TRIM(SSTRAIL2), 
     *                 TABCHAR,TRIM(DBSCTITLE) 
C 
      END IF 
C 
      RETURN 
      END 
