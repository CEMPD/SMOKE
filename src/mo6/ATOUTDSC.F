      SUBROUTINE ATOUTDSC(ICY) 
C 
C  ATOUTDSC writes descriptive output for scenario-level Air Toxic and  
C  HAPS data. 
C 
C  Changes: (Last change first) 
C 
C  15 Jan 02 AIR Task 15: Removed Reformulated Gas Yes/No from the output 
C                         since RFG cannot be used with Air Toxic runs. 
C  12 Sep 01 AIR Task 11-22: Fixed HDGV output bug. 
C  11 Sep 01 AIR Task 18: Added HAPS output. 
C  10 Jul 01 AIR Task 15: New routine. 
C 
C  Called by ATCALX 
C 
C  Calls ATOUTLDT, ATOUTHDV, ATOUTBUS 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    modules : None 
C 
C  Local array subscripts: 
C 
C  NAMPOL(8)  -  NAMPOL ( IAT ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C  CETYPE(6)  -  CETYPE ( IXV ) 
C 
C  Local variable / array dictionary (see also below): 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CETYPE   C*4  E-type name 
C  DASH8    C*8  8 character string of hyphens 
C  NAMPOL   C*4  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C 
C 
      USE DATABASE, ONLY: DBSELTOX 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'   ! MAXIAT,ATOXY,E200,E300,AROM,OLEF,BENZ 
      INCLUDE 'ATOX2.I'   ! ATWTOXY 
      INCLUDE 'ATOX4.I'   ! EFAT,VAT 
      INCLUDE 'ATOX5.I'   ! NMAP,ATHIV 
      INCLUDE 'CITRV1.I'  ! RVPICY 
      INCLUDE 'EVAL.I'    ! MONTH, MEVAL 
      INCLUDE 'IOUCOM.I'  ! IOUAT 
      INCLUDE 'IVTYPE.I'  ! MAXVEH,ILDGT12,ILDGT34,IHDGV,ILDDV,ILDDT,IHDDV,ILDGT (IVTYPE.I) 
      INCLUDE 'OPCNTL.I'  ! NODESC 
      INCLUDE 'STRING.I'  ! NOYES,COMMA,PERIOD 
      INCLUDE 'SULFUR.I'  ! S_PPM 
      INCLUDE 'TEMPS.I'   ! TEMMAX,TEMMIN 
      INCLUDE 'VMXCOM.I'  ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTHDD 
      INCLUDE 'WEEKLY.I'  ! IWEEK 
C 
C     Arguments: 
C 
      INTEGER, INTENT(IN) :: ICY 
C 
C     External functions: 
C 
      REAL ATEFCVT  
C 
C     Local variables:                          
C 
      CHARACTER*10 CEF(10) 
      INTEGER I               ! Temporary DO loop variable 
      INTEGER IAT             ! Air Toxic index 
      INTEGER IC              ! Count for a small vehicle group 
      INTEGER IMAP            ! HASP index, up to NMAP 
      INTEGER IN              ! Index to print the pollutant name (NAMPOL) 
      INTEGER IOXY            ! Oxygenate index 
      INTEGER IRB             ! Index for HAPS ratios and BERs 
      INTEGER ITER            ! Index to print out dashses 
      INTEGER IV              ! Vehicle class 
      INTEGER IVC             ! Count for all vehicles 
      INTEGER IXV             ! Exhaust (1), Evap (2-6) AT index 
      INTEGER JXV             ! Exhuast (1), Evap (2) HAPS index 
      LOGICAL PRNT            ! Printing switch 
      REAL    SUM             ! Summation variable 
      REAL    TOTAL(10)       ! Total Evap 
      REAL    XVTOT(10,2)     ! Exhaust + Evap total 
C 
      CHARACTER(LEN=9), DIMENSION(10) :: NAMVEH=(/ 
     *   '    LDGV', 
     *   '  LDGT12', 
     *   '  LDGT34', 
     *   '    LDGT', 
     *   '    HDGV', 
     *   '    LDDV', 
     *   '    LDDT', 
     *   '    HDDV', 
     *   '      MC', 
     *   ' All Veh'/) 
C 
      CHARACTER(LEN=17),DIMENSION(16)::NAMPOL=(/ 
     *'          Benzene', 
     *'             MTBE', 
     *'    1,3 Butadiene', 
     *'     Formaldehyde', 
     *'     Acetaldehyde', 
     *'         Acrolein', 
     *'Benzene Hot Soak ', 
     *'Benzene Diurnal  ', 
     *'Benzene Running  ', 
     *'Benzene Resting  ', 
     *'Benzene Refueling', 
     *'   MTBE Hot Soak ', 
     *'   MTBE Diurnal  ', 
     *'   MTBE Running  ', 
     *'   MTBE Resting  ', 
     *'   MTBE Refueling'/) 
C 
      CHARACTER(LEN=9) :: DASH8='  ------' 
C 
      CHARACTER(LEN=4),DIMENSION(6)::CETYPE=(/ 
     *'Exht', 
     *'HtSk', 
     *'Drnl', 
     *'Rung', 
     *'Rest', 
     *'Refu'/) 
C 
C 
C  If no descriptive output selected, return. 
C 
      IF(NODESC) RETURN 
C 
      WRITE(IOUAT, 400) 
  400 FORMAT() 
C 
C  Write header. Note that oxygenats are output in their original 
C  vol% format. However, the weighted oxygen level is in wt%. 
C 
   10 WRITE(IOUAT,210) ICY,MONTH(MEVAL), 
     *        ATWTOXY,S_PPM,TEMMAX,TEMMIN,RVPICY,E200,E300, 
     *        AROM,OLEF,BENZ,((ATOXY(I,IOXY),I=3,2,-1),IOXY=1,4) 
  210 FORMAT(T15,'                Calendar Year: ',I5/ 
     *       T15,'                        Month: ',1X,A4/ 
     *       T15,' Market Weighted Oxygen Level: ',F5.3,' wt%'/ 
     *       T15,' Gasoline Fuel Sulfur Content: ',F5.0,' ppm'/ 
C     *       T15,'   Diesel Fuel Sulfur Content: ',F5.0,' ppm'/ 
     *       T15,'          Maximum Temperature: ',F5.1,' F'/ 
     *       T15,'          Minimum Temperature: ',F5.1,' F'/ 
     *       T15,'                Weathered RVP: ',F5.1,' psi'/ 
     *       T15,'                         E200: ',F5.2,' %'/ 
     *       T15,'                         E300: ',F5.2,' %'/ 
     *       T15,'                    Aromatics: ',F5.2,' vol%'/ 
     *       T15,'                      Olefins: ',F5.2,' vol%'/ 
     *       T15,'                      Benzene: ',F5.2,' vol%'/ 
     *       T15,'                         MTBE: ',F5.3,' vol% ', 
     *                        '(market fraction: ',F5.3,')'/ 
     *       T15,'                         ETBE: ',F5.3,' vol% ', 
     *                        '(market fraction: ',F5.3,')'/ 
     *       T15,'                      Ethanol: ',F5.3,' vol% ', 
     *                        '(market fraction: ',F5.3,')'/ 
     *       T15,'                         TAME: ',F5.3,' vol% ', 
     *                        '(market fraction: ',F5.3,')'/) 
C 
      IF( IWEEK .EQ. 2 ) THEN 
        WRITE (IOUAT,217) 
      END IF 
  217 FORMAT(T1,'Emissions determined from WEEKEND hourly ', 
     *'vehicle activity fractions.'/) 
C 
C  Write VMTMIX summary 
C 
      SUM=VMTMIX(VTLDGV)+VMTGT12+VMTGT34+VMTHDG 
     *   +VMTMIX(VTLDDV)+VMTLDDT+VMTHDD+VMTMIX(VTMC) 
C 
      WRITE(IOUAT,230) 
     *   NAMVEH,(DASH8,ITER=1,10), 
     *   VMTMIX(VTLDGV),VMTGT12,VMTGT34,VMTHDG, 
     *   VMTMIX(VTLDDV),VMTLDDT,VMTHDD,VMTMIX(VTMC),SUM 
  230 FORMAT(T8,'Vehicle Type:',10(2X,A8)/ 
     *       T16,'GVWR:',T36,'<6000',T46,'>6000',T56,'(All)'/ 
     *       T21,10(2X,A8)/ 
     *       T4,'VMT Distribution:',3F10.4,10X,6F10.4) 
C 
C     Initialize Exh + Evap Total 
C 
      XVTOT = 0. 
C 
C********************************************************************** 
C 
C     Print the exhaust emission factor group. 
C 
      WRITE(IOUAT,240) 
  240 FORMAT(1X,119('-'),/T2, 
     *   'Exhaust Emission Factors (mg/mi):') 
C 
      IXV = 1 ! Exhaust 
C 
      DO IAT=1,MAXIAT 
C 
        IF(DBSELTOX(IAT).EQ.2) THEN 
          WRITE(IOUAT,250) NAMPOL(IAT), 
     *    EFAT(IXV,IAT,VTLDGV), 
     *    ATEFCVT(ILDGT12,IXV,IAT), 
     *    ATEFCVT(ILDGT34,IXV,IAT), 
     *    ATEFCVT(ILDGT,IXV,IAT), 
     *    ATEFCVT(IHDGV,IXV,IAT), 
     *    EFAT(IXV,IAT,VTLDDV), 
     *    ATEFCVT(ILDDT,IXV,IAT), 
     *    ATEFCVT(IHDDV,IXV,IAT), 
     *    EFAT(IXV,IAT,VTMC), 
     *    VAT(IXV,IAT) 
  250   FORMAT(T3,A17,':',9F10.2,F10.3) 
C 
      IF(IAT.LE.2) THEN 
       XVTOT( 1,IAT)=XVTOT( 1,IAT)+EFAT(IXV,IAT,VTLDGV) 
       XVTOT( 2,IAT)=XVTOT( 2,IAT)+ATEFCVT(ILDGT12,IXV,IAT) 
       XVTOT( 3,IAT)=XVTOT( 3,IAT)+ATEFCVT(ILDGT34,IXV,IAT) 
       XVTOT( 4,IAT)=XVTOT( 4,IAT)+ATEFCVT(ILDGT,IXV,IAT) 
       XVTOT( 5,IAT)=XVTOT( 5,IAT)+ATEFCVT(IHDGV,IXV,IAT) 
       XVTOT( 6,IAT)=XVTOT( 6,IAT)+EFAT(IXV,IAT,VTLDDV) 
       XVTOT( 7,IAT)=XVTOT( 7,IAT)+ATEFCVT(ILDDT,IXV,IAT) 
       XVTOT( 8,IAT)=XVTOT( 8,IAT)+ATEFCVT(IHDDV,IXV,IAT) 
       XVTOT( 9,IAT)=XVTOT( 9,IAT)+EFAT(IXV,IAT,VTMC) 
       XVTOT(10,IAT)=XVTOT(10,IAT)+VAT(IXV,IAT) 
      END IF 
C 
      ENDIF 
C 
      END DO 
C 
C********************************************************************** 
C 
C     Print the evap emission factor group. 
C 
      WRITE(IOUAT,260) 
  260 FORMAT(1X,119('-'),/T2, 
     *   'Evaporative Emission Factors (mg/mi):') 
C 
      DO IAT=1,2 
C 
        TOTAL = 0. 
C 
        IF(DBSELTOX(IAT).EQ.2) THEN 
          DO IXV=2,6 
            IN = MAXIAT + (IAT-1)*5 + IXV-1 
            WRITE(IOUAT,250) NAMPOL(IN), 
     *      EFAT(IXV,IAT,VTLDGV), 
     *      ATEFCVT(ILDGT12,IXV,IAT), 
     *      ATEFCVT(ILDGT34,IXV,IAT), 
     *      ATEFCVT(ILDGT,IXV,IAT), 
     *      ATEFCVT(IHDGV,IXV,IAT), 
     *      EFAT(IXV,IAT,VTLDDV), 
     *      ATEFCVT(ILDDT,IXV,IAT), 
     *      ATEFCVT(IHDDV,IXV,IAT), 
     *      EFAT(IXV,IAT,VTMC), 
     *      VAT(IXV,IAT) 
C 
      TOTAL( 1)=TOTAL( 1)+EFAT(IXV,IAT,VTLDGV) 
      TOTAL( 2)=TOTAL( 2)+ATEFCVT(ILDGT12,IXV,IAT) 
      TOTAL( 3)=TOTAL( 3)+ATEFCVT(ILDGT34,IXV,IAT) 
      TOTAL( 4)=TOTAL( 4)+ATEFCVT(ILDGT,IXV,IAT) 
      TOTAL( 5)=TOTAL( 5)+ATEFCVT(IHDGV,IXV,IAT) 
      TOTAL( 6)=TOTAL( 6)+EFAT(IXV,IAT,VTLDDV) 
      TOTAL( 7)=TOTAL( 7)+ATEFCVT(ILDDT,IXV,IAT) 
      TOTAL( 8)=TOTAL( 8)+ATEFCVT(IHDDV,IXV,IAT) 
      TOTAL( 9)=TOTAL( 9)+EFAT(IXV,IAT,VTMC) 
      TOTAL(10)=TOTAL(10)+VAT(IXV,IAT) 
C             
      END DO 
C 
      DO I=1,10 
      XVTOT(I,IAT)=XVTOT(I,IAT)+TOTAL(I) 
      END DO 
C 
      IF(IAT.EQ.1) THEN 
      WRITE(IOUAT,251) TOTAL 
  251   FORMAT(T3,'Benzene Total Evp:',9F10.2,F10.3) 
      ELSE 
      WRITE(IOUAT,252) TOTAL 
  252   FORMAT(T3,'   MTBE Total Evp:',9F10.2,F10.3) 
      END IF 
C       
      END IF 
C 
      END DO 
C 
C********************************************************************** 
C 
C     Print the exhaust+evap emission factor group. 
C 
      WRITE(IOUAT,254) 
  254 FORMAT(1X,119('-'),/T2, 
     *   'Exhaust + Evaporative Emission Factors (mg/mi):') 
C 
      IF(DBSELTOX(1).EQ.2)  
     *WRITE(IOUAT,255) (XVTOT(I,1),I=1,10) 
  255   FORMAT(T3,'Benzene Exh + Evp:',9F10.2,F10.3) 
      IF(DBSELTOX(2).EQ.2)  
     *WRITE(IOUAT,256) (XVTOT(I,2),I=1,10) 
  256   FORMAT(T3,'   MTBE Exh + Evp:',9F10.2,F10.3) 
C 
      WRITE(IOUAT,259) 
  259 FORMAT(1X,119('-')) 
C 
C********************************************************************** 
C 
C     Print the HAPS exhaust emission factor group. 
C 
      IF(ATHAPS) THEN 
C 
      WRITE(IOUAT,241) 
  241 FORMAT(T2, 
     *   'HAPS Exhaust Emission Factors (mg/mi):') 
C 
      JXV = 1 ! Exhaust index for HAPS 
      IXV = 1 ! Exhaust index 
C 
      DO IMAP=1,NMAP 
C 
      IAT=MAXIAT+IMAP 
C 
      PRNT = .FALSE. 
C 
      CEF='        NA' 
      IVC=0 
C 
      DO IRB=1,4 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDGV)) THEN 
       WRITE(CEF(1),20) EFAT(IXV,IAT,ILDGV) 
   20  FORMAT(F10.2) 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
C 
      IC=0 
      DO IV=VTLDGT1,VTLDGT2 
       IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.2) WRITE(CEF(2),20) ATEFCVT(ILDGT12,IXV,IAT) 
C 
      IC=0 
      DO IV=VTLDGT3,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.2) WRITE(CEF(3),20) ATEFCVT(ILDGT34,IXV,IAT) 
C 
      IC=0 
      DO IV=VTLDGT1,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) IC=IC+1 
      END DO 
      IF(IC.EQ.4) WRITE(CEF(4),20) ATEFCVT(ILDGT,IXV,IAT) 
C 
      IC=0 
      DO IV=VTHDGV2B,VTHDGV8B 
       IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.8) WRITE(CEF(5),20) ATEFCVT(IHDGV,IXV,IAT) 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDV)) THEN 
      WRITE(CEF(6),20) EFAT(IXV,IAT,VTLDDV) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      END IF 
C 
      IC=0 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT12)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT34)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      IF(IC.EQ.2) WRITE(CEF(7),20) ATEFCVT(ILDDT,IXV,IAT) 
C 
      IC=0 
      DO IV=VTHDDV2B,VTHDDV8B 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.8) WRITE(CEF(8),20) ATEFCVT(IHDDV,IXV,IAT) 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTMC)) THEN 
      WRITE(CEF(9),20) EFAT(IXV,IAT,VTMC) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      END IF 
C 
      IF(IVC.EQ.28) WRITE(CEF(10),25) VAT(IXV,IAT) 
   25 FORMAT(F10.3) 
C 
      END DO !IRB 
C 
      IF(PRNT) WRITE(IOUAT,232) TRIM(ATHSNAME(IMAP)),CETYPE(IXV),CEF 
  232 FORMAT(T3,A,' ',A,T20,':',10A10) 
C 
      END DO !IMAP 
C 
C********************************************************************** 
C 
C     Print the HAPS evap emission factor group. 
C 
      WRITE(IOUAT,242) 
  242 FORMAT(1X,119('-'),/T2, 
     *   'HAPS Evaporative Emission Factors (mg/mi):') 
C 
      JXV = 2 ! Evap index for HAPS 
C 
      DO IMAP=1,NMAP 
C 
      IAT=MAXIAT+IMAP 
C 
      DO IXV=2,6 
C 
      PRNT = .FALSE. 
C 
      CEF='        NA' 
      IVC=0 
C 
      DO IRB=1,4 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDGV)) THEN 
       WRITE(CEF(1),20) EFAT(IXV,IAT,ILDGV) 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
C 
      IC=0 
      DO IV=VTLDGT1,VTLDGT2 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.2) WRITE(CEF(2),20) ATEFCVT(ILDGT12,IXV,IAT) 
C 
      IC=0 
      DO IV=VTLDGT3,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.2) WRITE(CEF(3),20) ATEFCVT(ILDGT34,IXV,IAT) 
C 
      IC=0 
      DO IV=VTLDGT1,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) IC=IC+1 
      END DO 
      IF(IC.EQ.4) WRITE(CEF(4),20) ATEFCVT(ILDGT,IXV,IAT) 
C 
      IC=0 
      DO IV=VTHDGV2B,VTHDGV8B 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.8) WRITE(CEF(5),20) ATEFCVT(IHDGV,IXV,IAT) 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDV)) THEN 
      WRITE(CEF(6),20) EFAT(IXV,IAT,VTLDDV) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      END IF 
C 
      IC=0 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT12)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT34)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      IF(IC.EQ.2) WRITE(CEF(7),20) ATEFCVT(ILDDT,IXV,IAT) 
C 
      IC=0 
      DO IV=VTHDDV2B,VTHDDV8B 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
       IC=IC+1 
       IVC=IVC+1 
       PRNT = .TRUE. 
      END IF 
      END DO 
      IF(IC.EQ.8) WRITE(CEF(8),20) ATEFCVT(IHDDV,IXV,IAT) 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTMC)) THEN 
      WRITE(CEF(9),20) EFAT(IXV,IAT,VTMC) 
      IVC=IVC+1 
      PRNT = .TRUE. 
      END IF 
C 
      IF(IVC.EQ.28) WRITE(CEF(10),25) VAT(IXV,IAT) 
C 
      END DO !IRB 
C 
      IF(PRNT) WRITE(IOUAT,232) TRIM(ATHSNAME(IMAP)),CETYPE(IXV),CEF 
C 
      END DO !IXV 
C 
      END DO !IMAP 
C 
      WRITE(IOUAT,259) 
C 
      END IF !HAPS 
C 
C********************************************************************** 
C 
C  Check the expanded output requests. 
C 
   85 IF (OPTABL(1).EQ.2) CALL ATOUTLDT 
      IF (OPTABL(2).EQ.2) CALL ATOUTHDV(1) 
      IF (OPTABL(3).EQ.2) CALL ATOUTHDV(2) 
      IF (OPTABL(4).EQ.2) CALL ATOUTBUS 
C 
      RETURN 
      END 
