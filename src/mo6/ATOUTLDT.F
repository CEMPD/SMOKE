      SUBROUTINE ATOUTLDT() 
C 
C  OUTLDT writes out the expanded emission factor output table for 
C  LDTs, LDGT1, LDGT2, LDGT3, LDGT4, LDDT12 and LDDT34 - on unit 
C  IOUAT using descriptive output. 
C 
C  Called by ATOUTDSC. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  12 Mar 02 AIR Task Bug384 Output diesel evap 
C  08 Mar 02 AIR Task 15,18: Changed formats to ES to accomodate the 
C            wide range of ADDITIONAL HAPS values. Added USE ATHEAP 
C            due to ADD HAPS dynamic array expansion. 
C  11 Sep 01 AIR Task 18: Add HAPS output. 
C  10 Jul 01 AIR Task 15: New Routine 
C 
C  Input on call: 
C 
C    argument list: None. 
C    modules: DATABASE 
C    common blocks: 
C    /ATOX4/  EFAT 
C    /IOUCOM/ IOUAT 
C    /MAXIMA/ MAXPOL 
C    /VMXCOM/ VMTMIX 
C 
C  Local array subscripts: 
C 
C  NAMPOL(3)  -  NAMPOL ( IAT ) 
C  NAMVEH(6)  -  NAMVEH ( IVF3 ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NAMPOL   C*4  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  DASH8    C*8  8 character string of hyphens 
C 
C 
      USE DATABASE, ONLY : DBSELTOX 
      USE ATHEAP, ONLY : EFAT, VAT, ATHSNAME, ATHIV 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'VMXCOM.I' 
C 
      CHARACTER*10 CEF(6) 
      INTEGER I 
      INTEGER IAT 
      INTEGER IMAP 
      INTEGER IN 
      INTEGER IRB 
      INTEGER IV 
      INTEGER IXV 
      INTEGER JXV 
      LOGICAL PRNT 
      REAL    TOTAL(6) 
      REAL    XVTOT(6,2) 
C 
      CHARACTER(LEN=8), DIMENSION(6) :: NAMVEH=(/ 
     *   '  LDGT1 ', 
     *   '  LDGT2 ', 
     *   '  LDGT3 ', 
     *   '  LDGT4 ', 
     *   ' LDDT12 ', 
     *   ' LDDT34 '/) 
C 
      CHARACTER(LEN=9)::DASH8='  ------' 
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
      CHARACTER(LEN=4),DIMENSION(6)::CETYPE=(/ 
     *'Exht', 
     *'HtSk', 
     *'Drnl', 
     *'Rung', 
     *'Rest', 
     *'Refu'/) 
C 
C 
      WRITE(IOUAT,230) NAMVEH,(DASH8,I=1,6), 
     *                  (VMTMIX(IV),IV=VTLDGT1,VTLDGT4), 
     *                  VMTMIX(VTLDDT12),VMTMIX(VTLDDT34) 
  230 FORMAT(T11,'Veh. Type: ',6(2X,A8),/T21,6(2X,A8), 
     *      /T13,'VMT Mix:',6(2X,F8.4)) 
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
  240 FORMAT(1X,119('-'),/T2,'Exhaust Emission Factors (mg/mi):') 
C 
      IXV = 1  ! Exhaust 
C 
      DO IAT=1,MAXIAT 
C 
        IF(DBSELTOX(IAT).EQ.2) THEN 
        WRITE(IOUAT,250) NAMPOL(IAT), 
     *    (EFAT(IXV,IAT,IV),IV=VTLDGT1,VTLDGT4), 
     *     EFAT(IXV,IAT,VTLDDT12), 
     *     EFAT(IXV,IAT,VTLDDT34) 
  250   FORMAT(T3,A17,':',6F10.2) 
C 
      IF(IAT.LE.2) THEN 
      DO IV=VTLDGT1,VTLDGT4 
      I=IV-VTLDGT1+1 
      XVTOT(I,IAT)=XVTOT(I,IAT)+EFAT(IXV,IAT,IV) 
      END DO 
      XVTOT(5,IAT)=XVTOT(5,IAT)+EFAT(IXV,IAT,VTLDDT12) 
      XVTOT(6,IAT)=XVTOT(6,IAT)+EFAT(IXV,IAT,VTLDDT34) 
      END IF 
C 
      END IF 
C 
      END DO 
C 
C********************************************************************** 
C 
C     Print the evap emission factor group. 
C 
      WRITE(IOUAT,260) 
  260 FORMAT(1X,119('-'),/T2,'Evaporative Emission Factors (mg/mi):') 
C 
      DO IAT=1,2 
C 
        TOTAL = 0. 
C 
        IF(DBSELTOX(IAT).EQ.2) THEN 
          DO IXV=2,6 
            IN = MAXIAT + (IAT-1)*5 + IXV-1 
            WRITE(IOUAT,250) NAMPOL(IN), 
     *      (EFAT(IXV,IAT,IV),IV=VTLDGT1,VTLDGT4), 
     *       EFAT(IXV,IAT,VTLDDT12), 
     *       EFAT(IXV,IAT,VTLDDT34) 
C 
      DO IV=VTLDGT1,VTLDGT4 
      I=IV-VTLDGT1+1 
      TOTAL(I)=TOTAL(I)+EFAT(IXV,IAT,IV) 
      END DO 
      TOTAL(5)=TOTAL(5)+EFAT(IXV,IAT,VTLDDT12) 
      TOTAL(6)=TOTAL(6)+EFAT(IXV,IAT,VTLDDT34) 
C 
          END DO 
C 
      DO I=1,6 
      XVTOT(I,IAT)=XVTOT(I,IAT)+TOTAL(I) 
      END DO 
C 
      IF(IAT.EQ.1) THEN 
      WRITE(IOUAT,251) TOTAL 
  251   FORMAT(T3,'Benzene Total Evp:',6F10.2) 
      ELSE 
      WRITE(IOUAT,252) TOTAL 
  252   FORMAT(T3,'   MTBE Total Evp:',6F10.2) 
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
     *WRITE(IOUAT,255) (XVTOT(I,1),I=1,6) 
  255   FORMAT(T3,'Benzene Exh + Evp:',6F10.2) 
      IF(DBSELTOX(2).EQ.2)  
     *WRITE(IOUAT,256) (XVTOT(I,2),I=1,6) 
  256   FORMAT(T3,'   MTBE Exh + Evp:',6F10.2) 
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
C 
      DO IRB=1,4 
C 
      DO IV=VTLDGT1,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
      PRNT = .TRUE. 
      WRITE(CEF(IV-VTLDGT1+1),20) EFAT(IXV,IAT,IV) 
   20 FORMAT(ES10.3) 
      END IF 
      END DO 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT12)) THEN 
      PRNT = .TRUE. 
      WRITE(CEF(5),20) EFAT(IXV,IAT,VTLDDT12) 
      END IF 
C 
      IF(ATHIV(JXV,IRB,IMAP,VTLDDT34)) THEN 
      PRNT = .TRUE. 
      WRITE(CEF(6),20) EFAT(IXV,IAT,VTLDDT34) 
      END IF 
C 
      END DO !IRB 
C 
      IF(PRNT) WRITE(IOUAT,232) TRIM(ATHSNAME(IMAP)),CETYPE(IXV),CEF 
  232 FORMAT(T3,A,' ',A,T20,':',6A10) 
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
C 
      DO IRB=1,4 
C 
      DO IV=VTLDGT1,VTLDGT4 
      IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
      PRNT = .TRUE. 
      WRITE(CEF(IV-VTLDGT1+1),20) EFAT(IXV,IAT,IV) 
      END IF 
      END DO 
C 
      END DO !IRB 
C 
      WRITE(CEF(5),20) EFAT(IXV,IAT,VTLDDT12) 
C 
      WRITE(CEF(6),20) EFAT(IXV,IAT,VTLDDT34) 
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
   85 RETURN 
      END 
 
