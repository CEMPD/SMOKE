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
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER I 
      INTEGER IAT 
      INTEGER IN 
      INTEGER IV 
      INTEGER IXV 
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
   85 RETURN 
      END 
 
