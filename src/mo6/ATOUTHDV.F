      SUBROUTINE ATOUTHDV(HDTYPE) 
C 
C  OUTHDV writes out the expanded emission factor output table for HDGV (if 
C  HDTYPE=1) or HDDV (if HDTYPE=2) on unit IOUAT using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  10 Jul 01 AIR Task 15: New Routine 
C 
C  Called by ATOUTDSC 
C 
C  Calls None 
C 
C  Input on call: 
C 
C    argument list: HDTYPE 
C    modules: DATABASE 
C    common blocks: 
C    /IOUCOM/ IOUAT 
C    /VMXCOM/ VMTMIX 
C    /ATOX4 / EFAT 
C    /IDLCOM/ IDLEATFLG 
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
C  DASH8    C*8  8 character string of hyphens 
C  HDTYPE    I   heavy duty vehicle classes 
C  IAT       I   AT pollutant loop index 
C  IV        I   vehicle loop index 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  NAMPOL   C*5  AT pollutant names 
C  NAMVEH   C*8  vehicle class names 
C 
C 
      USE DATABASE, ONLY : DBSELTOX 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER  HDTYPE 
      INTEGER  I 
      INTEGER  IAT 
      INTEGER  IN 
      INTEGER  IV 
      INTEGER  IVFIRST 
      INTEGER  IVLAST 
      INTEGER  IXV 
      REAL     TOTAL(8) 
      REAL     XVTOT(8,2) 
C 
      CHARACTER(LEN=8) :: DASH8='  ------' 
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
      CHARACTER(LEN=8), DIMENSION(24) :: NAMVEH(24)=(/ 
     1   ' * N/A *', 
     2   ' * N/A *', 
     3   ' * N/A *', 
     4   ' * N/A *', 
     5   ' * N/A *', 
     6   '  HDGV2B', 
     7   '  HDGV3 ', 
     8   '  HDGV4 ', 
     9   '  HDGV5 ', 
     +   '  HDGV6 ', 
     1   '  HDGV7 ', 
     2   '  HDGV8A', 
     3   '  HDGV8B', 
     4   ' * N/A *', 
     5   ' * N/A *', 
     6   '  HDDV2B', 
     7   '  HDDV3 ', 
     8   '  HDDV4 ', 
     9   '  HDDV5 ', 
     +   '  HDDV6 ', 
     1   '  HDDV7 ', 
     2   '  HDDV8A', 
     3   '  HDDV8B', 
     4   ' * N/A *'/) 
C 
C 
C  Determine output range 
C 
      IF (HDTYPE.EQ.1) THEN 
        IVFIRST = 6 
        IVLAST  = 13 
      ELSEIF (HDTYPE.EQ.2) THEN 
        IVFIRST = 16 
        IVLAST  = 23 
      ELSE 
        GOTO 85 
      ENDIF 
C 
      WRITE(IOUAT,230) (NAMVEH(IV),IV=IVFIRST,IVLAST), 
     * (DASH8,IV=1,8),(VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',8(2X,A8),/T21,8(2X,A8), 
     *      /T13,'VMT Mix:',8(2X,F8.4)) 
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
        IF(DBSELTOX(IAT).EQ.2) 
     *  WRITE(IOUAT,250) NAMPOL(IAT), 
     *    (EFAT(IXV,IAT,IV),IV=IVFIRST,IVLAST) 
  250   FORMAT(T3,A17,':',8F10.2) 
C 
      IF(IAT.LE.2) THEN 
      DO IV=IVFIRST,IVLAST 
      I=IV-IVFIRST+1 
      XVTOT(I,IAT)=XVTOT(I,IAT)+EFAT(IXV,IAT,IV) 
      END DO 
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
     *      (EFAT(IXV,IAT,IV),IV=IVFIRST,IVLAST) 
C 
      DO IV=IVFIRST,IVLAST 
      I=IV-IVFIRST+1 
      TOTAL(I)=TOTAL(I)+EFAT(IXV,IAT,IV) 
      END DO 
C 
          END DO 
C 
      DO IV=IVFIRST,IVLAST 
      I=IV-IVFIRST+1 
      XVTOT(I,IAT)=XVTOT(I,IAT)+TOTAL(I) 
      END DO 
C 
      IF(IAT.EQ.1) THEN 
      WRITE(IOUAT,251) TOTAL 
  251   FORMAT(T3,'Benzene Total Evp:',8F10.2) 
      ELSE 
      WRITE(IOUAT,252) TOTAL 
  252   FORMAT(T3,'   MTBE Total Evp:',8F10.2) 
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
     *WRITE(IOUAT,255) (XVTOT(IV-IVFIRST+1,1),IV=IVFIRST,IVLAST) 
  255   FORMAT(T3,'Benzene Exh + Evp:',8F10.2) 
      IF(DBSELTOX(2).EQ.2)  
     *WRITE(IOUAT,256) (XVTOT(IV-IVFIRST+1,2),IV=IVFIRST,IVLAST) 
  256   FORMAT(T3,'   MTBE Exh + Evp:',8F10.2) 
C 
      WRITE(IOUAT,259) 
  259 FORMAT(1X,119('-')) 
C 
   85 RETURN 
      END 
