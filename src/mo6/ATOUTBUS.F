      SUBROUTINE ATOUTBUS() 
C 
C  OUTBUS writes out the expanded emission factor output table for the BUS 
C  vehicles on unit IOUAT using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  12 Mar 02 AIR Task Bug384: Print diesel evap output 
C  08 Mar 02 AIR Task 15,18: Changed formats to ES to accomodate the 
C            wide range of ADDITIONAL HAPS values. Added USE ATHEAP 
C            due to ADD HAPS dynamic array expansion. 
C  20 Sep 01 AIR task 15,18: Comments 
C  11 Sep 01 AIR Task 18: Add HAPS output. 
C  10 Jul 01 AIR Task 15: New routine 
C 
C 
C  Local array subscripts: 
C 
C  NAMPOL(17)  -  NAMPOL ( IAT ) 
C  NAMVEH(3)   -  NAMVEH ( IVBUS) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DASH8    C*8  8 character string of hyphens 
C  NAMPOL   C*4  AT pollutant names 
C  NAMVEH   C*8  vehicle class names 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS  modified multispeed to a single speed. 
C                   SPD changed from an array to a single varibale. 
C 
      USE DATABASE, ONLY: DBSELTOX 
      USE ATHEAP, ONLY: EFAT, VAT, ATHIV, ATHSNAME 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I'  ! MAXIAT     
      INCLUDE 'ATOX5.I'  ! ATHSNAME 
      INCLUDE 'IOUCOM.I' ! IOUAT 
      INCLUDE 'IVTYPE.I' ! VTGASBUS,VTSCHOOL,VVDSL 
      INCLUDE 'VMXCOM.I' ! VMTMIX 
C 
      CHARACTER*10  CEF(3)      ! Character representation of EF 
      INTEGER       I           ! Loop index 
      INTEGER       IAT         ! Air Toxic pollutant index 
      INTEGER       IMAP        ! HAPS pollutant index 
      INTEGER       IN          ! Pollutant name index 
      INTEGER       IRB         ! Ratio/BER index 
      INTEGER       IV          ! Vehicle class 
      INTEGER       IVFIRST     ! First vehicle class 
      INTEGER       IVLAST      ! Last vehicle class 
      INTEGER       IXV         ! eXhaust/eVap index 
      INTEGER       JXV         ! HAPS eXhaust/eVap index 
      LOGICAL       PRNT        ! Printing togle 
      REAL          TOTAL(3)    ! Totals 
      REAL          XVTOT(3,2)  ! eXhaust/eVap totals 
C 
      CHARACTER(LEN=8) :: DASH8 = '  ------' 
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
      CHARACTER(LEN=8), DIMENSION(3) ::  NAMVEH=(/ 
     *   '  GasBUS','   URBAN','  SCHOOL'/) 
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
C  Determine output range 
C 
      IVFIRST = VTGASBUS 
      IVLAST  = VTSCHOOL 
C 
      WRITE(IOUAT,230) (NAMVEH(IV),IV=1,3), 
     *  (DASH8,IV=IVFIRST,IVLAST), 
     *  (VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',3(2X,A8), 
     *      /T21,3(2X,A8), 
     *      /T13,'VMT Mix:',3(2X,F8.4)) 
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
      IXV = 1 ! Exhaust only 
C 
      DO IAT=1,MAXIAT 
C 
        IF(DBSELTOX(IAT).EQ.2) THEN 
        WRITE(IOUAT,250) NAMPOL(IAT), 
     *    (EFAT(IXV,IAT,IV),IV=IVFIRST,IVLAST) 
  250   FORMAT(T3,A17,':',3F10.2) 
C 
      IF(IAT.LE.2) THEN 
       DO IV=IVFIRST,IVLAST 
        I=IV-IVFIRST+1 
        XVTOT(I,IAT)=XVTOT(I,IAT)+EFAT(IXV,IAT,IV) 
       END DO 
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
  251   FORMAT(T3,'Benzene Total Evp:',3F10.2) 
      ELSE 
      WRITE(IOUAT,252) TOTAL 
  252   FORMAT(T3,'   MTBE Total Evp:',3F10.2) 
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
     *WRITE(IOUAT,255) (XVTOT(IV-IVFIRST+1,1),IV=IVFIRST,IVLAST) 
  255   FORMAT(T3,'Benzene Exh + Evp:',39F10.2) 
      IF(DBSELTOX(2).EQ.2)  
     *WRITE(IOUAT,256) (XVTOT(IV-IVFIRST+1,2),IV=IVFIRST,IVLAST) 
  256   FORMAT(T3,'   MTBE Exh + Evp:',3F10.2) 
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
      DO IV=IVFIRST,IVLAST 
       IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
        PRNT = .TRUE. 
        WRITE(CEF(IV-IVFIRST+1),20) EFAT(IXV,IAT,IV) 
   20   FORMAT(ES10.3) 
       END IF 
      END DO 
C 
      END DO !IRB 
C 
      IF(PRNT) WRITE(IOUAT,232) TRIM(ATHSNAME(IMAP)),CETYPE(IXV),CEF 
  232 FORMAT(T3,A,' ',A,T20,':',3A10) 
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
        DO IV=IVFIRST,IVLAST 
        IF(VVDSL(IV).EQ.1) WRITE(CEF(IV-IVFIRST+1),20) EFAT(IXV,IAT,IV) 
        END DO 
C 
        DO IRB=1,4 
C 
         DO IV=IVFIRST,IVLAST 
          IF(ATHIV(JXV,IRB,IMAP,IV)) THEN 
           PRNT = .TRUE. 
           WRITE(CEF(IV-IVFIRST+1),20) EFAT(IXV,IAT,IV) 
          END IF 
        END DO 
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
   85 RETURN 
      END 
