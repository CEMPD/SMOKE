      SUBROUTINE PMOUTLDT() 
C 
C  OUTLDT writes out the expanded emission factor output table for 
C  LDTs, LDGT1, LDGT2, LDGT3, LDGT4, LDDT12 and LDDT34 - on unit 
C  IOUPM using descriptive output. 
C 
C  Called by PMOUTDSC 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  10 May 01 AIR Task 11/5: Added Idle PM. 
C  17 Apr 01 AIR Task 11: Changed code to always print all PM 
C            components, regardless of database selections. 
C  19 Feb 01 AIR Task 11: New Routine 
C 
C  Input on call: 
C 
C    argument list: None. 
C    modules: DATABASE 
C    common blocks: 
C    /IOUCOM/ IOUPM 
C    /MAXIMA/ MAXPOL 
C    /PART1/  MAXIPPM 
C    /PART3/  EFPM 
C    /VMXCOM/ VMTMIX 
C    /IDLCOM/ IDLEPMFLG 
C 
C  Local array subscripts: 
C 
C  NAMPOL(3)  -  NAMPOL ( IPPM ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C  TOTPM(6)   -  TOTPM  ( IPPM ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NAMPOL   C*4  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  DASH8     C*8  8 character string of hyphens 
C  EMISS    C*5  holds name of pollutant type for output 
C  MAPOUT   I    maps the PM output into the proper order 
C  TOTPM    R    temporary PM totals 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS modified multiple speed to  a single speed. 
C                   SPD changed from an array to a single variable. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      USE DATABASE, ONLY : DBSELPART 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IDLE.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER I 
      INTEGER IMAP 
      INTEGER IPPM 
      INTEGER IV 
      INTEGER, SAVE, DIMENSION(11) :: MAPOUT=(/5,4,3,2,1,8,9,6,7,10,11/) 
      REAL    TOTPM(6) 
C 
      CHARACTER(LEN=8) EMISS 
      CHARACTER(LEN=8), DIMENSION(6) :: NAMVEH=(/ 
     *   '  LDGT1 ', 
     *   '  LDGT2 ', 
     *   '  LDGT3 ', 
     *   '  LDGT4 ', 
     *   ' LDDT12 ', 
     *   ' LDDT34 '/) 
C 
      CHARACTER(LEN=9)::DASH8='  ------' 
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
C 
C  Initialize the total exhaust and total PM arrays. 
C 
      TOTPM=0. 
C 
      WRITE(IOUPM,230) NAMVEH,(DASH8,I=1,6), 
     *                  (VMTMIX(IV),IV=VTLDGT1,VTLDGT4), 
     *                  VMTMIX(VTLDDT12),VMTMIX(VTLDDT34) 
  230 FORMAT(T11,'Veh. Type: ',6(2X,A8),/T21,6(2X,A8), 
     *      /T13,'VMT Mix:',6(2X,F8.4)) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUPM,240) 
  240 FORMAT(1X,119('-'),/T2,'Composite Emission Factors (g/mi):') 
C 
      DO 30 IMAP=1,MAXIPPM 
        IPPM=MAPOUT(IMAP) 
C 
C     Skip output of EMFAC_PM 
C 
      IF(IPPM.EQ.10) CYCLE 
C 
C     Write output if selected. 
C 
        IF(DBSELPART(IPPM).EQ.2) THEN 
C 
          EMISS = NAMPOL(IPPM) 
C 
      SELECT CASE(IMAP) 
C 
      CASE (5,6,7,8,9,10) !S04,Tire,Brake,Tire,SO2,NH3,EMFAC_PM 
C 
        WRITE(IOUPM,250) EMISS, 
     *    (EFPM(IPPM,IV),IV=VTLDGT1,VTLDGT4), 
     *     EFPM(IPPM,VTLDDT12), 
     *     EFPM(IPPM,VTLDDT34) 
  250   FORMAT(T12,A8,':',6F10.4) 
C 
      CASE (1,2) !LEAD,GASPM 
C 
        WRITE(IOUPM,251) EMISS, 
     *    (EFPM(IPPM,IV),IV=VTLDGT1,VTLDGT4), 
     *    (DASH8,I=1,2) 
  251   FORMAT(T12,A8,':',4F10.4,2(2X,A8)) 
C 
      CASE (3,4) !ECARBON,OCARBON 
C 
        WRITE(IOUPM,252) EMISS, 
     *    (DASH8,I=VTLDGT1,VTLDGT4), 
     *     EFPM(IPPM,VTLDDT12), 
     *     EFPM(IPPM,VTLDDT34) 
  252   FORMAT(T12,A8,':',4(2X,A8),2F10.4) 
C 
      CASE (11) !Idle 
C 
        IF(IDLEPMFLG) WRITE(IOUPM,253) EMISS, 
     *    (DASH8,I=VTLDGT1,VTLDGT4), 
     *    (DASH8,I=1,2) 
  253   FORMAT(T2,'Idle Emissions (g/hr)'/ 
     *         T12,A8,':',6(2X,A8)) 
C 
      END SELECT 
C 
      ENDIF 
C 
C  Compute and write total exhaust and total PM values. 
C 
      IF(IPPM.LE.9) THEN 
        DO IV=VTLDGT1,VTLDGT4 
          TOTPM(IV-VTLDGT1+1)=TOTPM(IV-VTLDGT1+1)+EFPM(IPPM,IV) 
        END DO 
        TOTPM(5)=TOTPM(5)+EFPM(IPPM,VTLDDT12) 
        TOTPM(6)=TOTPM(6)+EFPM(IPPM,VTLDDT34) 
      ENDIF 
C 
      IF(IMAP.EQ.5) WRITE(IOUPM,254) TOTPM 
  254 FORMAT(T4,'Total Exhaust PM:',6F10.4) 
      IF(IMAP.EQ.7) WRITE(IOUPM,255) TOTPM 
  255 FORMAT(T12,'Total PM:',6F10.4) 
C 
   30 CONTINUE 
C 
      WRITE(IOUPM,256) 
  256 FORMAT(1X,119('-')) 
C 
   85 RETURN 
      END 
 
