      SUBROUTINE PMOUTBUS() 
C 
C  OUTBUS writes out the expanded emission factor output table for the BUS 
C  vehicles on unit IOUPM using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  10 May 01 AIR Task 11/5: added Idle PM. 
C  17 Apr 01 AIR Task 11: Changed code to always print all PM 
C            components, regardless of database selections. 
C  19 Feb 01 AIR Task 11: New routine 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    module: DATABASE 
C    common blocks: 
C    /IOUCOM/ IOUPM 
C    /PART1/ MAXIPPM 
C    /PART3/ EFPM,VPM 
C    /VMXCOM/ VMTMIX 
C    /VVTYPS/ VTGASBUS, VTSCHOOL  (IVTYPE.I) 
C    /IDLCOM/ IDLEPMFLG 
C 
C  Local array subscripts: 
C 
C  NAMPOL(3)  -  NAMPOL ( IPPM ) 
C  NAMVEH(4)  -  NAMVEH ( IVBUS) 
C  TOTPM(3)   -  TOTPM  ( IPPM ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DASH8    C*8  8 character string of hyphens 
C  EMISS    C*5  holds name of pollutant type for output 
C  IPPM      I   PM pollutant loop variable 
C  IV        I   vehicle type index 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  MAPOUT    I   maps the PM output into the proper order 
C  NAMPOL   C*4  PM pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  TOTPM     R   temporary PM totals 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS  modified multispeed to a single speed. 
C                   SPD changed from an array to a single varibale. 
C 
      USE DATABASE, ONLY: DBSELPART 
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
      INTEGER       IMAP 
      INTEGER       IPPM 
      INTEGER       IV 
      INTEGER       IVFIRST 
      INTEGER       IVLAST 
      INTEGER,SAVE, DIMENSION(11) :: MAPOUT=(/5,4,3,2,1,8,9,6,7,10,11/) 
      REAL          TOTPM(3) 
C 
      CHARACTER(LEN=8) :: EMISS 
      CHARACTER(LEN=8) :: DASH8 = '  ------' 
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
      CHARACTER(LEN=8), DIMENSION(3) ::  NAMVEH=(/ 
     *   '  GasBUS','   URBAN','  SCHOOL'/) 
C 
C 
C  Initialize the total exhaust and total PM arrays. 
C 
      TOTPM=0. 
C 
C  Determine output range 
C 
      IVFIRST = VTGASBUS 
      IVLAST  = VTSCHOOL 
C 
      WRITE(IOUPM,230) (NAMVEH(IV),IV=1,3), 
     *  (DASH8,IV=IVFIRST,IVLAST), 
     *  (VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',3(2X,A8), 
     *      /T21,3(2X,A8), 
     *      /T13,'VMT Mix:',3(2X,F8.4)) 
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
     *    (EFPM(IPPM,IV),IV=IVFIRST,IVLAST) 
  250   FORMAT(T12,A8,':',3F10.4) 
C 
      CASE (1,2) !LEAD,GASPM 
C 
        WRITE(IOUPM,251) EMISS, 
     *    EFPM(IPPM,VTGASBUS), 
     *    (DASH8,IV=VTURBAN,VTSCHOOL)    
  251   FORMAT(T12,A8,':',F10.4,2(2X,A8)) 
C 
      CASE (3,4) !ECARBON,OCARBON 
C 
        WRITE(IOUPM,252) EMISS, 
     *    DASH8, 
     *    (EFPM(IPPM,IV),IV=VTURBAN,VTSCHOOL) 
  252   FORMAT(T12,A8,':',(2X,A8),2F10.4) 
C 
      CASE (11) !Idle 
C 
        IF(IDLEPMFLG) WRITE(IOUPM,253) EMISS, 
     *    DASH8, 
     *    (EFPM(IPPM,IV),IV=VTURBAN,VTSCHOOL) 
  253   FORMAT(T2,'Idle Emissions (g/hr)'/ 
     *         T12,A8,':',(2X,A8),2F10.4) 
C 
      END SELECT 
C 
      ENDIF 
C 
C  Compute and write total exhaust and total PM values. 
C  Never include EMFAC_PM or Idle. 
C 
      IF(IPPM.LE.9) THEN 
        DO IV=IVFIRST,IVLAST 
          TOTPM(IV-IVFIRST+1)=TOTPM(IV-IVFIRST+1)+EFPM(IPPM,IV) 
        END DO 
      ENDIF 
C 
      IF(IMAP.EQ.5) WRITE(IOUPM,254) TOTPM 
  254 FORMAT(T4,'Total Exhaust PM:',3F10.4) 
      IF(IMAP.EQ.7) WRITE(IOUPM,255) TOTPM 
  255 FORMAT(T12,'Total PM:',3F10.4) 
C 
   30 CONTINUE 
C 
      WRITE(IOUPM,256) 
  256 FORMAT(1X,119('-')) 
C 
   85 RETURN 
      END 
