      SUBROUTINE PMOUTHDV(HDTYPE) 
C 
C  OUTHDV writes out the expanded emission factor output table for HDGV (if 
C  HDTYPE=1) or HDDV (if HDTYPE=2) on unit IOUPM using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  17 Apr 01 AIR Task 11: Changed code to always print all PM 
C            components, regardless of database selections. 
C  19 Feb 01 AIR Task 11: New Routine 
C 
C  Called by PMOUTDSC 
C 
C  Calls None 
C 
C  Input on call: 
C 
C    argument list: HDTYPE 
C    modules: DATABASE 
C    common blocks: 
C    /IOUCOM/ IOUPM 
C    /VMXCOM/ VMTMIX 
C    /PART1 / MAXIPPM 
C    /PART3 / EFPM 
C    /IDLCOM/ IDLEPMFLG 
C 
C  Local array subscripts: 
C 
C  NAMPOL(3)  -  NAMPOL ( IPPM ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C  TOTPM(8)   -  TOTPM  ( IPPM ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DASH8     C*8  8 character string of hyphens 
C  EMISS    C*5  holds name of pollutant type for output 
C  HDTYPE    I   heavy duty vehicle classes 
C  IPPM       I  PM pollutant loop index 
C  IV       I    vehicle loop index 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  MAPOUT   I    maps the PM output into the proper order 
C  NAMPOL   C*5  PM pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  TOTPM    R    temporary PM totals 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS modified multispeed to single speed. SPD  
C                   changed from an array to a single variable. 
C 
      USE DATABASE, ONLY : DBSELPART 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IDLE.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER  HDTYPE 
      INTEGER  IMAP 
      INTEGER  IPPM 
      INTEGER  IV 
      INTEGER  IVFIRST 
      INTEGER  IVLAST 
      INTEGER, SAVE, DIMENSION(11) :: MAPOUT=(/5,4,3,2,1,8,9,6,7,10,11/) 
      REAL     TOTPM(8) 
C 
      CHARACTER(LEN=8) :: EMISS 
      CHARACTER(LEN=8) :: DASH8='  ------' 
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
C  Initialize the total exhaust and total PM arrays. 
C 
      TOTPM=0. 
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
      WRITE(IOUPM,230) (NAMVEH(IV),IV=IVFIRST,IVLAST), 
     * (DASH8,IV=1,8),(VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',8(2X,A8),/T21,8(2X,A8), 
     *      /T13,'VMT Mix:',8(2X,F8.4)) 
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
  250   FORMAT(T12,A8,':',8F10.4) 
C 
      CASE (1,2) !LEAD,GASPM 
C 
        IF(HDTYPE.EQ.1) THEN 
        WRITE(IOUPM,251) EMISS, 
     *    (EFPM(IPPM,IV),IV=IVFIRST,IVLAST) 
  251   FORMAT(T12,A8,':',8F10.4) 
        ELSE 
        WRITE(IOUPM,252) EMISS, 
     *    (DASH8,IV=IVFIRST,IVLAST) 
  252   FORMAT(T12,A8,':',8(2X,A8)) 
        ENDIF 
C 
      CASE (3,4) !ECARBON,OCARBON 
C 
        IF(HDTYPE.EQ.1) THEN 
        WRITE(IOUPM,252) EMISS, 
     *    (DASH8,IV=IVFIRST,IVLAST) 
        ELSE 
        WRITE(IOUPM,251) EMISS, 
     *    (EFPM(IPPM,IV),IV=IVFIRST,IVLAST) 
        ENDIF 
C 
      CASE (11) !Idle 
C 
        IF(IDLEPMFLG) THEN 
          IF(HDTYPE.EQ.1) THEN 
            WRITE(IOUPM,253) EMISS, 
     *      (DASH8,IV=IVFIRST,IVLAST) 
  253       FORMAT(T2,'Idle Emissions (g/hr)'/ 
     *         T12,A8,':',8(2X,A8)) 
          ELSE 
            WRITE(IOUPM,254) EMISS, 
     *      (EFPM(IPPM,IV),IV=IVFIRST,IVLAST) 
  254       FORMAT(T2,'Idle Emissions (g/hr)'/ 
     *         T12,A8,':',8F10.4) 
          END IF 
        END IF 
C 
      END SELECT 
C 
      ENDIF 
C 
C  Compute and write total exhaust and total PM values. 
C 
      IF(IPPM.LE.9) THEN 
        DO IV=IVFIRST,IVLAST 
          TOTPM(IV-IVFIRST+1)=TOTPM(IV-IVFIRST+1)+EFPM(IPPM,IV) 
        END DO 
      ENDIF 
C 
      IF(IMAP.EQ.5) WRITE(IOUPM,255) TOTPM 
  255 FORMAT(T4,'Total Exhaust PM:',8F10.4) 
      IF(IMAP.EQ.7) WRITE(IOUPM,256) TOTPM 
  256 FORMAT(T12,'Total PM:',8F10.4) 
C 
   30 CONTINUE 
C 
      WRITE(IOUPM,257) 
  257 FORMAT(1X,119('-')) 
C 
   85 RETURN 
      END 
