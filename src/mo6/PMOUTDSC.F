      SUBROUTINE PMOUTDSC(ICY) 
C 
C  OUTDSC writes descriptive output for scenario-level data. 
C 
C  Called by OUTPUT. 
C 
C  Calls EFCVT, OUTBUS, OUTDTG, OUTHD3, OUTHDR, OUTHDV, and OUTLDT. 
C 
C  Changes: (Last change first) 
C 
C  17 Apr 01 AIR Task 11: Changed code to always print all PM 
C            components, regardless of database selections. 
C  19 Feb 01 AIR Task 11: New routine. 
C 
C  Called by PMCALX 
C 
C  Calls PMOUTLDT, PMOUTHDV, PMOUTBUS 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    modules : None 
C    common blocks: 
C    /CITPAR/     SCNAME 
C    /CVNAMS/     ILDGT12,ILDGT34,IHDGV,ILDDV,ILDDT,IHDDV,ILDGT (IVTYPE.I) 
C    /EVAL/       MONTH, MEVAL 
C    /IOUCOM/     IOUPM 
C    /PART1/      MAXIPPM 
C    /PART3/      EFPM 
C    /RFORM1/     RFGFLG 
C    /STRING/     NOYES,COMMA,PERIOD 
C    /SULFUR/     S_PPM 
C    /VMXCOM/     VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTHDD 
C    /YEARS4/     IY1994 
C 
C  Local array subscripts: 
C 
C  NAMPOL(8)  -  NAMPOL ( IPPM ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C  TOTPM(10)  -  TOTPM  ( IPPM ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DASH8    C*8  8 character string of hyphens 
C  EMISS    C*5  holds name of pollutant type for output 
C  IRFG      I   index for NOYES(IRFG) refueling loss 
C  ITER      I   write statement index 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  MAPOUT    I   maps the PM output into the proper order 
C  NAMPOL   C*4  pollutant names 
C  NAMUNT   C*12 component units, HCFLAG=3,5: mixed, HCFLAG=4: g/mi 
C  NAMVEH   C*8  vehicle class names 
C  SUM       R   summation variable 
C  TOTPM     R   temporary PM totals 
C 
C 
      USE DATABASE, ONLY: DBSELPART 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IDLE.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'STRING.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'VMXCOM.I' 
      INCLUDE 'WEEKLY.I' 
      INCLUDE 'YEARS4.I' 
C 
C     Arguments: 
C 
      INTEGER, INTENT(IN) :: ICY 
C 
C     External functions: 
C 
      REAL EFCVT  
C 
C     Local variables:                          
C 
      INTEGER IMAP 
      INTEGER IPPM 
      INTEGER IRFG 
      INTEGER ITER 
      INTEGER, SAVE, DIMENSION(11) :: MAPOUT=(/5,4,3,2,1,8,9,6,7,10,11/) 
      REAL    SUM 
      REAL    TOTPM(10) 
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
      CHARACTER(LEN=9) :: DASH8='  ------' 
      CHARACTER(LEN=8) :: EMISS 
C 
C 
C  If no descriptive output selected, return. 
C 
      IF(NODESC) RETURN 
C 
C  Initialize the total exhaust and total PM arrays. 
C 
      TOTPM=0. 
C 
      WRITE(IOUPM, 400) 
  400 FORMAT() 
C 
C  Determine reformulated gas reponse. 
C 
      IF (RFGFLG.EQ.1) THEN 
        IRFG = 1 
      ELSE 
        IRFG = 2 
      ENDIF 
C 
C  Write header 
C 
   10 WRITE(IOUPM,210) ICY,MONTH(MEVAL),S_PPM,D_PPM,PSZCUT,NOYES(IRFG) 
  210 FORMAT(T15,'                Calendar Year: ',I5/ 
     *       T15,'                        Month: ',1X,A4/ 
     *       T15,' Gasoline Fuel Sulfur Content: ',F5.0,' ppm'/ 
     *       T15,'   Diesel Fuel Sulfur Content: ',F5.0,' ppm'/ 
     *       T15,'         Particle Size Cutoff: ',F5.2,' Microns'/ 
     *       T15,'             Reformulated Gas: ',1X,A4/) 
C 
      IF( IWEEK .EQ. 2 ) THEN 
        WRITE (IOUPM,217) 
      END IF 
  217 FORMAT(T1,'Emissions determined from WEEKEND hourly ', 
     *'vehicle activity fractions.'/) 
C 
C  Write VMTMIX summary 
C 
      SUM=VMTMIX(VTLDGV)+VMTGT12+VMTGT34+VMTHDG 
     *   +VMTMIX(VTLDDV)+VMTLDDT+VMTHDD+VMTMIX(VTMC) 
C 
      WRITE(IOUPM,230) 
     *   NAMVEH,(DASH8,ITER=1,10), 
     *   VMTMIX(VTLDGV),VMTGT12,VMTGT34,VMTHDG, 
     *   VMTMIX(VTLDDV),VMTLDDT,VMTHDD,VMTMIX(VTMC),SUM 
  230 FORMAT(T8,'Vehicle Type:',10(2X,A8)/ 
     *       T16,'GVWR:',T36,'<6000',T46,'>6000',T56,'(All)'/ 
     *       T21,10(2X,A8)/ 
     *       T4,'VMT Distribution:',3F10.4,10X,6F10.4) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUPM,240) 
  240 FORMAT(1X,119('-'),/T2, 
     *   'Composite Emission Factors (g/mi):') 
C 
      DO IMAP=1,MAXIPPM 
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
     *    EFPM(IPPM,VTLDGV), 
     *    EFCVT(ILDGT12,IEFPM,IPPM), 
     *    EFCVT(ILDGT34,IEFPM,IPPM), 
     *    EFCVT(ILDGT,IEFPM,IPPM), 
     *    EFPM(IPPM,IHDGV), 
     *    EFPM(IPPM,VTLDDV), 
     *    EFCVT(ILDDT,IEFPM,IPPM), 
     *    EFCVT(IHDDV,IEFPM,IPPM), 
     *    EFPM(IPPM,VTMC), 
     *    VPM(IPPM) 
  250   FORMAT(T12,A8,':',9F10.4,F10.4) 
C 
      CASE (1,2) !LEAD,GASPM 
C 
          WRITE(IOUPM,251) EMISS, 
     *    EFPM(IPPM,VTLDGV), 
     *    EFCVT(ILDGT12,IEFPM,IPPM), 
     *    EFCVT(ILDGT34,IEFPM,IPPM), 
     *    EFCVT(ILDGT,IEFPM,IPPM), 
     *    EFPM(IPPM,IHDGV), 
     *    (DASH8,ITER=1,3), 
     *    EFPM(IPPM,VTMC), 
     *    VPM(IPPM) 
  251   FORMAT(T12,A8,':',5F10.4,3(2X,A8),2F10.4) 
C 
      CASE (3,4) !ECARBON,OCARBON 
C 
          WRITE(IOUPM,252) EMISS, 
     *    (DASH8,ITER=1,5), 
     *    EFPM(IPPM,VTLDDV), 
     *    EFCVT(ILDDT,IEFPM,IPPM), 
     *    EFCVT(IHDDV,IEFPM,IPPM), 
     *    DASH8, 
     *    VPM(IPPM) 
  252   FORMAT(T12,A8,':',5(2X,A8),3F10.4,(2X,A8),F10.4) 
C 
      CASE (11) !Idle 
C 
        IF(IDLEPMFLG) WRITE(IOUPM,253) EMISS, 
     *    (DASH8,ITER=1,7), 
     *    EFCVT(IHDDV,IEFPM,IPPM), 
     *    DASH8, 
     *    VPM(IPPM) 
  253   FORMAT(T2,'Idle Emissions (g/hr)'/ 
     *         T12,A8,':',7(2X,A8),F10.4,2X,A8,F10.4) 
C 
      END SELECT 
C 
      ENDIF 
C 
C  Compute and write total exhaust and total PM values. 
C 
      IF(IPPM.LE.9) THEN 
        TOTPM( 1)=TOTPM( 1)+EFPM(IPPM,VTLDGV) 
        TOTPM( 2)=TOTPM( 2)+EFCVT(ILDGT12,IEFPM,IPPM) 
        TOTPM( 3)=TOTPM( 3)+EFCVT(ILDGT34,IEFPM,IPPM) 
        TOTPM( 4)=TOTPM( 4)+EFCVT(ILDGT,IEFPM,IPPM) 
        TOTPM( 5)=TOTPM( 5)+EFPM(IPPM,IHDGV) 
        TOTPM( 6)=TOTPM( 6)+EFPM(IPPM,VTLDDV) 
        TOTPM( 7)=TOTPM( 7)+EFCVT(ILDDT,IEFPM,IPPM) 
        TOTPM( 8)=TOTPM( 8)+EFCVT(IHDDV,IEFPM,IPPM) 
        TOTPM( 9)=TOTPM( 9)+EFPM(IPPM,VTMC) 
        TOTPM(10)=TOTPM(10)+VPM(IPPM) 
      ENDIF 
C 
      IF(IMAP.EQ.5) WRITE(IOUPM,254) TOTPM 
  254 FORMAT(T4,'Total Exhaust PM:',10F10.4) 
      IF(IMAP.EQ.7) WRITE(IOUPM,255) TOTPM 
  255 FORMAT(T12,'Total PM:',10F10.4) 
C 
      END DO 
C 
      WRITE(IOUPM,256) 
  256 FORMAT(1X,119('-')) 
C 
C  Check the expanded output requests. 
C 
   85 IF (OPTABL(1).EQ.2) CALL PMOUTLDT 
      IF (OPTABL(2).EQ.2) CALL PMOUTHDV(1) 
      IF (OPTABL(3).EQ.2) CALL PMOUTHDV(2) 
      IF (OPTABL(4).EQ.2) CALL PMOUTBUS 
C 
      RETURN 
      END 
