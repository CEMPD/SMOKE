      SUBROUTINE OUTHDV(HDTYPE,ICY) 
C 
C  OUTHDV writes out the expanded emission factor output table for HDGV (if 
C  HDTYPE=1) or HDDV (if HDTYPE=2) on unit IOUREP using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  17 Aug 01 EPA bag  Expand the number of significant digits 
C  11 Oct 00 AIR Task 02: Bug Fix 197: Expand VMTMIX precision 
C  27 Jun 00 AIR Task 04: Added EXPANDED EVAP output. Commented out the 
C                        "natural" units evap code in case it needs to 
C                        be reinstated later. 
C  11 May 00 AIR Task 02: Changed Gm/mile to g/mi in format statement #240 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  16 Apr 00 AIR Task 02: Commented out evap messages. 
C  10 Apr 00 AIR Task 02: Replaced PRTFLG with DBSELPOL Removed OUTPOL calls. 
C  05 Apr 00 AIR Task 02: Rearranged and commented out the code which wrote 
C                         the Start and Running EFs. Since these values may 
C                         be output in the future, the code was left in. 
C  10 Nov 99 @DynTel-BG 1-006  added index parameter constants, IDIURNAL, IHOTSOAK, 
C            ICRNCASE for the evaporative emissions array, GREVP 
C  11 Mar 99 @DynTel-MLA 2-701 Modified to use new descriptive output format. 
C  16 Mar 99 @DynTel-bg 2-664 Removed evaporative temperature output, TEMEVP 
C  09 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Input on call: 
C 
C    argument list: HDTYPE 
C    common blocks: 
C    /CITRV1/ RVPICY 
C    /EVAPGR/ GREVP 
C    /FLAGS4/ NMHFLG,HCFLAG 
C    /IOUCOM/ IOUREP 
C    /MAXIMA/ MAXPOL 
C    /OPCNTL/ EXPEXH 
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL 
C    /RESUL2/ RLGGAL,RSTGPH,EFIDLE,G_MI_EVP 
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE 
C    /RESUL4/ RUNEF,STARTEF,TOTEVAP 
C    /TEMPS/  TEMEXH 
C    /VMXCOM/ VMTMIX 
C 
C  Local array subscripts: 
C 
C  NAMEVP(5)  -  NAMEVP ( IEVP ) 
C  NAMNMH(5)  -  NAMNMH ( NMHFLG ) 
C  NAMPAR(3)  -  NAMPAR ( IPARTL ) 
C  NAMPOL(3)  -  NAMPOL ( IP ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IV       I    vehicle loop index 
C  IP       I    pollutant loop index 
C  NAMNMH   C*5  prefix indicating type of Exh.-HC (total, non-meth., VOC, TOG) 
C  NAMPAR   C*9  partial diurnal cases 
C  NAMPOL   C*5  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  DASH8     C*8  8 character string of hyphens 
C  FMT1     I    format statement number assignment variable 
C  EMISS    C*5  holds name of pollutant type for output 
C  BUFFER1  C*80 holds character string for output 
C  BUFFER2  C*80 reformatted BUFFER1 from LJCOMP 
C  HEADER   C*80 message character string 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  HDTYPE    I   heavy duty vehicle classes 
C  IUNT      I   HCFLAG mapping for HEADER output 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS modified multispeed to single speed. SPD  
C                   changed from an array to a single variable. 
C 
      USE DATABASE, ONLY : DBSELPOL 
      IMPLICIT NONE 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RESUL2.I' 
      INCLUDE 'RESUL3.I' 
      INCLUDE 'RESUL4.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'VMXCOM.I' 
C 
CCC      CHARACTER*80  BUFFER1 
CCC      CHARACTER*80  BUFFER2 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: HDTYPE 
      INTEGER  FMT0,FMT1,FMT2,FMT3 
      INTEGER  IVFIRST, IVLAST 
      INTEGER  IP, IV 
CCC      INTEGER  IUNT 
      INTEGER, PARAMETER :: ICRNCASE=4,IRESTLOSS=5 
      INTEGER, PARAMETER :: IDIURNAL=2,IHOTSOAK=1 
C 
      CHARACTER(LEN=5) :: EMISS 
CCC      CHARACTER(LEN=80) :: HEADER 
      CHARACTER(LEN=8) :: DASH8='  ------' 
C 
      CHARACTER(LEN=5),DIMENSION(5)::NAMNMH=(/' THC ',' NMHC',' VOC ', 
     *            ' TOG ',' NMOG'/) 
C 
      CHARACTER(LEN=5),DIMENSION(3)::NAMPOL=(/' HC  ',' CO  ',' NOX '/) 
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
      WRITE(IOUREP,230) (NAMVEH(IV),IV=IVFIRST,IVLAST), 
     * (DASH8,IV=1,8),(VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',8(2X,A8),/T21,8(2X,A8), 
     *      /T13,'VMT Mix:',8(2X,F8.4)) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUREP,240) 
  240 FORMAT(118('-'),/T2,'Composite Emission Factors (g/mi):') 
C 
      DO 30 IP=1,MAXPOL 
       IF(DBSELPOL(IP).EQ.1) CYCLE 
C 
        IF (IP.EQ.1) THEN 
          EMISS = NAMNMH(NMHFLG) 
        ELSE 
          EMISS = NAMPOL(IP) 
        ENDIF 
C 
        ASSIGN 250 TO FMT0 
        IF(ICY>1999)  THEN 
           ASSIGN 251 TO FMT0 
          IF(IP==2.AND.HDTYPE==1) ASSIGN 252 TO FMT0 
        ENDIF 
        WRITE(IOUREP,FMT0) EMISS,(EFFTP(IP,IV),IV=IVFIRST,IVLAST) 
  250   FORMAT(T6,'Composite',A5,':',8F10.2) 
  251   FORMAT(T6,'Composite',A5,':',8F10.3) 
  252   FORMAT(T6,'Composite',A5,':',F9.2,7F10.2) 
C 
   30 CONTINUE 
C 
CC      WRITE(BUFFER1,252) NAMNMH(NMHFLG) 
CC  252 FORMAT('(Composite',A5,' emissions do not include', 
CC     *       ' refueling emissions.)') 
CC      CALL LJCOMP(BUFFER1,BUFFER2) 
CC      WRITE(IOUREP,253) BUFFER2 
CC  253 FORMAT(T2,A80,/118('-')) 
      WRITE(IOUREP,253) 
  253 FORMAT(118('-')) 
C 
C  Print EXPANDED EXHAUST, if necessary 
C 
      IF(EXPEXH) THEN 
C 
      WRITE(IOUREP,260) 
  260 FORMAT('Exhaust emissions (g/mi):') 
      IF (NMHFLG.EQ.2.OR.NMHFLG.EQ.5) THEN 
        ASSIGN 262 TO FMT1 
        IF(ICY>1999) ASSIGN 264 TO FMT1 
      ELSE 
        ASSIGN 263 TO FMT1 
        IF(ICY>1999) ASSIGN 265 TO FMT1 
      ENDIF 
  262 FORMAT( 
C     *       T9,A5,T15,'Start:',8F10.2/ 
C     *       T7,A5,T13,'Running:',8F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',8F10.2) 
  263 FORMAT( 
C     *       T8,A5,T15,'Start:',8F10.2/ 
C     *       T6,A5,T13,'Running:',8F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',8F10.2) 
  264 FORMAT( 
C     *       T9,A5,T15,'Start:',8F10.2/ 
C     *       T7,A5,T13,'Running:',8F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',8F10.3) 
  265 FORMAT( 
C     *       T8,A5,T15,'Start:',8F10.2/ 
C     *       T6,A5,T13,'Running:',8F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',8F10.3) 
      IF(DBSELPOL(1).EQ.2) WRITE(IOUREP,FMT1)  
C     * NAMNMH(NMHFLG),(STARTEF(1,IV),IV=IVFIRST,IVLAST), 
C     * NAMNMH(NMHFLG),(RUNEF(1,IV),IV=IVFIRST,IVLAST), 
     * NAMNMH(NMHFLG),(EFEXH(IV),IV=IVFIRST,IVLAST) 
C 
C 
      ASSIGN 270 TO FMT2 
      IF(ICY>1999)  THEN 
        ASSIGN 271 TO FMT2 
        IF(HDTYPE==1) ASSIGN 272 TO FMT2 
      ENDIF 
      IF(DBSELPOL(2).EQ.2) WRITE(IOUREP,FMT2) 
C     * (STARTEF(2,IV),IV=IVFIRST,IVLAST), 
C     * (RUNEF(2,IV),IV=IVFIRST,IVLAST), 
     * (EFFTP(2,IV),IV=IVFIRST,IVLAST) 
C 
  270 FORMAT( 
C     *        T12,'CO Start:',8F10.2/ 
C     *        T10,'CO Running:',8F10.2/ 
     *        T4,'CO Total Exhaust:',8F10.2) 
  271 FORMAT( 
C     *        T12,'CO Start:',8F10.2/ 
C     *        T10,'CO Running:',8F10.2/ 
     *        T4,'CO Total Exhaust:',8F10.3) 
  272 FORMAT( 
C     *        T12,'CO Start:',8F10.2/ 
C     *        T10,'CO Running:',8F10.2/ 
     *        T4,'CO Total Exhaust:',F9.2,7F10.2) 
C 
 
      ASSIGN 276 TO FMT3 
      IF(ICY>1999)  THEN 
        ASSIGN 277 TO FMT3 
      ENDIF 
      IF(DBSELPOL(3).EQ.2) WRITE(IOUREP,FMT3) 
C     * (STARTEF(3,IV),IV=IVFIRST,IVLAST), 
C     * (RUNEF(3,IV),IV=IVFIRST,IVLAST), 
     * (EFFTP(3,IV),IV=IVFIRST,IVLAST) 
C 
  276 FORMAT( 
C     *       T11,'NOx Start:',8F10.2/ 
C     *        T9,'NOx Running:',8F10.2/ 
     *        T3,'NOx Total Exhaust:',8F10.2) 
  277 FORMAT( 
C     *       T11,'NOx Start:',8F10.2/ 
C     *        T9,'NOx Running:',8F10.2/ 
     *        T3,'NOx Total Exhaust:',8F10.3) 
C 
      WRITE(IOUREP,282) 
  282 FORMAT(118('-')) 
C 
C  END EXPANDED EXHAUST SECTION 
      ENDIF 
C 
CCC  IN THE EVENT THAT "NATUALAL" EMISSIONS ARE REINSTATED, THE FOLLOWING 
CCC  CODE HAS BEEN COMMENTED OUT, RATHER THAN DELETED. 
CCC   50 IF(HCFLAG.LT.3.OR.DBSELPOL(2).EQ.2.OR.DBSELPOL(3).EQ.2) GOTO 85 
CCC      IF (HDTYPE.NE.1) GOTO 85 
CCC      IUNT=HCFLAG-2 
CCC      IF(IUNT.GT.2) IUNT=1 
CCC      IF (IUNT.EQ.1) THEN 
CCC        HEADER =  '(Hot Soak: g/trip; Diurnals: g; Crankcase:'// 
CCC     &            'g/mi.; Refuel: g/gal.; Resting: g/hr.' 
CCC      ELSE 
CCC        HEADER = '(All Components in Grams per Mile)' 
CCC      ENDIF 
CCC      WRITE(IOUREP,285) HEADER,RVPICY 
CCC  285 FORMAT('Non-Exhaust Emissions (g/mi):',T40,A80, 
CCC     *      /T55,'Weathered RVP:',F5.1/) 
C 
CCC      WRITE(IOUREP,290) (GREVP(IHOTSOAK,IV),IV=IVFIRST,IVLAST) 
CCC  290 FORMAT(T7,'Hot Soak Loss:',8F10.2) 
C 
CCC      IF(HCFLAG.GE.2) WRITE(IOUREP,295) 
CCC     *  (G_MI_EVP(IRESTLOSS,IV),IV=IVFIRST,IVLAST) 
CCC  295 FORMAT(T8,'Resting Loss:',8F10.2) 
C 
CCC      WRITE(IOUREP,300) (GREVP(ICRNCASE,IV),IV=IVFIRST,IVLAST) 
CCC  300 FORMAT(T6,'Crankcase Loss:',8F10.2) 
C 
CCC      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
CCC        WRITE(IOUREP,305) (REFL_GM_MI(IV),IV=IVFIRST,IVLAST), 
CCC     *    (EFRUNL(IV),IV=IVFIRST,IVLAST) 
CCC  305   FORMAT(T6,'Refueling Loss:',8F10.2, 
CCC     *      /T8,'Running Loss:',8F10.2) 
CCC      ELSE 
CCC        WRITE(IOUREP,310) (RLGGAL(IV),IV=IVFIRST,IVLAST), 
CCC     *    (RSTGPH(IV),IV=IVFIRST,IVLAST) 
CCC  310   FORMAT(T6,'Refueling Loss:',8F10.2, 
CCC     *        /T8,'Resting Loss:',8F10.2) 
CCC      ENDIF 
C 
CCC      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
CCC        WRITE(IOUREP,315) (GREVP(IDIURNAL,IV),IV=IVFIRST,IVLAST) 
CCC  315   FORMAT(T1,'Diurnal Loss (Vap.):',8F10.2) 
CCC      ELSE 
CCC        WRITE(IOUREP,316) (GREVP(IDIURNAL,IV),IV=IVFIRST,IVLAST) 
CCC  316   FORMAT(T6,'WtDiurnal Loss:',8F10.2) 
CCC      ENDIF 
C 
CCC      WRITE(IOUREP,320) (GREVP(3,IV),IV=IVFIRST,IVLAST) 
CCC  320 FORMAT(T1,'Diurnal Loss (Liq.):',8F10.2) 
C 
CCC      IF(HCFLAG.EQ.4) THEN 
CCC        DO 70 IDU=2,4 
CCC          IPARTL=IDU-1 
CCC          WRITE(IOUREP,325) NAMPAR(IPARTL), 
CCC     *      (GRPD(IPARTL,IV),IV=IVFIRST,IVLAST) 
CCC  325     FORMAT(T11,A9,':',8F10.2) 
CCC   70   CONTINUE 
CCC      ENDIF 
C 
CCC      IF (HCFLAG.GT.2) THEN 
CCC        WRITE(IOUREP,330) (TOTEVAP(IV),IV=IVFIRST,IVLAST) 
CCC  330   FORMAT(T4,'Total Evap. Loss:',8F10.2) 
CCC      ENDIF 
C 
CCC   80 WRITE(IOUREP,335) 
CCC  335 FORMAT(118('-'),/'Non-Exaust Emissions (alternate units):', 
CCC     *      /'        Diurnal (g):', 
CCC     *      /'     Resting (g/hr):', 
CCC     *      /'  Refueling (g/gal):',/118('-')) 
C 
C     Expanded evaporative/HC output, HCFLAG > 1 
C 
      IF(HCFLAG.GE.2 .AND. DBSELPOL(1).EQ.2.AND.HDTYPE.EQ.1) THEN 
C 
          WRITE(IOUREP,285) 
  285     FORMAT('Non-Exhaust Emissions (g/mi):') 
C 
          WRITE(IOUREP,290) (G_MI_EVP(IHOTSOAK,IV),IV=IVFIRST,IVLAST) 
  290     FORMAT(T7,'Hot Soak Loss:',8F10.3) 
C 
          WRITE(IOUREP,295) (G_MI_EVP(IDIURNAL,IV),IV=IVFIRST,IVLAST) 
  295     FORMAT(T8,'Diurnal Loss:',8F10.3) 
C 
          WRITE(IOUREP,300) (G_MI_EVP(IRESTLOSS,IV),IV=IVFIRST,IVLAST) 
  300     FORMAT(T8,'Resting Loss:',8F10.3) 
C 
          WRITE(IOUREP,315) (EFRUNL(IV),IV=IVFIRST,IVLAST) 
  315     FORMAT(T8,'Running Loss:',8F10.3) 
C 
          WRITE(IOUREP,305) (G_MI_EVP(ICRNCASE,IV),IV=IVFIRST,IVLAST) 
  305     FORMAT(T6,'Crankcase Loss:',8F10.3) 
C 
          WRITE(IOUREP,310) (REFL_GM_MI(IV),IV=IVFIRST,IVLAST) 
  310     FORMAT(T6,'Refueling Loss:',8F10.3) 
C 
        WRITE(IOUREP,320) (TOTEVAP(IV),IV=IVFIRST,IVLAST) 
  320   FORMAT(T3,'Total Non-Exhaust:',8F10.3/2X,118('-')) 
C 
      ENDIF 
C 
   85 RETURN 
      END 
