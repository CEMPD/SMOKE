      SUBROUTINE OUTBUS(ICY) 
C 
C  OUTBUS writes out the expanded emission factor output table for the BUS 
C  vehicles on unit IOUREP using output format 3. 
C 
C  Called by OUTDT3. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  03 Dec 01 AIR Task 28: Added CO2 and MPG output 
C  17 Aug 01 EPA bag  Expand the number of significant digits 
C  11 Oct 00 AIR Task 02: Bug Fix 197: Expand VMTMIX precision 
C  27 Jun 00 AIR Task 04: Added EXPANDED EVAP output. Commented out the 
C                        "natural" units evap code in case it needs to 
C                        be reinstated later. 
C  11 May 00 AIR Task 02: Changed Gm/mile to g/mi in format statement #240 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  15 Apr 00 AIR Task 02: Commented out evap messages. 
C  10 Apr 00: Replaced PRTFLG with DBSELPOL Removed OUTPOL calls. 
C  05 Apr 00 AIR Task 02: Rearranged and commented out the code which output the  
C                         Start and Running EFs. These values are not computed 
C                         separately now, but may be in the future. 
C  10 Nov 99 @DynTel-BG 1-006  added index parameter constants, IDIURNAL, IHOTSOAK, 
C                              ICRNCASE for the evaporative emissions array, GREVP 
C  11 Mar 99 @DynTel-MLA 2-701 Modified to use new descriptive output format. 
C  16 Mar 99 @DynTel-bag 2-664 Removed TEMEVP. Renamed multiple diurnals to 
C            diurnals from liquid leaks. 
C  02 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  22 Oct. 98 @Dyntel-HXQ 2-680 delete the vtcomml data. 
C 
C  Input on call: 
C 
C    argument list: ICY 
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
C    /VVTYPS/ VTGASBUS, VTSCHOOL  (IVTYPE.I) 
C 
C  Local array subscripts: 
C 
C  NAMEVP(5)  -  NAMEVP ( IEVP ) 
C  NAMNMH(5)  -  NAMNMH ( NMHFLG ) 
C  NAMPOL(3)  -  NAMPOL ( IP ) 
C  NAMVEH(4)  -  NAMVEH ( IVBUS ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IV       I    vehicle type index 
C  IP       I    pollutant loop variable 
C  NAMNMH   C*8  prefix indicating type of Exh.-HC (total, non-meth., VOC, TOG) 
C  NAMPOL   C*4  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  DASH8    C*8  8 character string of hyphens 
C  FMT1     I    format statement number assignment variable 
C  EMISS    C*5  holds name of pollutant type for output 
C  BUFFER1  C*80 holds character string for output 
C  BUFFER2  C*80 reformatted BUFFER1 from LJCOMP 
C  HEADER   C*80 message character string 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C 
C  Notes: 
C 
C  This routine was borrowed from OUTDT3. 
C  16 December 1997 Dyntel@MMS  modified multispeed to a single speed. 
C                   SPD changed from an array to a single varibale. 
C 
      USE DATABASE, ONLY: DBSELPOL 
      IMPLICIT NONE 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CO2COM.I' 
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
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
      INTEGER       FMT0,FMT1,FMT2,FMT3 
      INTEGER       IVFIRST, IVLAST, IV, IP 
CCC      INTEGER       IUNT 
      INTEGER, PARAMETER :: IDIURNAL=2,IHOTSOAK=1 
      INTEGER, PARAMETER :: ICRNCASE=4,IRESTLOSS=5 
      CHARACTER(LEN=5) :: EMISS 
CCC      CHARACTER(LEN=80) :: HEADER 
      CHARACTER(LEN=8) :: DASH8 = '  ------' 
      CHARACTER(LEN=4), DIMENSION(4) :: NAMPOL=(/ 
     *      ' HC  ',' CO  ',' NOX ',' CO2 '/) 
C 
      CHARACTER(LEN=5), DIMENSION(5) :: NAMNMH=(/ 
     * ' THC ',' NMHC',' VOC ', ' TOG ',' NMOG'/) 
C 
      CHARACTER(LEN=8), DIMENSION(3) ::  NAMVEH=(/ 
     *   '  GasBUS','   URBAN','  SCHOOL'/) 
C 
      IVFIRST = VTGASBUS 
      IVLAST  = VTSCHOOL 
C 
      WRITE(IOUREP,230) (NAMVEH(IV),IV=1,3), 
     *  (DASH8,IV=IVFIRST,IVLAST), 
     *  (VMTMIX(IV),IV=IVFIRST,IVLAST) 
  230 FORMAT(T11,'Veh. Type:',3(2X,A8), 
     *      /T21,3(2X,A8), 
     *      /T13,'VMT Mix:',3(2X,F8.4)) 
C 
C  If CO2 Is selected, print the harmonic mean MPG values. 
C 
      IF(DBSELPOL(4).EQ.2) WRITE(IOUREP,235)  
     *   (MPG(IV),IV=IVFIRST,IVLAST) 
  235 FORMAT(T2,'Fuel Economy (mpg):',10F10.1) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUREP,240) 
  240 FORMAT(118('-'),/T2,'Composite Emission Factors (g/mi):') 
C 
      DO 30 IP=1,MAXPOL+1 ! include CO2 
       IF(DBSELPOL(IP).EQ.1) CYCLE 
C 
        IF (IP.EQ.1) THEN 
          EMISS = NAMNMH(NMHFLG) 
        ELSE 
          EMISS = NAMPOL(IP) 
        ENDIF 
C 
        ASSIGN 250 TO FMT0 
        IF(ICY>1999) THEN 
          ASSIGN 251 TO FMT0 
          IF(IP==2) ASSIGN 252 TO FMT0 
        ENDIF 
C 
        IF(IP.EQ.4) THEN ! CO2 
        IF(ICY<2000) THEN 
        ASSIGN 253 TO FMT0  
        ELSE 
        ASSIGN 254 TO FMT0 
        END IF 
        END IF 
C 
        WRITE(IOUREP,FMT0) EMISS,(EFFTP(IP,IV),IV=IVFIRST,IVLAST) 
  250   FORMAT(T6,'Composite',A5,':',3F10.2) 
  251   FORMAT(T6,'Composite',A5,':',3F10.3) 
  252   FORMAT(T6,'Composite',A5,':',F9.2,1X,2F10.3) 
  253   FORMAT(T6,'Composite',A5,':',F9.1,2F10.1) 
  254   FORMAT(T6,'Composite',A5,':',F8.1,2F10.1) 
C 
   30 CONTINUE 
C 
CC      WRITE(BUFFER1,252) NAMNMH(NMHFLG) 
CC  252 FORMAT('(Composite',A5,' emissions do not include', 
CC     *       ' refueling emissions.)') 
CC      CALL LJCOMP(BUFFER1,BUFFER2) 
CC      WRITE(IOUREP,253) BUFFER2 
CC  253 FORMAT(T2,A80,/118('-')) 
      WRITE(IOUREP,255) 
  255 FORMAT(118('-')) 
C 
C  Print EXPANDED EXHAUST, if necessary 
C  NOTE: AS OF THIS VERSION, START AND RUNNING EMISSIONS 
C  FOR HDGV AND HDDV ARE NOT CALCULATED, SO BLANK OUT THE OUTPUT 
C  FOR NOW, 
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
C     *       T9,A5,T15,'Start:',3F10.2/ 
C     *       T7,A5,T13,'Running:',3F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',3F10.2) 
  263 FORMAT( 
C     *       T8,A5,T15,'Start:',3F10.2/ 
C     *       T6,A5,T13,'Running:',3F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',3F10.2) 
  264 FORMAT( 
C     *       T9,A5,T15,'Start:',3F10.2/ 
C     *       T7,A5,T13,'Running:',3F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',3F10.3) 
  265 FORMAT( 
C     *       T8,A5,T15,'Start:',3F10.2/ 
C     *       T6,A5,T13,'Running:',3F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',3F10.3) 
      IF(DBSELPOL(1).EQ.2) WRITE(IOUREP,FMT1)  
C     * NAMNMH(NMHFLG),(STARTEF(1,IV),IV=IVFIRST,IVLAST), 
C     * NAMNMH(NMHFLG),(RUNEF(1,IV),IV=IVFIRST,IVLAST), 
     * NAMNMH(NMHFLG),(EFEXH(IV),IV=IVFIRST,IVLAST) 
C 
      ASSIGN 270 TO FMT2 
      IF(ICY>1999) ASSIGN 271 TO FMT2 
      IF(DBSELPOL(2).EQ.2) WRITE(IOUREP,FMT2) 
C     * (STARTEF(2,IV),IV=IVFIRST,IVLAST), 
C     * (RUNEF(2,IV),IV=IVFIRST,IVLAST), 
     * (EFFTP(2,IV),IV=IVFIRST,IVLAST) 
C 
  270 FORMAT( 
C     *        T12,'CO Start:',3F10.2/ 
C     *        T10,'CO Running:',3F10.2/ 
     *        T4,'CO Total Exhaust:',3F10.2) 
  271 FORMAT( 
C     *        T12,'CO Start:',3F10.2/ 
C     *        T10,'CO Running:',3F10.2/ 
     *        T4,'CO Total Exhaust:',F9.2,1X,2F10.3) 
C 
      ASSIGN 272 TO FMT3 
      IF(ICY>1999) ASSIGN 273 TO FMT3 
      IF(DBSELPOL(3).EQ.2) WRITE(IOUREP,FMT3) 
C     * (STARTEF(3,IV),IV=IVFIRST,IVLAST), 
C     * (RUNEF(3,IV),IV=IVFIRST,IVLAST), 
     * (EFFTP(3,IV),IV=IVFIRST,IVLAST) 
C 
  272 FORMAT( 
C     *       T11,'NOx Start:',3F10.2/ 
C     *        T9,'NOx Running:',3F10.2/ 
     *        T3,'NOx Total Exhaust:',3F10.2) 
  273 FORMAT( 
C     *       T11,'NOx Start:',3F10.2/ 
C     *        T9,'NOx Running:',3F10.2/ 
     *        T3,'NOx Total Exhaust:',3F10.3) 
C 
      WRITE(IOUREP,282) 
  282 FORMAT(118('-')) 
C 
C  END OF EXPANDED EXHAUST 
      ENDIF 
C 
CCC  IN THE EVENT THAT "NATUALAL" EMISSIONS ARE REINSTATED, THE FOLLOWING 
CCC  CODE HAS BEEN COMMENTED OUT, RATHER THAN DELETED. 
ccc   50 IF(HCFLAG.LT.3.OR.DBSELPOL(2).EQ.2.OR.DBSELPOL(3).EQ.2) GOTO 85 
ccc      IUNT=HCFLAG-2 
ccc      IF(IUNT.GT.2) IUNT=1 
ccc      IF (IUNT.EQ.1) THEN 
ccc        HEADER =  '(Hot Soak: g/trip; Diurnals: g; Crankcase:'// 
ccc     &            'g/mi.; Refuel: g/gal.; Resting: g/hr.' 
ccc      ELSE 
ccc        HEADER = '(All Components in Grams per Mile)' 
ccc      ENDIF 
ccc      WRITE(IOUREP,285) HEADER,RVPICY 
ccc  285 FORMAT('Non-Exhaust Emissions (g/mi):',T40,A80, 
ccc     *      /T55,'Weathered RVP:',F5.1/) 
C 
ccc      WRITE(IOUREP,290) GREVP(IHOTSOAK,VTGASBUS) 
ccc  290 FORMAT(T7,'Hot Soak Loss:',F10.2) 
C 
ccc      IF(HCFLAG.GE.2) WRITE(IOUREP,295) G_MI_EVP(IRESTLOSS,VTGASBUS) 
ccc  295 FORMAT(T8,'Resting Loss:',F10.2) 
C 
ccc      WRITE(IOUREP,300) GREVP(ICRNCASE,VTGASBUS) 
ccc  300 FORMAT(T6,'Crankcase Loss:',F10.2) 
C 
ccc      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
ccc        WRITE(IOUREP,305) REFL_GM_MI(VTGASBUS),EFRUNL(VTGASBUS) 
ccc  305   FORMAT(T6,'Refueling Loss:',F10.2, 
ccc     *      /T8,'Running Loss:',F10.2) 
ccc      ELSE 
ccc        WRITE(IOUREP,310) RLGGAL(VTGASBUS),RSTGPH(VTGASBUS) 
ccc  310   FORMAT(T6,'Refueling Loss:',F10.2, 
ccc     *        /T8,'Resting Loss:',F10.2) 
ccc      ENDIF 
C 
ccc      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
ccc        WRITE(IOUREP,315) GREVP(IDIURNAL,VTGASBUS) 
ccc  315   FORMAT(T1,'Diurnal Loss (Vap.):',F10.2) 
ccc      ELSE 
ccc        WRITE(IOUREP,316) GREVP(IDIURNAL,VTGASBUS) 
ccc  316   FORMAT(T6,'WtDiurnal Loss:',F10.2) 
ccc      ENDIF 
C 
ccc      WRITE(IOUREP,320) GREVP(3,VTGASBUS) 
ccc  320 FORMAT(T1,'Diurnal Loss (Liq.):',F10.2) 
C 
CCC      IF(HCFLAG.EQ.4) THEN 
CCC        DO 70 IDU=2,4 
CCC          IPARTL=IDU-1 
CCC          WRITE(IOUREP,325) NAMPAR(IPARTL),GRPD(IPARTL,VTGASBUS) 
CCC  325     FORMAT(T11,A9,':',F10.2) 
CCC   70   CONTINUE 
CCC      ENDIF 
C 
ccc      IF (HCFLAG.GT.2) THEN 
ccc        WRITE(IOUREP,330) TOTEVAP(VTGASBUS) 
ccc  330   FORMAT(T4,'Total Evap. Loss:',F10.2) 
ccc      ENDIF 
C 
ccc   80 WRITE(IOUREP,335) 
ccc  335 FORMAT(118('-'),/'Non-Exaust Emissions (alternate units):', 
ccc     *      /'        Diurnal (g):', 
ccc     *      /'     Resting (g/hr):', 
ccc     *      /'  Refueling (g/gal):',/118('-')) 
C 
C     Expanded evaporative/HC output, HCFLAG > 1 
C 
      IF(HCFLAG.GE.2 .AND. DBSELPOL(1).EQ.2) THEN 
C 
          WRITE(IOUREP,285) 
  285     FORMAT('Non-Exhaust Emissions (g/mi):') 
C 
          WRITE(IOUREP,290) G_MI_EVP(IHOTSOAK,VTGASBUS) 
  290     FORMAT(T7,'Hot Soak Loss:',F10.3) 
C 
          WRITE(IOUREP,295) G_MI_EVP(IDIURNAL,VTGASBUS) 
  295     FORMAT(T8,'Diurnal Loss:',F10.3) 
C 
          WRITE(IOUREP,300) G_MI_EVP(IRESTLOSS,VTGASBUS) 
  300     FORMAT(T8,'Resting Loss:',F10.3) 
C 
          WRITE(IOUREP,315) EFRUNL(VTGASBUS) 
  315     FORMAT(T8,'Running Loss:',F10.3) 
C 
          WRITE(IOUREP,305) G_MI_EVP(ICRNCASE,VTGASBUS) 
  305     FORMAT(T6,'Crankcase Loss:',F10.3) 
C 
          WRITE(IOUREP,310) REFL_GM_MI(VTGASBUS) 
  310     FORMAT(T6,'Refueling Loss:',F10.3) 
C 
        WRITE(IOUREP,320) (TOTEVAP(IV),IV=IVFIRST,IVLAST) 
  320   FORMAT(T3,'Total Non-Exhaust:',3F10.3/2X,118('-')) 
C 
      ENDIF 
C 
   85 RETURN 
      END 
