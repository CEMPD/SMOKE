      SUBROUTINE OUTLDT(ICY) 
C 
C  OUTLDT writes out the expanded emission factor output table for 
C  LDTs, LDGT1, LDGT2, LDGT3, LDGT4, LDDT12 and LDDT34 - on unit 
C  IOUREP using descriptive output. 
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
C  16 Apr 00 AIR Task 02: Commented out evap messages. 
C  13 Apr 00 AIR Task 02: Included LDDT12 and LDDT34 values in output. 
C  10 Apr 00 AIR Task 02: Replaced PRTFLG with DBSELPOL Removed OUTPOL calls. 
C  25 Feb 00 @DynTel-ddj 1-032, Expanded output to include the light 
C            duty diesel trucks (LDDT12 and LDDT34). 
C  13 Jan 00 @DynTel-ddj 1-013 Removed wetbulb, drybulb temperatures 
C            and towing and trailer parameters. 
C  10 Nov 99 @DynTel-BG 1-006  added index parameter constants, IDIURNAL, IHOTSOAK, 
C            ICRNCASE for the evaporative emissions array, GREVP 
C  22 Jun 99 @DynTel-MLA 2-694 Deleted obsolete LEV variables. 
C  11 Mar 99 @DynTel-MLA 2-701 Modified to use new descriptive output format. 
C  16 Mar 99 @DynTel-bag 2-664 Removed TEMEVP. Renamed multiple diurnals to 
C            diurnals from liquid leaks. 
C  10-Dec-98 @DynTel-ZK 2-000 Explicit Typing 
C  14-Sep-98 @Dyntel-bag 2-000 Added call to OUTZTT. 
C 
C  Input on call: 
C 
C    argument list: None. 
C    common blocks: 
C    /ALUHIN/ AC,XLOAD,TRAILR,ABSHUM,DB,WB 
C    /BASE10/ ZEVTOT 
C    /CITRV1/ RVPICY 
C    /EVAPGR/ GREVP,GRPD 
C    /FLAGS2/ ALHFLG 
C    /FLAGS4/ NMHFLG,HCFLAG 
C    /IOUCOM/ IOUREP 
C    /MAXIMA/ MAXPOL 
C    /OPCNTL/ EXPEXH 
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL 
C    /RESUL2/ RLGGAL,RSTGPH,EFIDLE,G_MI_EVP 
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE 
C    /RESUL4/ RUNEF,START3F,TOTEVAP 
C    /TEMPS/  TEMEXH 
C    /VMXCOM/ VMTMIX 
C 
C  Local array subscripts: 
C 
C  NAMNMH(5)  -  NAMNMH ( NMHFLG ) 
C  NAMPOL(3)  -  NAMPOL ( IP ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NAMNMH   C*8  prefix indicating type of Exh.-HC (total, non-meth., VOC, TOG) 
C  NAMPOL   C*4  pollutant names 
C  NAMVEH   C*8  vehicle class names 
C  DASH8     C*8  8 character string of hyphens 
C  FMT1     I    format statement number assignment variable 
C  EMISS    C*5  holds name of pollutant type for output 
C  BUFFER1  C*80 holds character string for output 
C  BUFFER2  C*80 reformatted BUFFER1 from LJCOMP 
C  HEADER   C*80 message character string 
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
      USE DATABASE, ONLY : DBSELPOL 
      IMPLICIT NONE 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'BASE10.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CO2COM.I' 
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'FLAGS2.I' 
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
 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER IP,I 
CCC   INTEGER IUNT 
C 
CCC      CHARACTER*80  BUFFER1 
CCC      CHARACTER*80  BUFFER2 
      INTEGER FMT0,FMT1,FMT2,FMT3 
      INTEGER IV 
      INTEGER, PARAMETER :: IDIURNAL=2,IHOTSOAK=1 
      INTEGER, PARAMETER :: ICRNCASE=4,IRESTLOSS=5 
CCC      CHARACTER*80  BUFFER1 
CCC      CHARACTER*80  BUFFER2 
CCC      CHARACTER(LEN=80) HEADER 
      CHARACTER(LEN=5) EMISS 
      CHARACTER(LEN=8), DIMENSION(6) :: NAMVEH=(/ 
     *   '  LDGT1 ', 
     *   '  LDGT2 ', 
     *   '  LDGT3 ', 
     *   '  LDGT4 ', 
     *   ' LDDT12 ', 
     *   ' LDDT34 '/) 
C 
      CHARACTER(LEN=9)::DASH8='  ------' 
      CHARACTER(LEN=5), DIMENSION(5) :: NAMNMH=( 
     *             /' THC ',' NMHC',' VOC ',' TOG ',' NMOG'/) 
      CHARACTER(LEN=5),DIMENSION(4)::NAMPOL=(/' HC  ',' CO  ',' NOX ', 
     *                                        ' CO2 '/) 
C 
      WRITE(IOUREP,230) NAMVEH,(DASH8,I=1,6),(VMTMIX(IV),IV=2,5), 
     *                  VMTMIX(15),VMTMIX(28) 
  230 FORMAT(T11,'Veh. Type: ',6(2X,A8),/T21,6(2X,A8), 
     *      /T13,'VMT Mix:',6(2X,F8.4)) 
C 
C  If CO2 Is selected, print the harmonic mean MPG values. 
C 
      IF(DBSELPOL(4).EQ.2) WRITE(IOUREP,235)  
     *   (MPG(IV),IV=2,5),MPG(15),MPG(28) 
  235 FORMAT(T2,'Fuel Economy (mpg):',10F10.1) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUREP,240) 
  240 FORMAT(118('-'),/T2,'Composite Emission Factors (g/mi):') 
C 
      DO 30 IP=1,MAXPOL+1   ! CO2 
       IF(DBSELPOL(IP).EQ.1) GOTO 30 
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
        WRITE(IOUREP,FMT0) EMISS,(EFFTP(IP,IV),IV=2,5), 
     *                    EFFTP(IP,15),EFFTP(IP,28) 
  250   FORMAT(T6,'Composite',A5,':',6F10.2) 
  251   FORMAT(T6,'Composite',A5,':',6F10.3) 
  252   FORMAT(T6,'Composite',A5,':',F9.2,3F10.2,1X,2F10.3) 
  253   FORMAT(T6,'Composite',A5,':',6(F9.1,1X)) 
  254   FORMAT(T6,'Composite',A5,':',F8.1,3F10.1,2F10.1) 
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
C 
      IF(EXPEXH) THEN 
C 
      WRITE(IOUREP,260) 
  260 FORMAT('Exhaust emissions (g/mi):') 
      IF (NMHFLG.EQ.2.OR.NMHFLG.EQ.5) THEN 
        ASSIGN 262 TO FMT1 
        IF(ICY>1999)  ASSIGN 264 TO FMT1 
      ELSE 
        ASSIGN 263 TO FMT1 
        IF(ICY>1999)  ASSIGN 265 TO FMT1 
      ENDIF 
  262 FORMAT(T9,A5,T15,'Start:',6F10.2/ 
     *       T7,A5,T13,'Running:',6F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',6F10.2) 
  263 FORMAT(T8,A5,T15,'Start:',6F10.2/ 
     *       T6,A5,T13,'Running:',6F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',6F10.2) 
  264 FORMAT(T9,A5,T15,'Start:',6F10.3/ 
     *       T7,A5,T13,'Running:',6F10.3/ 
     *       T1,A5,T7 ,'Total Exhaust:',6F10.3) 
  265 FORMAT(T8,A5,T15,'Start:',6F10.3/ 
     *       T6,A5,T13,'Running:',6F10.3/ 
     *       T2,A5,T7 ,'Total Exhaust:',6F10.3) 
      IF(DBSELPOL(1).EQ.2) WRITE(IOUREP,FMT1)  
     * NAMNMH(NMHFLG), 
     *(STARTEF(1,IV),IV=2,5),STARTEF(1,15),STARTEF(1,28), 
     * NAMNMH(NMHFLG), 
     *(RUNEF(1,IV),IV=2,5),RUNEF(1,15),RUNEF(1,28), 
     * NAMNMH(NMHFLG), 
     *(EFEXH(IV),IV=2,5),EFEXH(15),EFEXH(28) 
C 
        ASSIGN 270 TO FMT2 
        IF(ICY>1999)  ASSIGN 271 TO FMT2 
      IF(DBSELPOL(2).EQ.2) WRITE(IOUREP,FMT2) 
     * (STARTEF(2,IV),IV=2,5),STARTEF(2,15),STARTEF(2,28), 
     * (RUNEF(2,IV),IV=2,5),RUNEF(2,15),RUNEF(2,28), 
     * (EFFTP(2,IV),IV=2,5),EFFTP(2,15),EFFTP(2,28) 
C 
  270 FORMAT(/T12,'CO Start:',6F10.2/ 
     *        T10,'CO Running:',6F10.2/ 
     *        T4,'CO Total Exhaust:',6F10.2) 
  271 FORMAT(/T12,'CO Start:',F9.2,3F10.2,1X,2F10.3/ 
     *        T10,'CO Running:',F9.2,3F10.2,1X,2F10.3/ 
     *        T4,'CO Total Exhaust:',F9.2,3F10.2,1X,2F10.3) 
C 
        ASSIGN 272 TO FMT3 
        IF(ICY>1999)  ASSIGN 273 TO FMT3 
      IF(DBSELPOL(3).EQ.2) WRITE(IOUREP,FMT3) 
     * (STARTEF(3,IV),IV=2,5),STARTEF(3,15),STARTEF(3,28), 
     * (RUNEF(3,IV),IV=2,5),RUNEF(3,15),RUNEF(3,28), 
     * (EFFTP(3,IV),IV=2,5),EFFTP(3,15),EFFTP(3,28) 
C 
  272 FORMAT(/T11,'NOx Start:',6F10.2/ 
     *        T9,'NOx Running:',6F10.2/ 
     *        T3,'NOx Total Exhaust:',6F10.2) 
  273 FORMAT(/T11,'NOx Start:',6F10.3/ 
     *        T9,'NOx Running:',6F10.3/ 
     *        T3,'NOx Total Exhaust:',6F10.3) 
C 
      WRITE(IOUREP,282) 
  282 FORMAT(118('-')) 
C 
C  END OF EXPANDED EXHAUST SECTION 
      ENDIF 
C 
CCC  IN THE EVENT THAT "NATUALAL" EMISSIONS ARE REINSTATED, THE FOLLOWING 
CCC  CODE HAS BEEN COMMENTED OUT, RATHER THAN DELETED. 
CCC   50 IF(HCFLAG.LT.3.OR.DBSELPOL(2).EQ.2.OR.DBSELPOL(3).EQ.2) GOTO 85 
CCC      IUNT=HCFLAG-2 
CCC      IF(IUNT.GT.2) IUNT=1 
CCC      IF (IUNT.EQ.1) THEN 
CCC       HEADER =  '(Hot Soak: g/trip; Diurnals: g; Crankcase:'// 
CCC     &            'g/mi.; Refuel: g/gal.; Resting: g/hr.' 
CCC      ELSE 
CCC        HEADER = '(All Components in Grams per Mile)' 
CCC      ENDIF 
CCC      WRITE(IOUREP,285) HEADER,RVPICY 
CCC  285 FORMAT('Non-Exhaust Emissions (g/mi):',T40,A80, 
CCC     *      /T55,'Weathered RVP:',F5.1/) 
C 
CCC      WRITE(IOUREP,290) (GREVP(IHOTSOAK,IV),IV=2,5) 
CCC  290 FORMAT(T7,'Hot Soak Loss:',4F10.2) 
C 
CCC      IF(HCFLAG.GE.2) WRITE(IOUREP,295) (G_MI_EVP(IRESTLOSS,IV),IV=2,5) 
CCC  295 FORMAT(T8,'Resting Loss:',4F10.2) 
C 
CCC      WRITE(IOUREP,300) (GREVP(ICRNCASE,IV),IV=2,5) 
CCC  300 FORMAT(T6,'Crankcase Loss:',4F10.2) 
C 
CCC      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
CCC        WRITE(IOUREP,305) (REFL_GM_MI(IV),IV=2,5),(EFRUNL(IV),IV=2,5) 
CCC  305   FORMAT(T6,'Refueling Loss:',4F10.2, 
CCC     *      /T8,'Running Loss:',4F10.2) 
CCC      ELSE 
CCC        WRITE(IOUREP,310) (RLGGAL(IV),IV=2,5),(RSTGPH(IV),IV=2,5) 
CCC  310   FORMAT(T6,'Refueling Loss:',4F10.2, 
CCC     *        /T8,'Resting Loss:',4F10.2) 
CCC      ENDIF 
C 
CCC      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5) THEN 
CCC        WRITE(IOUREP,315) (GREVP(IDIURNAL,IV),IV=2,5) 
CCC  315   FORMAT(T1,'Diurnal Loss (Vap.):',4F10.2) 
CCC      ELSE 
CCC        WRITE(IOUREP,316) (GREVP(IDIURNAL,IV),IV=2,5) 
CCC  316   FORMAT(T6,'WtDiurnal Loss:',4F10.2) 
CCC      ENDIF 
C 
CCC      WRITE(IOUREP,320) (GREVP(3,IV),IV=2,5) 
CCC  320 FORMAT(T1,'Diurnal Loss (Liq.):',4F10.2) 
C 
CCC      IF (HCFLAG.GT.2) THEN 
CCC        WRITE(IOUREP,330) (TOTEVAP(IV),IV=2,5) 
CCC  330   FORMAT(T4,'Total Evap. Loss:',4F10.2) 
CCC      ENDIF 
C 
CCC   80 WRITE(IOUREP,335) 
CCC  335 FORMAT(118('-'),/'Non-Exaust Emissions (alternate units):', 
CCC     *      /'        Diurnal (g):', 
CCC     *      /'     Resting (g/hr):', 
CCC     *      /'  Refueling (g/gal):',/118('-')) 
C 
C 
C     Expanded evaporative/HC output, HCFLAG > 1 
C 
      IF(HCFLAG.GE.2 .AND. DBSELPOL(1).EQ.2) THEN 
C 
          WRITE(IOUREP,285) 
  285     FORMAT('Non-Exhaust Emissions (g/mi):') 
C 
          WRITE(IOUREP,290) (G_MI_EVP(IHOTSOAK,IV),IV=2,5) 
  290     FORMAT(T7,'Hot Soak Loss:',4F10.3) 
C 
          WRITE(IOUREP,295) (G_MI_EVP(IDIURNAL,IV),IV=2,5) 
  295     FORMAT(T8,'Diurnal Loss:',4F10.3) 
C 
          WRITE(IOUREP,300) (G_MI_EVP(IRESTLOSS,IV),IV=2,5) 
  300     FORMAT(T8,'Resting Loss:',4F10.3) 
C 
          WRITE(IOUREP,315) (EFRUNL(IV),IV=2,5) 
  315     FORMAT(T8,'Running Loss:',4F10.3) 
C 
          WRITE(IOUREP,305) (G_MI_EVP(ICRNCASE,IV),IV=2,5) 
  305     FORMAT(T6,'Crankcase Loss:',4F10.3) 
C 
          WRITE(IOUREP,310) (REFL_GM_MI(IV),IV=2,5) 
  310     FORMAT(T6,'Refueling Loss:',4F10.3) 
C 
        WRITE(IOUREP,320) (TOTEVAP(IV),IV=2,5), 
     *   TOTEVAP(15),TOTEVAP(28) 
  320   FORMAT(T3,'Total Non-Exhaust:',6F10.3/2X,118('-')) 
C 
      ENDIF 
   85 RETURN 
      END 
 
