      SUBROUTINE OUTDSC(ICY) 
C 
C  OUTDSC writes descriptive output for scenario-level data. 
C 
C  Called by OUTPUT. 
C 
C  Calls EFCVT, OUTBUS, OUTDTG, OUTHD3, OUTHDR, OUTHDV, and OUTLDT. 
C 
C  Changes: (Last change first) 
C 
C  27 Nov 01 AIR Task 27: Added CO2 and MPG output 
C  17 Aug 01 EPA bag  Expand the number of significant digits 
C  23 Oct 00 @EPA-djb Bug288 Added decimal to temperture output. 
C  06 Oct 00 @EPA-bag BUG155 Remove FTP_FLAG 
C  20 Sep 00 AIR Task 07: Bug Fix 197: insufficient VMT mix precision, 
C            Bug Fix 192: Changed AC output 
C  15 Sep 00 AIR Task 08: Updated output to display Exh and Evap I/M separately 
C  21 Jul 00 @EPA-RAG running loss; removed pass/fail fractions, PFRATE's, 
C            for purge/pressure evap. tests 
C  27 JUN 00 AIR Task 04: Blanked out the start and running emissions from 
C            the All Veh columns. Commented out the alternate evap units 
C            outputs. Removed the by-etype values from All Veh. 
C  16 May 00 @EPA-ddj BUG87, Removed affect of STARFG = GPS. 
C  11 May 00 AIR Task 02: Changed mileage sums 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  15 Apr 00 AIR Task 02: Commented out evap messages 
C  10 Apr 00 AIR Task 02: Replaced PRTFLG with DBSELPOL Removed OUTPOL calls. 
C  03 Apr 00 AIR Task 02: Removed output of ZEV fractions, blanked out 
C            start emissions for HDGV and HDDV, put in total vmtmix, coded 
C            the mandatory output section, added EXPANDED EXHAUST code. 
C  20 Feb 00 @DynTel-BG 1-005  added notification for liquid leaks not bing included 
C            in the evaporative (and exhaust) emissions sums and calculations 
C   3 Feb 00 @DynTel-ddj 1-032, Modified output for LDDTs.  LDDT12 
C            and LDDT34 combined as a composite vehicle group. 
C  13 Jan 00 @DynTel-ddj 1-013, Removed wetbulb and drybulb temperatues 
C            and towing and trailer parameters. Added provisions for 
C            AC parameters in descriptive output. 
C  10 Nov 99 @DynTel-BG 1-006  added index parameter constants, IDIURNAL, 
C            IHOTSOAK, ICRNCASE for the evaporative emissions array, 
C            GREVP. Added grams/start index, GPS, for start emissions. 
C   3 Nov 99 @DynTel-bg 1-006 Removed partial diurnals and liquid leak diurnals 
C   8 Oct 99 @DynTel-bg 2-000 Deleted output of '0' at the beginning of 
C            the weekend and user supplied temperature messages 
C   6 Sep 99 @DynTel-HXQ 2-684 The label is changed when user wants 
C             to express start emission as g/start 
C  13 Jul 99 @DynTel-ddj 2-699 Changed implementation logic of FTP_FLAG 
C            back so that FTP_FLAG is equal to 1 when invoked and 0 
C            otherwise. 
C  22 Jun 99 @DynTel-MLA 2-694 Deleted obsolete LEV variables. 
C   1 Jun 99 @DynTel-ddj 2-700 Changed implementation logic of FTP_FLAG 
C            to be consistent with other flags. 
C  26 May 99 @Dyntel-HQ 2-698 delete code of including scene1.i 
C  17 May 99 @DynTel-ddj 2-700 Added Message indicating FTP results when 
C            FTP FLAG set. 
C  14 May 99 @DynTel-bg  2-701 Adjusted write statements for refueling and running loss 
C   6 May 99 @DynTel-MLA 2-674 Changed output formats for the min and max 
C            temperatures from F4.0 to F5.1. TOTEVAP includes running, resting, 
C            diurnal, hot soak, refueling, and crankcase emissions. 
C   3 Mar 99 @DynTel-MLA 2-659 Restored call to OUTDTG. Merged code from 
C            OUTLAP into this routine (WRITE statement for TEMMIN, TEMMAX). 
C  16 Mar 99 @DynTel-bag 2-664 Removed TEMEVP and added output of whether temperature 
C            is from user input of 24 temperatures for hourly intervals. Renamed multiple diurnals 
C            to diurnals from liquid leaks. For runs which call for use of weekend vehicle 
C            activity fractions a note of this will be output (see format statement 217). 
C  03 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  06 Nov 98 @Dyntel-BAG 2-000 Renamed VMTGT1 and VMTGT2 to VMTGT12 and VMTGT34, respectively 
C  14-Sep-98 @Dyntel-bag 2-000 Added call to OUTZTT. 
C 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    modules : 
C     LEAKYDAT - ILIQUID_LEAKS 
C    common blocks: 
C    /ALUHIN/     AC,ABSHUM 
C    /CITPAR/     SCNAME 
C    /CITRV1/     RVPICY,RVPBAS 
C    /CVNAMS/     ILDGT12,ILDGT34,IHDGV,ILDDV,ILDDT,IHDDV,ILDGT (IVTYPE.I) 
C    /EVAL/       MONTH, MEVAL 
C    /EVAPGR/     GREVP,IGREVP 
C    /FLAGS1/     OXYFLG,DTGFLG 
C    /FLAGS2/     IMFLAG,ALHFLG,EVIMFLAG 
C    /FLAGS3/     ATPFLG 
C    /FLAGS4/     NMHFLG,HCFLAG,COLDFG 
C    /FLAGS5/     FTP_FLAG 
C    /IOUCOM/     IOUREP 
C    /OPCNTL/     OPTABL,EXPEXH 
C    /OXY1/       SHRMKT,OXYCNT,IGASHW 
C    /OXY2/       NFUEL 
C    /MAXIMA/     MAXPOL 
C    /P94INF/     P94AR 
C    /REGION/     IREJN,ALT,INITPR 
C    /RESUL1/     EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL 
C                 IEFFTP,IEFEXH,IEFEVAP,IEFREFL,IEFRUNL 
C    /RESUL2/     RLGGAL,RSTGPH,EFIDLE,G_MI_EVP, 
C                 IRLGGAL,IRSTGPH,IEFIDLE,IGMIEVP 
C    /RESUL3/     VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE 
C    /RESUL4/     IRUNEF,ISTARTEF,RUNEF,STARTEF 
C    /RFORM1/     RFGFLG 
C    /STRING/     NOYES,COMMA,PERIOD 
C    /SULFUR/     S_PPM 
C    /TEMPC5/     NCC,MYCOL,CCSTD 
C    /TEMPC6/     PIVPCT 
C    /TEMPS/      TEMEXH,TEMP_DIST 
C    /UNIT_PARAM/ GPS 
C    /VMXCOM/     VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTHDD 
C    /YEARS4/     IY1994 
C 
C  Local array subscripts: 
C 
C  NAMNMH(5)  -  NAMNMH ( NMHFLG ) 
C  NAMPOL(3)  -  NAMPOL ( IP ) 
C  NAMREG(3)  -  NAMREG ( IREJN ) 
C  NAMVEH(10) -  NAMVEH ( IVF3 ) 
C  NAMUNT(7,2) - NAMUNT ( ICH, IUNT ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IDLNAM   I    pointer for idle names 
C  IPIV     I    Phase-In model year pointer MYCOL=1, MYCOL+1=2, etc. 
C  NAMNMH   C*8  prefix indicating type of Exh.-HC (total, non-meth., VOC, TOG) 
C  NAMPOL   C*4  pollutant names 
C  NAMREG   C*4  region names 
C  NAMVEH   C*8  vehicle class names 
C  NAMUNT   C*12  component units, HCFLAG=3,5: mixed, HCFLAG=4: g/mi 
C  DASH8    C*8  8 character string of hyphens 
C  FMT1     I    format statement number assignment variable 
C  EMISS    C*5  holds name of pollutant type for output 
C  BUFFER1  C*80 holds character string for output 
C  BUFFER2  C*80 reformatted BUFFER1 from LJCOMP 
C  HEADER   C*80 message character string 
C  IVFIRST   I   print formatting index variable; sequences gas bus/schoolbus 
C  IVLAST    I   print formatting index variable; sequences gas bus/schoolbus 
C  BLANKLINE I   format statement number assignment variable 
C  IATP      I   index for NOYES(IATP) for ATP program 
C  IMF       I   index for NOYES(IMF) for Exh IM program 
C  EVIMF     I   index for NOYES(EVIMF) for Evap IM program 
C  IFUEL     I   index for writing out oxygentaed fuel percentages,SHRMKT 
C  ITER      I   write statement index 
C  IRFG      I   index for NOYES(IRFG) refueling loss 
C  IDIURNAL  I    evaporative emissions index (2=IDIURNAL) -  diurnal 
C  IHOTSOAK  I    evaporative emissions index (1=IHOTSOAK) -  hotsoak 
C  IRUNLOSS  I    evaporative emissions index (3=IRUNLOSS) -  running loss 
C  VMTSUM       R    summation variable 
C 
C  Notes: 
C 
C  OUTDT3 (of OUTFMT 3) prints the All Vehicle Classes average EFs to 3 dp. 
C 
C  OUTDT3 was modified for MOBILE4 In-House Version 02 to update the array name 
C  EFLOSS to REFL_GM_MI. 
C 
C  OUTDT3 was modified for MOBILE4,1 to print detailed HC evap components. 
C  in mixed or g/mi units. 
C 
C  22-April-1994 @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's 
C  September-30-1994 @CSC-pme request 2-464. Include message that indicate  
C  that the effects of detergent gasoline additive rulemaking have been  
C  altered by the user. (One LAP record for each scenario) 
C  Subroutine OUTDTG was created to print the  message. 
C  A call to subroutine OUTDTG was added.  
C 
C October-5-1994 @CSC pme request 2-472. A check for APPLEV variable to  
C see whether or not to apply LEV stds to LDGT2 was included. In case the 
C standards are not applied, the ZEV FRACTIONS are not printed. 
C 
C  08-Feb-1996 @ DynTel Request 2-615 Detergent gasoline effects are disabled 
C                when NEWFLG is set to 5 or 6. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE and altered this routine to use 
C            the new database index numbers for Mobile6. They have changed 
C            since Mobile5 because new vehicle types were introduced. 
C  12 November 1997 Dyntel@MMS 2-645 removed index of TEMEXH(IP) 
C            Moved up the line that prints out "Speed" before printing out the  
C            Ambient Temp. Previous version printed out multiple speed values.  
C            This version will print out only one value. 
C  9  January 1998 Dyntel@MMS 2-645 changed the format of the output 
C  30 Jan 98 @DynTel-MLA 2-622  Added the emission factors for BUS to the 
C            Heavy Duty Gasoline and Heavy Duty Diesel categories. 
C  10 Mar 98 @DynTel-RJD 2-622  Added calls to EFCVT so that emission 
C            factors for composite vehicles are not calculated in Write 
C            statements. Moved calculation of VMTMIX for composite vehicles 
C            to TFCALX. 
C  20 Apr 98 @DynTel-BAG 2-661 By model year output is only done in the  
C            spreadsheet output. So, this program no longer calls OUTHDR. 
C            Code from OUTPUT which writes header information (OUTHD3)  
C            and month and year information to IOUREP has been moved to  
C            the beginning of this code. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C 24-Aug-98 @Dyntel-nh Request 0-632 Removed dead code related to special I/M 
C 
      USE DATABASE, ONLY: DBSELPOL 
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS 
      IMPLICIT NONE 
      INCLUDE 'AC.I' 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'BASE10.I' 
      INCLUDE 'CITPAR.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'CO2COM.I' 
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'OXY2.I' 
      INCLUDE 'P94INF.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RESUL2.I' 
      INCLUDE 'RESUL3.I' 
      INCLUDE 'RESUL4.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'STRING.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'TEMPC5.I' 
      INCLUDE 'TEMPC6.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'UNIT_PARAM.I' 
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
      REAL MPGCVT 
C 
C     Local variables:                          
C 
      INTEGER       BLANKLINE, FMT0, FMT1, FMT2, FMT3 
      INTEGER       IATP, IFUEL, IMF, IP, ITER, EVIMF 
      INTEGER       IRFG 
      INTEGER, SAVE :: IDIURNAL=2,IHOTSOAK=1 
      INTEGER, SAVE :: ICRNCASE=4,IRESTLOSS=5 
      INTEGER, SAVE :: EXPND_HC=2 
      INTEGER, SAVE :: INO=0 
      REAL     VMTSUM 
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
      CHARACTER(LEN=5), DIMENSION(5) :: NAMNMH=(/ 
     *                   ' THC ',' NMHC',' VOC ',' TOG ',' NMOG'/) 
C 
      CHARACTER(LEN=5),DIMENSION(4)::NAMPOL=(/' HC  ',' CO  ',' NOX ', 
     *                                        ' CO2 '/) 
      CHARACTER(LEN=4), DIMENSION(3) :: NAMREG=(/'Low ','High','Mid '/) 
      CHARACTER(LEN=9) :: DASH8='  ------' 
      CHARACTER(LEN=5) :: EMISS 
C 
C 
      IF (ILIQUID_LEAKS==INO) THEN 
        WRITE(IOUREP,200) 
  200   FORMAT(//,T10,'Evaporative emissions sums DO NOT include fuel',/ 
     &         T10,'losses from vehicles with liquid leaks.'///) 
      ENDIF 
C 
  400 FORMAT() 
      ASSIGN 400 TO BLANKLINE 
      WRITE(IOUREP, BLANKLINE) 
C 
C     Print the name of file with LEV phase-in data. 
C 
      IF (P94AR(1).NE.' ') THEN 
        WRITE(IOUREP,204) P94AR(1)  
  204   FORMAT(' LEV phase-in data read from file ',A80) 
      ENDIF 
C 
      IF (IMFLAG.EQ.1) THEN 
        IMF = 1 
      ELSE 
        IMF = 2 
      ENDIF 
C 
      IF (EVIMFLAG.EQ.1) THEN 
        EVIMF = 1 
      ELSE 
        EVIMF = 2 
      ENDIF 
C 
      IF (ATPFLG.EQ.1) THEN 
        IATP = 1 
      ELSE 
        IATP = 2 
      ENDIF 
C 
      IF (RFGFLG.EQ.1) THEN 
        IRFG = 1 
      ELSE 
        IRFG = 2 
      ENDIF 
C 
   10 WRITE(IOUREP,210) ICY,MONTH(MEVAL),NAMREG(IREJN),TEMMIN,TEMMAX, 
     *                  ABSHUM,RVPBAS,RVPICY,S_PPM,NOYES(IMF), 
     *                  NOYES(EVIMF),NOYES(IATP),NOYES(IRFG) 
  210 FORMAT(T15,'      Calendar Year: ',I5/ 
     *       T15,'              Month: ',1X,A4/ 
     *       T15,'           Altitude: ',1X,A4/ 
     *       T15,'Minimum Temperature: ',F5.1,' (F)'/ 
     *       T15,'Maximum Temperature: ',F5.1,' (F)'/ 
     *       T15,'  Absolute Humidity: ',F5.0,' grains/lb'/ 
     *       T15,'   Nominal Fuel RVP: ',F5.1,' psi'/ 
     *       T15,'      Weathered RVP: ',F5.1,' psi'/ 
     *       T15,'Fuel Sulfur Content: ',F5.0,' ppm'// 
     *       T15,'Exhaust I/M Program: ',1X,A4/ 
     *       T15,'   Evap I/M Program: ',1X,A4/ 
     *       T15,'        ATP Program: ',1X,A4/ 
     *       T15,'   Reformulated Gas: ',1X,A4/) 
C 
      IF( TEMP_DIST ) THEN 
        WRITE (IOUREP,216) 
      END IF 
  216 FORMAT(T15,'User supplied hourly temperatures.'/) 
      IF( IWEEK .EQ. 2 ) THEN 
        WRITE (IOUREP,217) 
      END IF 
  217 FORMAT(T1,'Emissions determined from WEEKEND hourly ', 
     *'vehicle activity fractions.'/) 
C      
      IF(DTGFLG.EQ.2 .AND. NEWFLG.LT.5) CALL OUTDTG 
C 
C  Oxy Fuels 
C 
   15 IF(OXYFLG.EQ.2) WRITE(IOUREP,228) (SHRMKT(IFUEL),IFUEL=2,NFUEL), 
     *                                   OXYCNT,NOYES(IGASHW) 
  228 FORMAT(T4,'Ether Blend Market Share:',F6.3, 
     *       7X,'Alcohol Blend Market Share:',F6.3/ 
     *       T4,'Ether Blend Oxygen Content:',F6.3, 
     *       5X,'Alcohol Blend Oxygen Content:',F6.3/ 
     *       T43,'Alcohol Blend RVP Waiver: ',A4) 
C 
      VMTSUM=VMTMIX(1)+VMTGT12+VMTGT34+VMTHDG 
     *   +VMTMIX(14)+VMTLDDT+VMTHDD+VMTMIX(24) 
      WRITE(IOUREP,230) 
     *   NAMVEH,(DASH8,ITER=1,10), 
     *   VMTMIX(1),VMTGT12,VMTGT34,VMTHDG, 
     *   VMTMIX(14),VMTLDDT,VMTHDD,VMTMIX(24),VMTSUM 
  230 FORMAT(T8,'Vehicle Type:',10(2X,A8)/ 
     *       T16,'GVWR:',T36,'<6000',T46,'>6000',T56,'(All)'/ 
     *       T21,10(2X,A8)/ 
     *       T4,'VMT Distribution:',3F10.4,10X,6F10.4) 
C 
C  If CO2 Is selected, print the harmonic mean MPG values. 
C 
      IF(DBSELPOL(4).EQ.2) WRITE(IOUREP,235) MPG(1),MPGCVT(ILDGT12), 
     *    MPGCVT(ILDGT34),MPGCVT(ILDGT),MPGCVT(IHDGV),MPG(14), 
     *    MPGCVT(ILDDT),MPGCVT(IHDDV),MPG(24),FMPG 
  235 FORMAT(T2,'Fuel Economy (mpg):',10F10.1) 
C 
C     Print the composite emission factor group. 
C 
      WRITE(IOUREP,240) 
  240 FORMAT(2X,118('-'),/T2, 
     *   'Composite Emission Factors (g/mi):') 
C 
      DO IP=1,MAXPOL+1   ! +1 to include CO2. 
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
        WRITE(IOUREP,FMT0) EMISS,EFFTP(IP,1), 
     *    EFCVT(ILDGT12,IEFFTP,IP),EFCVT(ILDGT34,IEFFTP,IP), 
     *    EFCVT(ILDGT,IEFFTP,IP),EFCVT(IHDGV,IEFFTP,IP),EFFTP(IP,14), 
     *    EFCVT(ILDDT,IEFFTP,IP),EFCVT(IHDDV,IEFFTP,IP),EFFTP(IP,24), 
     *    VFTP(IP) 
  250   FORMAT(T6,'Composite',A5,':',9F10.2,F10.3) 
  251   FORMAT(T6,'Composite',A5,': ',5F10.3,F9.3,2F10.3,F10.2,F10.3) 
  252   FORMAT(T6,'Composite',A5,':',5F10.2,3F10.3,F10.2,F10.3) 
  253   FORMAT(T6,'Composite',A5,':',9(F9.1,1X),F9.2) 
  254   FORMAT(T6,'Composite',A5,':',5(F9.1,1X),F8.1,1X,2(F9.1,1X), 
     *                                F10.1,F10.2) 
C 
      END DO 
C 
CC      WRITE(BUFFER1,252) NAMNMH(NMHFLG) 
CC  252 FORMAT('(Composite',A5,' emissions do not include', 
CC     *       ' refueling emissions.)') 
CC      CALL LJCOMP(BUFFER1,BUFFER2) 
CC      WRITE(IOUREP,253) BUFFER2 
CC  253 FORMAT(T2,A80,/2X,118('-')) 
      WRITE(IOUREP,255) 
  255 FORMAT(2X,118('-')) 
C 
C    Print Expanded Exhaust, if requested 
C 
      IF(EXPEXH) THEN 
C       
        WRITE(IOUREP,260) 
  260   FORMAT('Exhaust emissions (g/mi):') 
        IF (NMHFLG.EQ.2.OR.NMHFLG.EQ.5) THEN 
          ASSIGN 262 TO FMT1 
          IF(ICY>1999)  ASSIGN 264 TO FMT1 
        ELSE 
          ASSIGN 263 TO FMT1 
          IF(ICY>1999)  ASSIGN 265 TO FMT1 
        ENDIF 
  262   FORMAT(T9,A5,T15,'Start:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *       T7,A5,T13,'Running:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *       T1,A5,T7 ,'Total Exhaust:',9F10.2,F10.3) 
  263   FORMAT(T8,A5,T15,'Start:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *       T6,A5,T13,'Running:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *       T2,A5,T7 ,'Total Exhaust:',9F10.2,F10.3) 
  264   FORMAT(T9,A5,T15,'Start:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *       T7,A5,T13,'Running:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *       T1,A5,T7 ,'Total Exhaust:',8F10.3,F10.2,F10.3) 
  265   FORMAT(T8,A5,T15,'Start:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *       T6,A5,T13,'Running:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *       T2,A5,T7 ,'Total Exhaust:',8F10.3,F10.2,F10.3) 
        IF(DBSELPOL(1).EQ.2) WRITE(IOUREP,FMT1) 
     *    NAMNMH(NMHFLG),STARTEF(1,1),EFCVT(ILDGT12,ISTARTEF,1), 
     *    EFCVT(ILDGT34,ISTARTEF,1),STARTEF(1,ILDGT), 
     *    STARTEF(1,14),EFCVT(ILDDT,ISTARTEF,1),STARTEF(1,24), 
     *    NAMNMH(NMHFLG),RUNEF(1,1),EFCVT(ILDGT12,IRUNEF,1), 
     *    EFCVT(ILDGT34,IRUNEF,1),RUNEF(1,ILDGT), 
     *    RUNEF(1,14),EFCVT(ILDDT,IRUNEF,1),RUNEF(1,24), 
     *    NAMNMH(NMHFLG), EFEXH(1),EFCVT(ILDGT12,IEFEXH,1), 
     *    EFCVT(ILDGT34,IEFEXH,1),EFEXH(ILDGT), 
     *    EFCVT(IHDGV,IEFEXH,1),EFEXH(14), 
     *    EFCVT(ILDDT,IEFEXH,1),EFCVT(IHDDV,IEFEXH,1),EFEXH(24),VEXH 
C 
        ASSIGN 270 TO FMT2 
        IF(ICY>1999)  ASSIGN 271 TO FMT2 
        IF(DBSELPOL(2).EQ.2) WRITE(IOUREP,FMT2) 
     *    STARTEF(2,1),EFCVT(ILDGT12,ISTARTEF,2), 
     *    EFCVT(ILDGT34,ISTARTEF,2),STARTEF(2,ILDGT), 
     *    STARTEF(2,14),EFCVT(ILDDT,ISTARTEF,2), 
     *    STARTEF(2,24), 
     *    RUNEF(2,1),EFCVT(ILDGT12,IRUNEF,2), 
     *    EFCVT(ILDGT34,IRUNEF,2),RUNEF(2,ILDGT), 
     *    RUNEF(2,14),EFCVT(ILDDT,IRUNEF,2), 
     *    RUNEF(2,24), 
     *    EFFTP(2,1),EFCVT(ILDGT12,IEFFTP,2),EFCVT(ILDGT34,IEFFTP,2), 
     *    EFFTP(2,ILDGT),EFCVT(IHDGV,IEFFTP,2),EFFTP(2,14), 
     *    EFCVT(ILDDT,IEFFTP,2),EFCVT(IHDDV,IEFFTP,2),EFFTP(2,24), 
     *    VFTP(2) 
C 
  270     FORMAT(/T12,'CO Start:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *        T10,'CO Running:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *        T4,'CO Total Exhaust:',9F10.2,F10.3) 
  271     FORMAT(/T12,'CO Start:',4F10.2,10X,2F10.3,11X,F10.3/ 
     *        T10,'CO Running:',4F10.2,10X,2F10.3,11X,F10.3/ 
     *        T4,'CO Total Exhaust:',5F10.2,3F10.3,F10.2,F10.3) 
C 
        ASSIGN 272 TO FMT3 
        IF(ICY>1999)  ASSIGN 273 TO FMT3 
        IF(DBSELPOL(3).EQ.2) WRITE(IOUREP,FMT3) 
     *    STARTEF(3,1),EFCVT(ILDGT12,ISTARTEF,3), 
     *    EFCVT(ILDGT34,ISTARTEF,3),STARTEF(3,ILDGT), 
     *    STARTEF(3,14),EFCVT(ILDDT,ISTARTEF,3), 
     *    STARTEF(3,24), 
     *    RUNEF(3,1),EFCVT(ILDGT12,IRUNEF,3), 
     *    EFCVT(ILDGT34,IRUNEF,3),RUNEF(3,ILDGT), 
     *    RUNEF(3,14),EFCVT(ILDDT,IRUNEF,3), 
     *    RUNEF(3,24), 
     *    EFFTP(3,1),EFCVT(ILDGT12,IEFFTP,3),EFCVT(ILDGT34,IEFFTP,3), 
     *    EFFTP(3,ILDGT),EFCVT(IHDGV,IEFFTP,3),EFFTP(3,14), 
     *    EFCVT(ILDDT,IEFFTP,3),EFCVT(IHDDV,IEFFTP,3),EFFTP(3,24), 
     *    VFTP(3) 
C 
  272   FORMAT(/T11,'NOx Start:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *        T9,'NOx Running:',4F10.2,10X,2F10.2,10X,F10.2/ 
     *        T3,'NOx Total Exhaust:',9F10.2,F10.3) 
  273   FORMAT(/T11,'NOx Start:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *        T9,'NOx Running:',4F10.3,10X,2F10.3,11X,F10.3/ 
     *        T3,'NOx Total Exhaust:',8F10.3,F10.2,F10.3) 
C 
        WRITE(IOUREP,282) 
  282   FORMAT(2X,118('-')) 
C 
C  END OF EXPANDED EXHAUST SECTION 
      ENDIF 
C 
C     Expanded evaporative/HC output, HCFLAG > 1 
   50 IF(HCFLAG.GE.EXPND_HC .AND. DBSELPOL(1).EQ.2) 
     *THEN 
C 
          WRITE(IOUREP,285) 
  285     FORMAT('Non-Exhaust Emissions (g/mi):') 
C 
          WRITE(IOUREP,290) 
     *     G_MI_EVP(IHOTSOAK,1),EFCVT(ILDGT12,IGMIEVP,IHOTSOAK), 
     *     EFCVT(ILDGT34,IGMIEVP,IHOTSOAK),G_MI_EVP(IHOTSOAK,ILDGT), 
     *     EFCVT(IHDGV,IGMIEVP,IHOTSOAK),G_MI_EVP(IHOTSOAK,24) 
  290     FORMAT(T7,'Hot Soak Loss:',5F10.3,30X,F10.3) 
C 
          WRITE(IOUREP,295) 
     *     G_MI_EVP(IDIURNAL,1),EFCVT(ILDGT12,IGMIEVP,IDIURNAL), 
     *     EFCVT(ILDGT34,IGMIEVP,IDIURNAL),G_MI_EVP(IDIURNAL,ILDGT), 
     *     EFCVT(IHDGV,IGMIEVP,IDIURNAL),G_MI_EVP(IDIURNAL,24) 
  295     FORMAT(T8,'Diurnal Loss:',5F10.3,30X,F10.3) 
C 
          WRITE(IOUREP,300) 
     *     G_MI_EVP(IRESTLOSS,1),EFCVT(ILDGT12,IGMIEVP,IRESTLOSS), 
     *     EFCVT(ILDGT34,IGMIEVP,IRESTLOSS),G_MI_EVP(IRESTLOSS,ILDGT), 
     *     EFCVT(IHDGV,IGMIEVP,IRESTLOSS),G_MI_EVP(IRESTLOSS,24) 
  300     FORMAT(T8,'Resting Loss:',5F10.3,30X,F10.3,F10.3) 
C 
          WRITE(IOUREP,301) 
     *     EFRUNL(1),EFCVT(ILDGT12,IEFRUNL,IP), 
     *     EFCVT(ILDGT34,IEFRUNL,IP),EFRUNL(ILDGT), 
     *     EFCVT(IHDGV,IEFRUNL,IP) 
  301     FORMAT(T8,'Running Loss:',5F10.3,40X,F10.3) 
C 
          WRITE(IOUREP,305) 
     *     G_MI_EVP(ICRNCASE,1),EFCVT(ILDGT12,IGMIEVP,ICRNCASE), 
     *     EFCVT(ILDGT34,IGMIEVP,ICRNCASE),G_MI_EVP(ICRNCASE,ILDGT), 
     *     EFCVT(IHDGV,IGMIEVP,ICRNCASE), G_MI_EVP(ICRNCASE,24) 
  305     FORMAT(T6,'Crankcase Loss:',5F10.3,30X,F10.3) 
C 
          WRITE(IOUREP,310) 
     *     REFL_GM_MI(1),EFCVT(ILDGT12,IEFREFL,IP), 
     *     EFCVT(ILDGT34,IEFREFL,IP),REFL_GM_MI(ILDGT), 
     *     EFCVT(IHDGV,IEFREFL,IP) 
  310     FORMAT(T6,'Refueling Loss:',5F10.3,40X,F10.3) 
C 
        WRITE(IOUREP,320) TOTEVAP(1),EFCVT(ILDGT12,ITEVAP,1), 
     *    EFCVT(ILDGT34,ITEVAP,IP),TOTEVAP(ILDGT), 
     *    EFCVT(IHDGV,ITEVAP,IP),TOTEVAP(ILDDV), 
     *    EFCVT(ILDDT,ITEVAP,IP),EFCVT(IHDDV,ITEVAP,IP), 
     *    TOTEVAP(24),VTEVAP 
  320   FORMAT(T3,'Total Non-Exhaust:',9F10.3,F10.3/2X,118('-')) 
C 
CCC        IF(HCFLAG.EQ.ALT_EVAP_UNITS) THEN 
CCC          WRITE(IOUREP,325) 
CCC  325     FORMAT('Non-Exhaust Emissions (alternate units):') 
C 
CCC          WRITE(IOUREP,330) 
CCC     *     GREVP(IHOTSOAK,1),EFCVT(ILDGT12,IGREVP,1), 
CCC     *     EFCVT(ILDGT34,IGREVP,1),GREVP(IHOTSOAK,ILDGT), 
CCC     *     EFCVT(IHDGV,IGREVP,1),GREVP(IHOTSOAK,24) 
CCC  330     FORMAT(T5,'Hot Soak (g/trp):',5F10.2,30X,F10.2) 
C 
CCC          WRITE(IOUREP,335) 
CCC     *     GREVP(IDIURNAL,1),EFCVT(ILDGT12,IGREVP,2), 
CCC     *     EFCVT(ILDGT34,IGREVP,2),GREVP(IDIURNAL,ILDGT), 
CCC     *     EFCVT(IHDGV,IGREVP,2), GREVP(IDIURNAL,24) 
CCC  335     FORMAT(T7,'Diurnal (g/hr):',5F10.2,30X,F10.2) 
C 
CCC          WRITE(IOUREP,340) 
CCC     *     RSTGPH(1),EFCVT(ILDGT12,IRSTGPH,IP), 
CCC     *     EFCVT(ILDGT34,IRSTGPH,IP),RSTGPH(ILDGT), 
CCC     *     EFCVT(IHDGV,IRSTGPH,IP),RSTGPH(24) 
CCC  340     FORMAT(T7,'Resting (g/hr):',5F10.2,30X,F10.2) 
C 
CCC          WRITE(IOUREP,345) 
CCC     *     RLGGAL(1),EFCVT(ILDGT12,IRLGGAL,IP), 
CCC     *     EFCVT(ILDGT34,IRLGGAL,IP),RLGGAL(ILDGT), 
CCC     *     EFCVT(IHDGV,IRLGGAL,IP) 
CCC  345     FORMAT(T4,'Refueling (g/gal):',5F10.2/2X,118('-')) 
C 
CCC        ENDIF 
C 
      ENDIF 
C 
C  Check the expanded output requests. 
C 
   85 IF (OPTABL(1).EQ.2) CALL OUTLDT(ICY) 
      IF (OPTABL(2).EQ.2) CALL OUTHDV(1,ICY) 
      IF (OPTABL(3).EQ.2) CALL OUTHDV(2,ICY) 
      IF (OPTABL(4).EQ.2) CALL OUTBUS(ICY) 
C 
      RETURN 
      END 
