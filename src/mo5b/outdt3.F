C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDT3(ICY)
C
C  OUTDT3 writes scenario parameters and output on unit IOUREP using output
C  format (OUTFMT) 3.
C
C  Called by OUTPUT.
C
C  Calls OUTHDR, OUTLAP, OUTDTG and OUTPOL.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /ALUHIN/ AC,XLOAD,TRAILR,ABSHUM,DB,WB
C    /CITCIN/ PFRATE(3,2),IPFCU
C    /CITRV1/ RVPICY
C    /EVAPGR/ GREVP,VGREVP,GRPD
C    /FLAGS1/ OXYFLG,DTGFLG
C    /FLAGS2/ IMFLAG,ALHFLG
C    /FLAGS3/ ATPFLG
C    /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVFLG,APPLEV
C    /MAXIMA/ MAXVEH
C    /OXY1/   SHRMKT,OXYCNT,IGASHW
C    /OXY2/   NFUEL
C    /REGION/ IREJN,ALT
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,EFREFL,EFRUNL
C    /RESUL2/ RLGGAL,RSTGPH,EFIDLE,EFRSTL
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE
C    /RFORM1/ RFGFLG
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /STRING/ NOYES,COMMA,PERIOD
C    /TEMPC5/ NCC,MYCOL,CCSTD
C    /TEMPC6/ PIVPCT
C    /TEMPS/  TEMEVP,TEMEXH
C    /USDATA/ USNAME,NUSD,IUSD
C    /VMXCOM/ VMTMIX
C
C  Local array subscripts:
C
C  NAMEVP(5)  -  NAMEVP ( IEVP )
C  NAMIDL(5)  -  NAMIDL ( IDLNAM )
C  NAMNMH(5)  -  NAMNMH ( NMHFLG )
C  NAMPAR(3)  -  NAMPAR ( IPARTL )
C  NAMPOL(3)  -  NAMPOL ( IP )
C  NAMREG(3)  -  NAMREG ( IREJN )
C  NAMVEH(10) -  NAMVEH ( IVF3 )
C  NAMUNT(7,2) -  NAMUNT ( ICH, IUNT )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IDLNAM   I    pointer for idle names
C  IP1      I    see parameter dictionary
C  IP2      I    see parameter dictionary
C  IPIV     I    Phase-In model year pointer MYCOL=1, MYCOL+1=2, etc.
C  NAMIDL   C*3  prefix indicating type of Idle (total, non-methane, VOC, TOG)
C  NAMNMH   C*8  prefix indicating type of Exh.-HC (total, non-meth., VOC, TOG)
C  NAMPAR   C*9  partial diurnal cases
C  NAMPOL   C*4  pollutant names
C  NAMREG   C*4  region names
C  NAMVEH   C*8  vehicle class names
C  NAMUNT   C*12  component units, HCFLAG=3,5: mixed, HCFLAG=4: g/mi
C  NWITHC   I    number of substrings to be written with a comma appended
C  USC8     C*8  8 character underscore
C
C  Notes:
C
C  OUTDT3 (of OUTFMT 3) prints the All Vehicle Classes average EFs to 3 dp.
C
C  OUTDT3 was modified for MOBILE4 In-House Version 02 to update the array name
C  EFLOSS to EFREFL.
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
      INCLUDE 'ALUHIN.I'
      INCLUDE 'BASE10.I'
      INCLUDE 'CITCIN.I'
      INCLUDE 'CITRV1.I'
      INCLUDE 'EVAPGR.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'OXY2.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'STRING.I'
      INCLUDE 'TEMPC5.I'
      INCLUDE 'TEMPC6.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'USDATA.I'
      INCLUDE 'VMXCOM.I'
C
      CHARACTER*9 NAMVEH(10),NAMNMH(5),USC8,NAMEVP(5),NAMPAR(3)
      CHARACTER*12 NAMUNT(7,2)
C
      CHARACTER*3 NAMIDL(5)
      CHARACTER*4 NAMREG(3),NAMPOL(3)
C
      DATA NAMVEH/
     *   '    LDGV',
     *   '   LDGT1',
     *   '   LDGT2',
     *   '    LDGT',
     *   '    HDGV',
     *   '    LDDV',
     *   '    LDDT',
     *   '    HDDV',
     *   '    MC  ',
     *   ' All Veh'/
C
      DATA NAMREG/'Low ','High','Mid '/
C
      DATA USC8/'  ______'/
C
      DATA NAMPOL/' HC:',' CO:','NOX:'/
C
      DATA NAMNMH/'Total    ','Non-Meth ','VOC      ',
     *            'TOG      ','NMOG     '/
C
      DATA NAMIDL/'   ','NM-','VOC','TOG','NMO'/
C
      DATA NAMEVP/'Hot Soak ','Diurnal  ','Multiple ','Crankcase',
     *            'WtDiurnal'/
C
      DATA NAMPAR/'8am-11am ','10am-3pm ','8am-2pm  '/
C
      DATA NAMUNT/
     3      '(Hot Soak: g','/trip, Diurn','als: g, Cran','kcase: g/mi,',
     3      ' Refuel: g/g','al, Resting:',' g/hr)      ',
     4      '(All Compone','nts in Grams',' per Mile)  ','            ',
     4      '            ','            ','            '/
C
C  If BYMY run, start each scenario after 1st on new page.
C
      IF(OUTFMT.EQ.5) CALL OUTHDR(IOUREP)
C
      CALL OUTPOL(IP1,IP2)
C
      IMF=1
      IF(IMFLAG.NE.1) IMF=2
      IF(NUSD.EQ.0) GOTO 10
      IF(NUSD.EQ.1) WRITE(IOUREP,200)(USNAME(ICH,IUSD(1)),ICH=1,4),
     *                               PERIOD
      NWITHC=NUSD-1
      IF(NWITHC.GT.0) WRITE(IOUREP,200)
     *  ((USNAME(ICH,IUSD(ICT)),ICH=1,4),COMMA,ICT=1,NWITHC),
     *  (USNAME(ICH,IUSD(NUSD)),ICH=1,4),PERIOD
  200 FORMAT('0','User supplied',4A8,T46,A1,4A8,T78,A1,4A8,T110,A1/
     *       ' ',T15,4A8,T46,A1)
C
   10 IF(PRTFLG.LT.4) WRITE(IOUREP,210)
     *   ICY,NOYES(IMF),TEMEXH(IP1),NAMREG(IREJN)
  210 FORMAT(
     *'0', T2,'Cal. Year:',I5,
     *    T25,'I/M Program: ',A4,
     *    T48,'Ambient Temp:',F5.1,' (F)',
     *    T85,'Region: ',A4)
C
      IF(PRTFLG.EQ.4) WRITE(IOUREP,215)
     *   ICY,NOYES(IMF),TEMEXH,NAMREG(IREJN)
  215 FORMAT(
     *'0', T2,'Cal. Year:',I5,
     *    T25,'I/M Program: ',A4,
     *    T48,'Ambient Temp:',F5.1,2(' /',F5.1),' (F)',
     *    T85,'Region: ',A4)
C
      IATP = ATPFLG
      IF(ATPFLG.NE.1) IATP = 2
      WRITE(IOUREP,220)
     *       NOYES(IATP),PCCN,PCHC,PCCC,ALT
  220 FORMAT(
     *' ',T19,'Anti-tam. Program: ',A4,
     *    T46,'Operating Mode:',F5.1,2(' /',F5.1),
     *    T83,'Altitude:',F6.0,' Ft.')
C
      IF(RFGFLG.EQ.2) THEN
        WRITE(IOUREP,222) NOYES(RFGFLG),ASTMCL
      ELSE
        WRITE(IOUREP,223) NOYES(RFGFLG)
      ENDIF
  222 FORMAT(
     *' ',T19,' Reformulated Gas: ',A4,T49,' ASTM Class: ',A1)
  223 FORMAT(
     *' ',T19,' Reformulated Gas: ',A4)
C
      IF(IPFCU.EQ.2) THEN
      DO 12 KEYEQ=1,3
      PFRATE(KEYEQ,IPFCU)=PFRATE(KEYEQ,IPFCU)*100.
   12 CONTINUE
      WRITE(IOUREP,224) (PFRATE(KEYEQ,IPFCU),KEYEQ=1,3)
  224 FORMAT(
     *' ',T19,'User-Supplied Pass, Purge Fail, ',
     *'and Pressure Fail Rates (%):',2(F5.0,','),F5.0)
      ENDIF
C
      IF(LOCFLG.EQ.1) THEN
         CALL OUTLAP
         IF(DTGFLG.EQ.2 .AND. NEWFLG.LT.5) CALL OUTDTG 
      ENDIF
C
      IF(ALHFLG.GT.1) WRITE(IOUREP,225)
     *   ABSHUM,
     *   AC,DB,WB
  225 FORMAT(
     *' ',T19,'Absolute Humidity:',F6.2,
     *    T48,'AC (DB / WB):',F5.1,' (',F5.1,' /',F5.1,')')
C
C  The COLD CO doc output is skipped in MOBILE4.1.
C
      IF(COLDFG.EQ.1) GOTO 15
      DO 13 IC=1,NCC
C     WRITE(IOUREP,227)
C    *    IC,MYCOL(IC),(CCSTD(IV,IC),IV=1,3),
C    *    ((PIVPCT(IPIV,IV,IC),IPIV=1,3),IV=1,3)
  227 FORMAT(
     *' ', T2,'Cold CO Standard #',I1,
     *    T23,'Start MY:',I5,
     *    T38,'Stds by class:',
     *    T53,'LDGV',F5.1,'  LDGT1',F5.1,'  LDGT2',F5.1/
     *    T23,'Phase-In %: LDGV:',3F6.1,' LDGT1:',3F6.1,' LDGT2:',3F6.1)
   13 CONTINUE
C
C  Oxy Fuels
C
   15 IF(OXYFLG.EQ.2) WRITE(IOUREP,228) (SHRMKT(IFUEL),IFUEL=2,NFUEL),
     *                                   OXYCNT,NOYES(IGASHW)
  228 FORMAT('0',2X,'Ether Blend Market Share:',F6.3,
     *           7X,'Alcohol Blend Market Share:',F6.3/
     *       ' ',   'Ether Blend Oxygen Content:',F6.3,
     *           5X,'Alcohol Blend Oxygen Content:',F6.3/
     *       ' ',42X,'Alcohol Blend RVP Waiver: ',A4)
C
      WRITE(IOUREP,230)
     *   NAMVEH,(USC8,ITER=1,10),
     *   SPD,
     *   VMTMIX
  230 FORMAT(
     *'0','  Veh. Type:',10(2X,A8)/'+',T15,9(2X,A8),1X,A8/
     *' ','Veh. Speeds:',F9.1,2F10.1,10X,5F10.1/
     *' ','    VMT Mix:',F11.3,2F10.3,10X,5F10.3)
C
      IF(ALHFLG.GT.1) WRITE(IOUREP,235)
     *   XLOAD,
     *   TRAILR
  235 FORMAT(
     *' ',' Extra Load:',F11.3,2F10.3/
     *' ','Trlr in Tow:',F11.3,2F10.3)
C
C  (22-April-1994) @ CSC ljn Request 2-420
C
      IF(LEVFLG.EQ.2) THEN
         IF(APPLEV.EQ.2)THEN
           WRITE(IOUREP,233) (ZEVTOT(I),I=1,3)
  233      FORMAT(
     *       ' ','  ZEV Fract:',F10.2,' %',F8.2,' %',F8.2,' %')
         ELSE
           WRITE(IOUREP,234) (ZEVTOT(I),I=1,2)
  234      FORMAT(
     *       ' ','  ZEV Fract:',F10.2,' %',F8.2,' %')
         ENDIF
      ENDIF
C
      WRITE(IOUREP,240)
  240 FORMAT(
     *'0','Composite Emission Factors (Gm/Mile)')
C
      DO 30 IP=IP1,IP2
      IF(IP.NE.1) GOTO 20
C
C  HC pass: print ef's scaled in grams.
C
      WRITE(IOUREP,250)
     *  NAMNMH(NMHFLG),NAMPOL(IP),
     *  (EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     *  VFTP(IP)
  250 FORMAT(
     *' ',A8,A4,
     *   9F10.2,
     *    F10.3)
C
C  HCFLAG >= 2 means write out HC components also.
C
C  MOBILE4 printed blanks for refueling HC, not zeroes.
C
      IF(HCFLAG.GE.2) WRITE(IOUREP,260)
     *   NAMPOL(IP),
     *   (EFEXH(IV),IV=1,3),EFEXH(9),(EFEXH(IV),IV=4,MAXVEH),VEXH,
     *   NAMPOL(IP),
     *   (EFEVAP(IV),IV=1,3),EFEVAP(9),EFEVAP(4),EFEVAP(8),VEVAP,
     *   NAMPOL(IP),
     *   (EFREFL(IV),IV=1,3),EFREFL(9),EFREFL(4),VLOSS,
     *   NAMPOL(IP),
     *   (EFRUNL(IV),IV=1,3),EFRUNL(9),EFRUNL(4),VRUNLS,
     *   NAMPOL(IP),
     *   (EFRSTL(IV),IV=1,3),EFRSTL(9),EFRSTL(4),EFRSTL(8),VRSTLS
  260 FORMAT(
     *' ','Exhaust ',A4,9F10.2,F10.3/
     *' ','Evaporat',A4,5F10.2,30X,F10.2,F10.3/
     *' ','Refuel L',A4,5F10.2,40X,F10.3/
     *' ','Runing L',A4,5F10.2,40X,F10.3/
     *' ','Rsting L',A4,5F10.2,30X,F10.2,F10.3)
      GOTO 30
C
   20 WRITE(IOUREP,265)
     *  NAMPOL(IP),
     *  (EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     *  VFTP(IP)
  265 FORMAT(
     *' ','Exhaust ',A4,
     *    9F10.2,
     *    F10.3)
   30 CONTINUE
C
      IF(IDLFLG.EQ.1) GOTO 50
      WRITE(IOUREP,270)
  270 FORMAT(
     *'0','Idle Emission Factors (Gm/Hr)')
C
      DO 40 IP=IP1,IP2
C
      IDLNAM=NMHFLG
      IF(IP.NE.1) IDLNAM=1
C
      WRITE(IOUREP,275)
     *  NAMIDL(IDLNAM),NAMPOL(IP),
     *  (EFIDLE(IP,IV),IV=1,3),EFIDLE(IP,9),(EFIDLE(IP,IV),IV=4,MAXVEH),
     *  VIDLE(IP)
  275 FORMAT(
     *' ',A3,'Idle ',A4,
     *    9F10.2,
     *    F10.3)
   40 CONTINUE
C
   50 IF(HCFLAG.LT.3.OR.PRTFLG.EQ.2.OR.PRTFLG.EQ.3) GOTO 80
      IUNT=HCFLAG-2
      IF(IUNT.GT.2) IUNT=1
      WRITE(IOUREP,285) RVPICY,TEMEVP(1),(NAMUNT(ICH,IUNT),ICH=1,7),
     *                 TEMEVP(2),TEMEVP(6)
  285 FORMAT(
     *'0','Evaporative Emissions by Component',
     *21X,'Weathered RVP:',F5.1,
     *T91,'Hot Soak Temp:',F5.1,' (F)'/
     *' ',7A12,
     *T87,'Running Loss Temp:',F5.1,' (F)'/
     *' ',
     *T87,'Resting Loss Temp:',F5.1,' (F)')
C
      DO 60 IEVP=1,4
C
      INEVP=IEVP
      IF((HCFLAG.EQ.3.OR.HCFLAG.EQ.5).AND.IEVP.EQ.2)INEVP=5
C
      IF(IEVP.NE.3) WRITE(IOUREP,290)
     *  NAMEVP(INEVP),(GREVP(IEVP,IV),IV=1,3),GREVP(IEVP,9),
     *  GREVP(IEVP,4),GREVP(IEVP,MAXVEH)
  290 FORMAT(' ',A9,3X,5F10.2,30X,F10.2)
C
      IF(IEVP.EQ.3) WRITE(IOUREP,295)
     *  NAMEVP(IEVP),(GREVP(IEVP,IV),IV=1,3),GREVP(IEVP,9),
     *  GREVP(IEVP,4)
  295 FORMAT(' ',A9,3X,5F10.2)
   60 CONTINUE
C
C   DEBUG
C     DO 41 I = 1, 4
C       WRITE(IOUREP,296) I, (GREVP(I,IV),IV=1,3)
C 296   FORMAT(
C    *' ','index =',I2,' Hot Soak 1 =',F10.2,
C    *' Hot Soak 2 =',F10.2,' Hot Soak 3 =',F10.2 )
C 41  CONTINUE
C
      IF(HCFLAG.EQ.4) THEN
      DO 70 IDU=2,4
      IPARTL=IDU-1
      WRITE(IOUREP,297) NAMPAR(IPARTL),(GRPD(IPARTL,IV),IV=1,3),
     *                 GRPD(IPARTL,9)
  297 FORMAT(' ',A9,3X,4F10.2)
   70 CONTINUE
      ENDIF
C
      IF (HCFLAG.EQ.3.OR.HCFLAG.EQ.5)
     *    WRITE(IOUREP,298) (RLGGAL(IV),IV=1,3),RLGGAL(9),RLGGAL(4),
     *    (RSTGPH(IV),IV=1,3),RSTGPH(9),RSTGPH(4),RSTGPH(8)
  298 FORMAT(' ','Refuel',6X,5F10.2/
     *       ' ','Resting',5X,5F10.2,30X,F10.2)
C
   80 WRITE(IOUREP,300)
  300 FORMAT(' ',112('_'))
C
      RETURN
      END
