C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDT4(ICY)
C
C  OUTDT4 writes scenario parameters and output on unit IOUREP using output
C  format (OUTFMT) 4.
C
C  Called by OUTPUT.
C
C  Calls OUTLAP and OUTPOL.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /ALUHIN/ AC,XLOAD,TRAILR,ABSHUM,DB,WB
C    /CITCIN/ PFRATE(3,2),IPFCU
C    /FLAGS1/ OXYFLG
C    /FLAGS2/ IMFLAG,ALHFLG
C    /FLAGS3/ ATPFLG
C    /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVFLG
C    /MAXIMA/ MAXVEH
C    /OXY1/ SHRMKT,OXYCNT,IGASHW
C    /OXY2/ NFUEL
C    /REGION/ IREJN,ALT
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,EFREFL,EFRUNL
C    /RESUL2/ EFIDLE,EFRSTL
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE
C    /RFORM1/ RFGFLG
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /TEMPS/  TEMEXH
C    /STRING/ NOYES,COMMA,PERIOD
C    /USDATA/ USNAME,NUSD,IUSD
C    /VMXCOM/ VMTMIX
C
C  Local array subscripts:
C
C  NAMIDL(5)  -  NAMIDL ( IDLNAM )
C  NAMNMH(5)  -  NAMNMH ( NMHFLG )
C  NAMPOL(3)  -  NAMNMH ( IP )
C  NAMREG(3)  -  NAMREG ( IREJN )
C  NAMVEH(10) -  NAMVEH ( IVF3 )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IDLNAM   I    pointer for idle names
C  IP1      I    see parameter dictionary
C  IP2      I    see parameter dictionary
C  NAMIDL   C*5  prefix indicating type of Idle (total, non-methane, VOC, TOG)
C  NAMNMH   C*8  prefix indicating type of HC (total, non-methane, VOC)
C  NAMPOL   C*4  pollutant names
C  NAMREG   C*4  region names
C  NAMVEH   C*8  vehicle class names
C  NWITHC   I    number of substrings to be written with a comma appended
C  USC8     C*8  8 character underscore
C
C  Notes:
C
C  OUTDT4 was modified for MOBILE4 In-House Version 02 to update the array name
C  EFLOSS to EFREFL.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
      INTEGER RFGFLG,WINFLG
      CHARACTER*8 USNAME
      CHARACTER*4 NOYES,COMMA,PERIOD
      CHARACTER*1 COLON,ASTMCL
      LOGICAL RFGON
C
C
      COMMON /ALUHIN/ AC,XLOAD(3),TRAILR(3),ABSHUM,DB,WB
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /OXY2/ NFUEL,MXOXMY,MAXOXY,OXYVAL(25,18,4),OXYBF(25,18,4)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /RESUL2/ RLGGAL(9),RSTGPH(9),EFRSTL(9),EFIDLE(3,9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
      COMMON /USDATA/ USNAME(4,4),NUSD,IUSD(4)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
      CHARACTER*8 NAMVEH(10),NAMNMH(5),USC8
      CHARACTER*4 NAMREG(3),NAMPOL(3)
      CHARACTER*5 NAMIDL(5)
C
      DATA NAMVEH/
     *   '  LDGV  ',
     *   ' LDGT1  ',
     *   ' LDGT2  ',
     *   '  LDGT  ',
     *   '  HDGV  ',
     *   '  LDDV  ',
     *   '  LDDT  ',
     *   '  HDDV  ',
     *   '  MC    ',
     *   'All Veh '/
C
      DATA USC8/' ______ '/
C
      DATA NAMREG/'Low ','High','Mid '/
C
      DATA NAMPOL/' HC:',' CO:','NOX:'/
C
      DATA NAMNMH/'Total   ','No-Mth  ','VOC     ',
     *            'TOG     ','NMOG    '/
      DATA NAMIDL/' Idle','NM-Id','VOCId','TOGId','NMOId'/
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
  200 FORMAT('0','User supplied',4A8,T46,A1,4A8,T78,A1/
     *       ' ',            T15,4A8,T46,A1,4A8,T78,A1)
C
C  Colons are aligned on columns 30 and 53 for the 80 column scenario header.
C
   10 IF(PRTFLG.LT.4) WRITE(IOUREP,210)
     *   ICY,NAMREG(IREJN),ALT,
     *   NOYES(IMF),TEMEXH(IP1)
  210 FORMAT(
     *'0', T2,'Cal. Year:',I5,
     *    T28,'Region: ',A4,
     *    T50,'Altitude:',F6.0,' Ft.'/
     *' ',T23,'I/M Program: ',A4,
     *    T46,'Ambient Temp:',F7.1,' (F)')
C
      IF(PRTFLG.EQ.4) WRITE(IOUREP,215)
     *   ICY,NAMREG(IREJN),ALT,
     *   NOYES(IMF),TEMEXH
  215 FORMAT(
     *'0', T2,'Cal. Year:',I5,
     *    T28,'Region: ',A4,
     *    T50,'Altitude:',F6.0,' Ft.'/
     *' ',T23,'I/M Program: ',A4,
     *    T46,'Ambient Temp:',F7.1,2(' /',F5.1),' F')
C
      IATP = ATPFLG
      IF(ATPFLG.NE.1) IATP = 2
      WRITE(IOUREP,220)
     *   NOYES(IATP),PCCN,PCHC,PCCC
  220 FORMAT(
     *' ',T17,'Anti-tam. Program: ',A4,
     *    T44,'Operating Mode:',F7.1,2(' /',F5.1))
C
      IF(RFGFLG.EQ.2) THEN
          WRITE(IOUREP,222) NOYES(RFGFLG),ASTMCL
      ELSE
          WRITE(IOUREP,223) NOYES(RFGFLG)
      ENDIF
  222 FORMAT(
     *' ',T17,' Reformulated Gas: ',A4,T47,' ASTM Class: ',A1)
  223 FORMAT(
     *' ',T17,' Reformulated Gas: ',A4)
C
      IF(IPFCU.EQ.2) THEN
        DO 12 KEYEQ=1,3
          PFRATE(KEYEQ,IPFCU)=PFRATE(KEYEQ,IPFCU)*100.
   12   CONTINUE
C
        WRITE(IOUREP,224) (PFRATE(KEYEQ,IPFCU),KEYEQ=1,3)
  224   FORMAT(
     *' ',T19,'User-Supplied Pass, Purge Fail, ',
     *'and Pressure Fail Rates (%):',2(F5.0,','),F5.0)
      ENDIF
C
      IF(LOCFLG.EQ.1) CALL OUTLAP
C
      IF(ALHFLG.GT.1) WRITE(IOUREP,225)
     *   ABSHUM,
     *   AC,DB,WB
  225 FORMAT(
     *' ',T17,'Absolute Humidity:',F6.2,
     *    T46,'AC (DB / WB):',F7.1,' (',F5.1,' /',F5.1,')')
C
C  Oxy Fuels
C
      IF(OXYFLG.EQ.2) WRITE(IOUREP,228) (SHRMKT(IFUEL),IFUEL=2,NFUEL),
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
     *'0','Veh. Type:',10A7/'+',T12,10A7/
     *' ','Veh. Spd.:',F5.1,2F7.1,7X,5F7.1/
     *' ','  VMT Mix:',3F7.3,7X,5F7.3)
C
      IF(ALHFLG.GT.1) WRITE(IOUREP,235)
     *   XLOAD,
     *   TRAILR
  235 FORMAT(
     *' ','Ext. Load:',3F7.3/
     *' ',' Trlr Tow:',3F7.3)
C
C
C Changed FORMAT 233 from F10.2,'%',F8.2,'%' - 3 DEC 92 - MLA.
      IF(LEVFLG.EQ.2) THEN
         WRITE(IOUREP,233) (ZEVTOT(I),I=1,2)
  233    FORMAT(
     *' ','ZEV Fract:',F6.2,'%',F6.2,'%')
      ENDIF
C
      WRITE(IOUREP,240)
  240 FORMAT(
     *'0','Composite Emission Factors (Gm/Mile)')
C
      DO 30 IP=IP1,IP2
        IF(IP.NE.1) GOTO 20
C
        WRITE(IOUREP,250)
     *  NAMNMH(NMHFLG),NAMPOL(IP),
     *  (EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     *  VFTP(IP)
  250 FORMAT(
     *' ',A6,A4,F6.2,9F7.2)
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
     *' ','Exhst ',A4,F6.2,9F7.2/
     *' ','Evap. ',A4,F6.2,4F7.2,21X,2F7.2/
     *' ','Refuel',A4,F6.2,4F7.2,28X,F7.2/
     *' ','Runing',A4,F6.2,4F7.2,28X,F7.2/
     *' ','Rsting',A4,F6.2,4F7.2,21X,2F7.2)
        GOTO 30
C
   20 WRITE(IOUREP,265)
     *  NAMPOL(IP),
     *  (EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     *  VFTP(IP)
  265 FORMAT(
     *' ','Exhst ',A4,
     *    F6.2,8F7.2,
     *    F7.2)
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
  275   FORMAT(
     *' ',A5,1X,A4,
     *    F6.2,8F7.2,
     *    F7.2)
   40 CONTINUE
C
   50 WRITE(IOUREP,290)
  290 FORMAT(' ',80('_'))
C
      RETURN
      END
