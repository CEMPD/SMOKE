
      SUBROUTINE MOBILE( INERR, ITEMP, MODEV, EFACT, DFACT )

cmh   SUBROUTINE MOBILE(INERR, *)
C
C  All "cmh" comments are modifications for the SMOKE modeling system
C  made by Marc Houyoux at MCNC.
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************
C
C  Called by user's program (example: Air Quality Analysis
C  System).
C
C  Calls ADJUST, BDSAVE, CONSEC, EFCALX, ONESEC, OUTNEW, OUTPUT,
C        PARSEC, QUITER, REINIT, RESTOR, REGMOD and SAVER.
C
C  The following variables are passed in argument lists or common
C  blocks and in MOBILE:
C
C    argument lists: ICY,INERR
C    common blocks:
C    /EVAL/   MEVAL
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUERR
C    /REGION/ IREJN
C    /SAVE01/ JCALL
C
C  Output on return:
C
C    parameter list: INERR
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ICYOLD   I    calendar year from previous pass of scenario loop
C  IROLD    I    region from previous pass of scenario loop
C                      0 = message has not yet been given
C                      1 = message has been given: not again!
C  IY       I    loop counting variable for interpolation runs
C  LEAST1   I    flag: 0 = OUTPUT has not yet been called.
C                      1 = OUTPUT has been called at least 1 time.
C  MEVOLD   I    calendar month from previous pass of scenario loop
C  NIY      I    number of year interations for by-month interpolation
C  NIYOLD   I    number of year interations from previous pass
C                of scenario loop
C
C
C  Notes:
C
C  For MOBILE4.1 the subroutine MOBILE was created from MAIN.
C  The /BYMYCd/ CBs were added.
C
C
      CHARACTER*40 TNAME
      CHARACTER*26 BYVTN
      CHARACTER*17 YFM
      CHARACTER*8 NAMFLG,NAMMMR
      CHARACTER*8 USNAME
      CHARACTER*7 PN,U7
      CHARACTER*8 NAMRVP
      CHARACTER*4 NAMTEM
      CHARACTER*4 NOYES,COMMA,PERIOD
      CHARACTER*4 PROJID,SCNAME
      CHARACTER*4 MONTH
      CHARACTER*3 PSEC
      CHARACTER*1 COLON,VB
      CHARACTER*1 ASTMCL
c
c  LMOS Modification 04/08/93
c
c  Add definition of ndi_di, ENVL, and M5ERR
c
      CHARACTER*1 ndi_di
cmh   CHARACTER*132 ENVL,M5ERR
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
      INTEGER ATPPGM,ATPFQT,DISTYP
      INTEGER PRSPGM,PRSFQT,PRGPGM,PRGFQT,PRSFLG,PRGFLG
      INTEGER ADDCSW,RFGSYR
      INTEGER RFGFLG,WINFLG
C
      REAL JULMYR,JANMYR,LDVIDL
      REAL LHDVM,MHDVM,HHDVM,LHDRG,MHDRG,HHDRG
      REAL LEVSTD, LEVIMP, LDTAB
      LOGICAL RFGON
C
cmh add includes
      INCLUDE 'MBDIMS3.EXT'
cmh
cmh add declarations
      REAL    EFACT( MXNDET, NVTYPE )
      REAL    DFACT( MXDET, NVTYPE )
      INTEGER MODEV                ! Unit number of MOBILE5a input file
      INTEGER ITEMP   ! When == 0, tmpr is only change from previous call
      LOGICAL WFLAG
cmh
cmh add commons
      COMMON /WARNING/ WFLAG
cmh

      COMMON /ALUHAC/ ACEQIP(8,3,2),ACCF(2,3,2,3,2)
      COMMON /ALUHIN/ AC,XLOAD(3),TRAILR(3),ABSHUM,DB,WB
      COMMON /ALUHLT/ XLCF(3,6,3,2),TWCF(3,6,3,2)
      COMMON /ALUOUT/ ALHRET(25,4,3)
      COMMON /ALUPAC/ MYGACE(8,3,2),MYGACC(2,3,2)
      COMMON /ALUPLT/ MYGXLC(6,3,2),MYGTWC(6,3,2)
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
      COMMON /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP(8)
      COMMON /ATPAR3/ AIR(2),CAT(2),ECK(3),TNK(2),EGR(2)
      COMMON /ATPAR4/ CAN(2),PCV(2),CAP(2),EFF(2,3,2,15)
      COMMON /BASEQ1/ ERBZML(25,3,8,2),ERBDR(25,3,8,2),ERB50K(18,3,3,2)
      COMMON /BASEQ2/ MYGERB(25,3,8,2),MAXERB
      COMMON /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW,KINKS
      COMMON /BASEQ4/NEWPAR(5,100),BERNEW(3,100),NEWFIT(100),NEWCT,MAXCT
      COMMON /BASEQ5/ ERUZML(18,3,8,2),ERUDR(18,3,8,2),ERU50K(18,3,3,2)
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASEQ7/ IPMOD(3)
      COMMON /BASEQ8/ IHDVMY(23,2),HDVZML(23,2),HDVBDR(23,2)
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /BASE11/ SLOIM1, SLOIM2
      COMMON /BASE12/ CAPIV1, CAPIV2
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC3/BYEVAP(25,8),BYRUNL(25,4),BYREFL(25,4),BYRSTL(25,8)
      COMMON /BYMYC4/ BYVTN(8),YFM,PN(10),U7,PSEC(3),VB
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
      COMMON /CEVCOM/ USRTPD(8),USRMPD(8)
      COMMON /CEVBMY/ BMYMPD(25,8),BMYTPD(25,8),BMTPDF(25,8)
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /CITPAR/ SCNAME(4)
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /CRACOM/ MAXCRA,CRAVAL(3,8,2),MYGCRA(3,8,2)
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
      COMMON /EMECAL/ EM3W(3,3,3),EMOX(3,3,3),EMEGR(4,3,3),EMI(7,3,3)
      COMMON /ETPPAR/ ETP(7,3),NEWETP,IETPYR(4),ETPPIV(4),NETPYR
      COMMON /EVADU1/ DUEQ(3,2,3),EFDU(3,2,5,5)
      COMMON /EVADU2/ HIDU(3,2),R2DU(2,2)
      COMMON /EVAHIA/ HIHS(3,2),HIADJ
      COMMON /EVAHS1/ HSEQ(3,2,9),EFHS(3,2,5,5)
      COMMON /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ(2)
      COMMON /EVAPGR/ EVP(4),GREVP(4,9),VGREVP(4),PARD(3),GRPD(3,9)
      COMMON /EVAPHD/ HDSAL,HDWGT(2),IVVMAP(5)
      COMMON /EVASTD/ MYGSTD(5,5),MAXSTD,KSTD81(4,3)
      COMMON /EVPDAT/ EVPAER(3,2),EVPTGS(3,2,4),CCEMI(7,4)
      COMMON /FID/ MAXFID,MYFID(16,5),TOGCF(16,5,3),VOCCF(16,5,3)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /GSFCOM/ MAXGSF,GSFRAC(22,8,2),MYGSF(22,8,2)
      COMMON /HDCCOM/ HDCFAC(37,2),LOWYR,MAXCIG
      COMMON /HDDMAR/ C2BDVM(25),LHDVM(25),MHDVM(25),HHDVM(25)
      COMMON /HDDREG/ C2BDRG(21),LHDRG(21),MHDRG(21),HHDRG(21)
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /IDLEQ2/ MYGIDL(15,3,8,2),MAXIDL
      COMMON /IDLEQ3/ UTIDOF(14,3,2),LDVIDL(25)
      COMMON /ICR4VA/ CRD4VA(24,3,18)
      COMMON /ICR4VB/ CRD4VB(24,3,18)
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR4/ CUTHC(2),CUTCO(2),CUTNO(2)
      COMMON /IMPAR5/ MYGIM2(3,3),LBIM4P,CRHDGV(2,2),DISCNT(3)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IMPAR7/ TNAME,IREAD
      COMMON /IMPAR8/ ICALC
      COMMON /IM12HC/ CR12HC(19,20,5,2)
      COMMON /IM12CO/ CR12CO(19,20,5,2)
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /INJECT/ TBI(13,4),PFI(13,4),MAXFIY,MINFIY
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /IVPCOM/ IVPTRL(8),IVPTRT(8),IVPTRA(8),IVPTRB(8)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MPDCOM/ MAXMPD,VALMPD(2,8),MYGMPD(2,8)
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
      COMMON /MYUB1/  MYTGS(11,6,4),MYEGR(4,4),MYVTS(3,2,4)
      COMMON /MYUB2/  MYCCEI(7,4)
      COMMON /NAMES1/ NAMFLG(20),NAMMMR(2,3,8),NAMRVP(4),NAMTEM(2)
      COMMON /NMETH1/ ZDMTH1(2,3,2,7,3,2),ZDMTH2(2,7,2,2),ZDMTH3(7,3,2)
      COMMON /NMETH2/ MYGMTH(7,8,2),MTHMYG
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OPMOD1/ BFRCOF(8,8,3,6),BFSTEP(8,12,3,3),BFR50K(4,12,3,3)
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
      COMMON /OPMOD3/ BFMILE,KINKBF,BFGT50
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /OXY2/ NFUEL,MXOXMY,MAXOXY,OXYVAL(25,18,4),OXYBF(25,18,4)
      COMMON /PROJEC/ PROJID(20)
      COMMON /QUITXQ/ N1QUIT
      COMMON /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG(4)
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS(4)
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /REGISF/ GSFVCT(8),USRGSF(25,2)
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /RESUL2/ RLGGAL(9),RSTGPH(9),EFRSTL(9),EFIDLE(3,9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
      COMMON /RLCOM1/ ROADFE(32,4),DISPL,SPILL,OBED,OBES,OBDF(4)
      COMMON /RLCOM2/ IOBMY,IVOB(4),IS2SY,NPHASE,S2EFF(4)
      COMMON /RLCOM3/ RLRATE(25,4)
      COMMON /RUNLS1/ TLEMI(6,3,4,4,4,4),TLVMT(6),TLVMTU(6),MAXTL
      COMMON /RUNLS2/ TLFAIL(6,4,4)
      COMMON /RUNLS3/ MYGRUL(4,4),RULRVP(4),RULTEM(4),RULSPD(3),MAXRUL
      COMMON /RESTLS/ PTGRST(13,5,2),RESTA(2),RESTB
      COMMON /RVPEX1/ R9G1(3,2),RG23(3,3,3,2),TLL(2)
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /RVPNAT/ RVPLIM(2),RVPDIU(2),PFUL
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SCENE2/ FCC,FCN,FHC,FHN
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,3),SK4C(2,14,3,3)
      COMMON /SPEED2/ MYGSP1(18,2,3),MYGSMC(18,2,3),MYGSP2(3,3)
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED5/ MYHSPC(2,3),MAXHSP(3)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /SPEED7/ GPBSCO(3,3,4)
      COMMON /SPEED8/ SIXT5(11,3,3)
      COMMON /SIZCAL/ ATR(5,2),TGSUSE,CSIZE(12,2),CSAE(12)
      COMMON /SIZPAR/ NTERMS(4),INDXCS(5,4)
      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
      COMMON /TAMEQ1/ TAMZML(9,4,3,2),TAMDR(9,4,3,2),TAMA50(9,4,3,2)
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
      COMMON /TAMEQ3/ CUMCAP,CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI
      COMMON /TAMEQ4/ BTR(9,2)
      COMMON /TAMID1/ EMI3W(3,3,3),EMIOX(3,3,3),EMIDL(3,4,3)
      COMMON /TAMID2/ MYIDLE(3,3),IEK
      COMMON /TAMID3/ OFFIDL(3,25,4)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TAMOU2/ TOB(25,4),GCONLY(25,4)
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
      COMMON /TECHFL/ CONSTF
      COMMON /TEMPC1/TTFC(2,3,9,3),TTFD(3,3,3),TT4(2,6,3),TT8(2,3,4,3)
      COMMON /TEMPC2/ MAXTCF,MYGTCF(9,3,5),MYGTFD(3),ISHIFT
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPC4/ MYGCOO(3),LOWCO,COOTFC(3),COOTFD(3)
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(3,2)
      COMMON /TEMPC6/ CSUB(3,2),CMUL(3,2),PIVPCT(3,3,2)
      COMMON /TEMPC7/ ICOLD,ADDCSW,MULCSW
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
      COMMON /TIER1S/ TR1SL1(6,3,3)
      COMMON /TPDCOM/ TPDCO(2,8)
      COMMON /USDATA/ USNAME(4,4),NUSD,IUSD(4)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
      COMMON /SAVE01/ JCALL,VA001(15,3,8,2),VA002(15,3,8,2)
      COMMON /SAVE04/ VA007(25,8),VA008(25,8)
      COMMON /SAVE05/ VA009(9,4,3,2),VA010(9,4,3,2)
      COMMON /SAVE06/ VA011(8),VA012(7,3),IVA013(4),VA014(4)
      COMMON /EVAL/ MEVAL,MONTH(12)
c
c  LMOS Modification 04/08/93
c
c  Add LMOS and SCENE3 common blocks
c
      COMMON /LMOS/ EFDNL(9),EFEVND(9),VDNL,DNLGPM
      COMMON /SCENE3/ ndi_di
c
C
cmh add 14 lines
      IOUGEN = MODEV
      IOUGET = 5
      IOUREP = 7
      IOUASK = 6
      IOUERR = 6
      IOBAT = 19

      WFLAG = .FALSE.     ! Initialize warning message write to off

      IF( ITEMP .EQ. 1 ) THEN
          IREAD = 0 ! Reset MOBILE5's indicator for reading TECH12 file
          ITEMP = 0 ! Indicate that only temperature is changing
          WFLAG = .TRUE.  ! Turn warning messages on
      ENDIF
cmh
      IF(JCALL.EQ.0) GOTO 10
C
C  Subroutine call, setup run.
C
      CALL REINIT
      IF(JCALL.EQ.1) CALL BDSAVE
      IF(JCALL.EQ.2) CALL RESTOR
      JCALL=2
C
C  Initialize non-COMMON variables (INERR & local variables).
C
   10 INERR=0
      ICYOLD=-1
      IROLD=-1
      MEVOLD=-1
      NIYOLD=-1
      LEAST1=0
C
C  Control and One-Time Data Sections:
C
C  Get the execution control flags and parameters, and then use them to read
C  in data, if any, that is read in once and applied to all scenarios.  Note:
C  alternate RETURN1 is the fatal read error exit: execution stops.
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,110)
  110 FORMAT(' Reading information.')
C
      CALL CONSEC(*99,*999)
      CALL ONESEC(INERR,*99)
C
C  ljn hc=3 3/5/93 Warning MSG in QUITER
C
      IF(HCFLAG.EQ.3.AND.OUTFMT.NE.3.AND.OUTFMT.NE.5)
     *                    CALL QUITER(0.,0,153,INERR)
C
C  only get error diagnostics out of the run.
C
c  LMOS Modification 04/08/93
c
c  Open the errors file for each input file separately to aviod creating a HUGE errors
c  file during testing. Also, eliminate the count of errors encountered.
c  This is needed because the temperature matrix is wide enough to cause errors.
c  Therefore, in this modified version, if this count exceeds 50, the run will
c  NOT abort.
c
cmh20 CLOSE(IOUERR)
cmh   CALL GETENV("EMS_LOC",ENVL)
cmh   M5ERR=ENVL(1:LNBLNK(ENVL))//'/M5ERRORS_WARNING_REPORT'
cmh   OPEN(IOUERR,FILE=M5ERR,STATUS='UNKNOWN')
cmh   CALL PARSEC(ICY,INERR,*98)
   20 CALL PARSEC(ICY,INERR,*98)
C
C  50+ errors is (arbritrary) run limit => stop execution.  Otherwise,
C  one or more input processing errors => abort current scenario & then do
C  input processing for error diagnostics on the remaining scenarios.  INERR
C  is never reset to 0, but instead accumulates until the run ends or the
C  cutoff value (50) is reached.
C
c
c  LMOS Modification 04/08/93
c  Comment out this check
c
c      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,115)
  115 FORMAT(' Performing calculations.')
C
C  For interpolated (non Jan 1st) emissions, loop through ICY and ICY+1,
C  save intermediate results, then interpolate to the evaluation month.
C
      NIY=2
      IF(MEVAL.EQ.1.OR.ICY.EQ.2020) NIY=1
C
      DO 35 IY=1,NIY
      ICY1=ICY+IY-1
 
      IF(NIY.GT.1 .OR.
     *   NIYOLD.NE.NIY .OR.
     *   ICYOLD.NE.ICY1 .OR.
     *   IROLD.NE.IREJN .OR.
     *   MEVOLD.NE.MEVAL) CALL REGMOD(ICY1,INERR)
c
c  LMOS Modification 04/08/93
c  Comment out this check
c
c      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      ICYOLD=ICY1
      IROLD=IREJN
      MEVOLD=MEVAL
      NIYOLD=NIY
C
      CALL EFCALX(ICY1,INERR,IY)
      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(NIY.GT.1) CALL SAVER(IY)
C
  35  CONTINUE
C
      IF(NIY.GT.1) CALL ADJUST
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,120)
  120 FORMAT(' Preparing output.')
C
cmh   CALL OUTPUT(ICY)
cmh
cmh add  6 lines

C.......  Initialize EFACT and DFACT
      DO JJ = 1, NVTYPE
          DO II = 1, MXNDET
              EFACT( II,JJ ) = 0.0
          ENDDO
          DO II = 1, MXDET
              DFACT( II,JJ ) = 0.0
          ENDDO
      ENDDO

      DO JJ = 1,NVTYPE

          EFACT( 1,JJ ) = EFFTP ( 2,JJ )   ! Exhaust CO
          EFACT( 2,JJ ) = EFFTP ( 3,JJ )   ! Exhaust NOX
          EFACT( 3,JJ ) = EFEXH ( JJ )     ! Exhaust VOC
          EFACT( 4,JJ ) = EFEVAP( JJ )     ! Evaporative VOC
          EFACT( 5,JJ ) = EFRUNL( JJ )     ! Running Loss VOC
          EFACT( 6,JJ ) = EFRSTL( JJ )     ! Resting Loss VOC

cmh       EFACT( 7,JJ ) = EFREFL( JJ )     ! Refueling Loss VOC
cmh       EFACT( 8, JJ ) = RLGGAL( JJ )     ! Refueling loss VOC g/gal
cmh       EFACT( 11,JJ ) = RSTGPH( JJ )     ! Resting Loss VOC g/hr

      ENDDO

      DO JJ = 1,4
          DFACT( 1, JJ ) = GREVP ( 2,JJ )   ! Weighted Diurnal Evap VOC
          DFACT( 2, JJ ) = EFDNL (   JJ )   ! LMOS Diurnal Evap VOC
          DFACT( 3, JJ ) = GREVP ( 1,JJ )   ! Hot Soak Evap VOC
          DFACT( 4, JJ ) = GREVP ( 4,JJ )   ! Crankcase Evap VOC
      ENDDO

      DFACT( 2, 8 ) = EFDNL ( 8 )           ! LMOS Diurnal for MC's
      
cmh end add
cmh 
      LEAST1=1
      GOTO 20
C
   30 CALL QUITER(0.0,0,28,INERR)
C
C  If user is supplying replacement basic FTP ef parameter records, always
C  print the records contents and disposition table, even when input errors
C  have prevented the normal output routines from being called.  The only
C  exception is when a READ error / end-of-file occurs in CONSEC or ONESEC.
C
C  This section can be expanded to print other optional input.
C
   98 IF(LEAST1.EQ.0.AND.
     *  (NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)) CALL OUTNEW
C
   99 RETURN
C
C  End-of-file on 1st read of call => input for call missing entirely.
C  Calling program can branch to attach a new input file or terminate
C  calls to MOBILE.
C
cmh  999 RETURN 1
  999 RETURN

94300 FORMAT ( A, I5, A, A )
      END
