C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EFCALX(ICY,INERR,IY)
C
C  EFCALX calculates composite emission factors.
C
C  Called by MOBILE.
C
C  Calls BEF, BIGCFX, GETCUM, GSFIDX, HCCALX, QUITER, REFUEL,
C        RVPEXH and TAMPER.
C
C  Input on call:
C
C    parameter list: ICY,INERR,IY
C    common blocks:
C    /CITRV2/ RVPX
C    /CUMCOM/ CUMMIL
C    /FLAGS1/ TAMFLG
C    /FLAGS3/ RLFLAG
C    /FLAGS4/ PRTFLG,IDLFLG
C    /GSFCOM/ GSFRAC
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /MYRCAL/ TF
C    /MYRSAV/ NEWCUM
C    /REGION/ IREJN
C    /RVPEX2/ RVPCF
C    /VMXCOM/ VMTMIX
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC5/ BYIMRE,REDSUM
C    /RESUL1/ EFFTP
C    /RESUL3/ VFTP
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COMPEF   R    composite emission factor for IDX: summed across IDX to get
C                EFFTP for the given IP & IV
C  EXHWGT   R    multiplicative weights to be applied to the exh ef returned
C                by BEF: speed & "BIGALH" cf, travel fraction & RVP cf
C                given vehicle type and model year
C  RVOPCF   R    exhaust HC & CO cf for RVP and open loop technology, IV=1,4
C  VMTAGE   R    see parameter dictionary
C
C
C  Notes:
C
C  For MOBILE4 EFCALX was modified to save by model year
C  computations for the by model year tables.
C  For MOBILE4.1 Oxygenated and certification fuels support was
C  removed and in particularly the OPENLP correction.
C  March 10, 1993 IY was added to the parameter list to guarantee
C  the 1990 inventory emission factors remain the same in newer
C  versions of MOBILE as in MOBILE5.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
      REAL JULMYR,JANMYR
C
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /GSFCOM/ MAXGSF,GSFRAC(22,8,2),MYGSF(22,8,2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
      IF(NEWCUM.NE.0) CALL GETCUM
      NEWCUM=0
C
C  CUMMIL, a parameter array, is recast here for the BYMY tables as BYCUM,
C  with zero replacing the potential CUMMIL, when TF=0.0.
C
      DO 10 IV=1,MAXVEH
      DO 10 JDX=1,MAXYRS
      IDX=(MAXYRS+1)-JDX
      IF(TF(IDX,IV).EQ.0.0) BYCUM(JDX,IV)=0.0
      IF(TF(IDX,IV).GT.0.0) BYCUM(JDX,IV)=CUMMIL(JDX,IV)
   10 CONTINUE
C
C  Setup for composite light duty gas truck.
C
      VMLDGT=VMTMIX(2)+VMTMIX(3)
      IF(VMLDGT.EQ.0.0) VMLDGT=1.
C
C  Calculate this scenario's tampering and refueling offsets and RVP correction
C  factors.
C
      IF(TAMFLG.GT.0) CALL TAMPER(ICY,IY)
      IF(TAMFLG.EQ.0) CALL QUITER(0.,0,119,INERR)
      IF(PRTFLG.EQ.1.OR.PRTFLG.EQ.4) CALL REFUEL(ICY)
      CALL RVPEXH(ICY)
C
C  Call "BIG" correction factor (cf) driver to get cf's for speed, temperature,
C  operating mode, tampering, air conditioning, trailer towing, humidity, extra
C  load and by-bag RVP, which cf's apply depending on user's input selections.
C
      CALL BIGCFX(ICY)
C
C  Finally ready for ef calculations - only for user selected pollutants.
C
C  First do idle efs (after interpolation zap, move to precede BIGCFX call).
C
      IF(IDLFLG.EQ.2) CALL IDLCAL(ICY,VMLDGT,INERR)
C
C  Now do FTP efs.
C
      IF(PRTFLG.EQ.1.OR.PRTFLG.EQ.4)
     *    CALL HCCALX(ICY,VMLDGT,INERR,IY)
C
      DO 90 IP=2,MAXPOL
      IF(PRTFLG.NE.IP.AND.PRTFLG.NE.4) GOTO 90
      VFTP(IP)=0.0
      EFFTP(IP,9)=0.0
C
      DO 80 IV=1,MAXVEH
      EFFTP(IP,IV)=0.0
      IF(IV.LE.4) REDSUM(IP,IV)=0.0
C
C  Zero out BYMY save areas before VMTMIX check.  Since VMTMIX can be altered
C  by scenario, a VMTMIX(i) > 0.0 case followed by a VMTMIX(i) = 0.0 case
C  would lead to the latter scenario using the former scenario's BYMY data for
C  IV=i.
C
      DO 30 JDX=1,MAXYRS
      BYBEF4(IP,JDX,IV)=0.0
      BYFER(IP,JDX,IV)=0.0
      IF(IV.GT.4) GOTO 30
      BYTAM(IP,JDX,IV)=0.0
      BYIMCR(IP,JDX,IV)=0.0
      BYIMRE(IP,JDX,IV)=0.0
   30 CONTINUE
C
      IF(VMTMIX(IV).EQ.0.0) GOTO 80
      MY=0
C
      DO 70 IDX=1,MAXYRS
C
C  If no mileage for this model year / vehicle class case, skip ef calculation.
C
      IF(TF(IDX,IV).LE.0.0) GOTO 70
C
C  Initialize emissions to zero.
C
      COMPEF=0.0
C
      JDX=(MAXYRS+1)-IDX
      VMTAGE=CUMMIL(JDX,IV)/10000.
      MY=ICY+IDX-MAXYRS
C
C  Combine speed and aluh correction factor and the travel fraction
C  into a single weight to correct the exhaust ef returned from BEF.
C
      EXHWGT=SALHCF(IDX,IP,IV)*TF(IDX,IV)
      COMPEF=BEF(MY,IDX,ICY,IP,IV,VMTAGE)*EXHWGT
      IF(IV.LE.4) BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*EXHWGT
C
C  Diesel classes have no further corrections.
C
      GOTO(40,40,40,40,60,60,60,60),IV
C
C  For LDGV/T & HDGV, the RVP CF applied may be 1.0, either because there is no
C  correction or because the correction is by bag (& has already been applied
C  via OMTCF).  Only RVP MYG 1 uses RVPCF.
C
C  LDGV/T and HDGV NOx have a RVPCF, but no open loop technology credit.
C
   40 IF(IP.EQ.2) GOTO 50
      BYBEF4(IP,JDX,IV)=BYBEF4(IP,JDX,IV)*RVPCF(IDX,IP,IV)
      BYTAM(IP,JDX,IV)=BYTAM(IP,JDX,IV)*RVPCF(IDX,IP,IV)
      BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*RVPCF(IDX,IP,IV)
      REDSUM(IP,IV)=REDSUM(IP,IV)+BYIMRE(IP,JDX,IV)
      COMPEF=COMPEF*RVPCF(IDX,IP,IV)
      GOTO 60
C
C  LDGV/T and HDGV CO get both RVP exhaust and open loop technology cfs
C  calculated and applied.
C
   50 RVOPCF=RVPCF(IDX,IP,IV)
      BYBEF4(IP,JDX,IV)=BYBEF4(IP,JDX,IV)*RVOPCF
      BYTAM(IP,JDX,IV)=BYTAM(IP,JDX,IV)*RVOPCF
      BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*RVOPCF
      REDSUM(IP,IV)=REDSUM(IP,IV)+BYIMRE(IP,JDX,IV)
      COMPEF=COMPEF*RVOPCF
C
   60 BYFER(IP,JDX,IV)=COMPEF
      EFFTP(IP,IV)=EFFTP(IP,IV)+COMPEF
   70 CONTINUE
C
C  Sales fractions = 1.0 for all 4 HDV (HDG,HDD x L,H alt) groups and all MC's.
C
      IF(MY.GT.0.AND.(IV.NE.4.AND.IV.LT.7).AND.
     *   (EFFTP(IP,IV).GT.0.0.AND.GSFIDX(ICY,MY,IV,IREJN).EQ.0.0))
     *   CALL QUITER(0.,IV,67,INERR)
C
      IF(EFFTP(IP,IV).LE.0.0) CALL QUITER(0.,IV,68,INERR)
C
C
      VFTP(IP)=VFTP(IP)+EFFTP(IP,IV)*VMTMIX(IV)
   80 CONTINUE
C
C  Assign composite light duty gas truck emission factor.
C
      IF(EFFTP(IP,2).GT.0.0.AND.EFFTP(IP,3).GT.0.0) EFFTP(IP,9)=
     *   (EFFTP(IP,2)*VMTMIX(2)+EFFTP(IP,3)*VMTMIX(3))/VMLDGT
   90 CONTINUE
C
      RETURN
      END
