C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EFCALX(ICY,INERR,IY)
C
C  EFCALX calculates composite emission factors.
C
C  Called by MOBILE.
C
C  Calls BEF, BIGCFX, EFNFTP, GETCUM, GSFIDX, HCCALX, QUITER, REFUEL,
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
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor 
C  calculation  function.
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4 
C  calculation to remove improperly making a SALHCF adjustment to by-model
C  Sep-21-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar
C  year removed IY from the call to TAMPER
C  Oct-11-1994 @ CSC-bsg Request 467 Made correction to second call to 
C  EFNFTP function, the order of parameters was incorrect
C
C
      INCLUDE 'BYMYC1.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC5.I'
      INCLUDE 'CITRV2.I'
      INCLUDE 'CUMCOM.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'GSFCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'MYRSAV.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'RVPEX2.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'VMXCOM.I'
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
      IF(TAMFLG.GT.0) CALL TAMPER(ICY)
      IF(TAMFLG.EQ.0) CALL QUITER(0.,0,119,INERR)
C
C  21-Nov-1994 @ CSC-ked Subtask 0-411
C     Eliminate the printing of warning message 154 twice by checking IY
C
      IF(PRTFLG.EQ.1.OR.PRTFLG.EQ.4) CALL REFUEL(ICY,IY)
      CALL RVPEXH(ICY)
C
C  Call "BIG" correction factor (cf) driver to get cf's for speed, temperature,
C  operating mode, tampering, air conditioning, trailer towing, humidity, extra
C  load and by-bag RVP, which cf's apply depending on user's input selections.
C
      CALL BIGCFX(ICY)
C
C  Finally ready for ef calculations - only for user selected pollutants.
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
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor
C  calculation function.  So had to move TF adjustment to end of section.
C
      EXHWGT=SALHCF(IDX,IP,IV)
      COMPEF=BEF(MY,IDX,ICY,IP,IV,VMTAGE,IY)*EXHWGT
      IF(IV.LE.4) BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*EXHWGT*TF(IDX,IV)
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
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor 
C  calculation function.  So had to move TF adjustment to end of section.
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4 
C  calculation to remove improperly making a SALHCF adjustment to by-model
C
   60 COMPEF=COMPEF+EFNFTP(MY,IV,IP)
C
C  Oct-11-1994 @ CSC-bsg Request 467 Made correction to second call to 
C  EFNFTP function, the order of parameters was incorrect
C
      BYBEF4(IP,JDX,IV)=BYBEF4(IP,JDX,IV)+EFNFTP(MY,IV,IP)
      BYFER(IP,JDX,IV)=COMPEF*TF(IDX,IV)
      EFFTP(IP,IV)=EFFTP(IP,IV)+COMPEF*TF(IDX,IV)
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
