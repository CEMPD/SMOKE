C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE TFCALX(ICY,INERR)
C
C  TFCALX computes the travel fractions, VMT mix, and registration mix.
C
C  Called by REGMOD.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: ICY,INERR
C    common blocks:
C    /FLAGS1/ VMFLAG
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ JANMYR,TFMYM
C    /REGISF/ GSFVCT
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /MYRCAL/ TF
C    /VMXCOM/ REGMIX,TFNORM,VMTMIX
C
C  Local array subscripts:
C
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  TOTREG   R    total registrations: used to normalize the registration mix
C  TOTVMT   R    total vmt: used to normalize the default vmt mix
C
C  Notes:
C
C  None.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      REAL JANMYR
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /REGISF/ GSFVCT(8),USRGSF(25,2)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
C
      DO 20 IV=1,MAXVEH
      TFNORM(IV)=0.0
C
      DO 10 JDX=1,MAXYRS
      TFNORM(IV)=TFNORM(IV)+TFMYM(JDX,IV)*JANMYR(JDX,IV)
   10 CONTINUE
   20 CONTINUE
C
C  Calculate normalized VMTMIX only if using default vmt data (VMFLAG = 1).
C
      IF(VMFLAG.GT.1) GOTO 50
      TOTVMT=0.0
      DO 30 IV=1,MAXVEH
C
C For all IV that are pre-1982 ICY, use the 1982 VCOUNT data.
C ICY indexes into VCOUNT, ICY requires to be offset by "1981"
C the start of the calander years + 1  = "1982"
C
      IF(ICY.LE.1981)
     *    VMTMIX(IV)=VCOUNT(1,IV)*GSFVCT(IV)*TFNORM(IV)
C
      IF(ICY.GT.1981)
     *    VMTMIX(IV)=VCOUNT(ICY-1981,IV)*GSFVCT(IV)*TFNORM(IV)
C
      TOTVMT=TOTVMT+VMTMIX(IV)
   30 CONTINUE
C
      DO 40 IV=1,MAXVEH
      VMTMIX(IV)=VMTMIX(IV)/TOTVMT
   40 CONTINUE
C
   50 DO 70 IV=1,MAXVEH
      IF(TFNORM(IV).LE.0.0.AND.VMTMIX(IV).GT.0.0)
     *   CALL QUITER(0.,IV,48,INERR)
      IF(TFNORM(IV).LE.0.0) TFNORM(IV)=1.
C
      DO 60 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      TF(IDX,IV)=(TFMYM(JDX,IV)*JANMYR(JDX,IV))/TFNORM(IV)
   60 CONTINUE
   70 CONTINUE
C
C  Calculate normalized REGMIX.
C
      TOTREG=0.0
      DO 80 IV=1,MAXVEH
      IF(ICY.LE.1981)
     *    REGMIX(IV)=VCOUNT(1,IV)*GSFVCT(IV)
      IF(ICY.GT.1981)
     *    REGMIX(IV)=VCOUNT(ICY-1981,IV)*GSFVCT(IV)
      TOTREG=TOTREG+REGMIX(IV)
   80 CONTINUE
C
      DO 90 IV=1,MAXVEH
      REGMIX(IV)=REGMIX(IV)/TOTREG
   90 CONTINUE
C
      RETURN
      END
