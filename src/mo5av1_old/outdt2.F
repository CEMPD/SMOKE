C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDT2(ICY)
C
C  OUTDT2 writes scenario parameters and output on unit IOUREP using output
C  format (OUTFMT) 2.
C
C  Called by OUTPUT.
C
C  Calls OUTPOL and FORTRAN library routine IFIX.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /FLAGS4/ IDLFLG,HCFLAG
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXVEH
C    /REGION/ IREJN,ALT
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,EFREFL,EFRUNL
C    /RESUL2/ EFIDLE,EFRSTL
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VIDLE
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /TEMPS/  TEMEXH
C    /VMXCOM/ VMTMIX
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IAMBT    I    integer equivalent of exhaust (ambient) temperature
C  ICY2     I    last 2 digits of calendar (scenario evaluation) year
C  IPIDLE   I    output pollutant code POL for idle efs = IP + 3
C
C  Notes:
C
C  OUTDT2 was modified for MOBILE4 In-House Version 02 to update the array name
C  EFLOSS to EFREFL.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /RESUL2/ RLGGAL(9),RSTGPH(9),EFRSTL(9),EFIDLE(3,9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
      CALL OUTPOL(IP1,IP2)
C
      ICY2=ICY-((ICY/100)*100)
C
      DO 10 IP=IP1,IP2
      IAMBT=IFIX(TEMEXH(IP)+.5)
      WRITE(IOUREP,200)
     * IREJN,ICY2,IAMBT,PCCN,PCHC,PCCC,IP,
     * (EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     * VFTP(IP),VMTMIX
C
  200 FORMAT(' ',I1,1X,I2,I4,1X,3F6.1,1X,I1,10F7.2,8(1X,F4.3))
C
C  Only write data for HC components, since already have id fields on total
C  HC line.  Whenever > 1 output line per scenario, postprocessor needs to
C  read id fields only once per scenario.  The ef fields will always be in the
C  same columns on each output record, resulting in only two postprocessor
C  read formats - 1 to get id + ef fields & 1 to get ef fields only.
C
      IF(IP.EQ.1.AND.HCFLAG.GE.2) WRITE(IOUREP,210)
     *   (EFEXH(IV),IV=1,3),EFEXH(9),(EFEXH(IV),IV=4,MAXVEH),VEXH,
     *   (EFEVAP(IV),IV=1,3),EFEVAP(9),EFEVAP(4),EFEVAP(8),VEVAP,
     *   (EFREFL(IV),IV=1,3),EFREFL(9),EFREFL(4),VLOSS,
     *   (EFRUNL(IV),IV=1,3),EFRUNL(9),EFRUNL(4),VRUNLS,
     *   (EFRSTL(IV),IV=1,3),EFRSTL(9),EFRSTL(4),EFRSTL(8),VRSTLS
C
  210 FORMAT(' ',28X,'X',10F7.2/
     *       ' ',28X,'V',5F7.2,21X,2F7.2/
     *       ' ',28X,'R',5F7.2,28X,F7.2/
     *       ' ',28X,'T',5F7.2,28X,F7.2/
     *       ' ',28X,'S',5F7.2,21X,2F7.2)
C
   10 CONTINUE
C
      IF(IDLFLG.EQ.1) RETURN
      DO 20 IP=IP1,IP2
      IPIDLE=IP+3
      WRITE(IOUREP,220)
     * IPIDLE,(EFIDLE(IP,IV),IV=1,3),EFIDLE(IP,9),
     * (EFIDLE(IP,IV),IV=4,MAXVEH),VIDLE(IP)
  220 FORMAT(' ',28X,I1,10F7.2)
   20 CONTINUE
C
      RETURN
      END
