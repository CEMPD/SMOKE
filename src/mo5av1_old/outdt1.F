C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTDT1(ICY)
C
C  OUTDT1 writes scenario parameters and output on unit IOUREP using output
C  format (OUTFMT) 1.
C
C  Called by OUTPUT.
C
C  Calls OUTPOL and FORTRAN library routine IFIX.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C
C    /CITPAR/ SCNAME
C    /FLAGS3/ OUTFMT
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
C  IDLEF    A    output idle pollutant string
C  IPIDLE   I    output idle pollutant code
C  NAMEF6   A    output pollutant string
C
C
C  Notes:
C
C  OUTDT1 was modified for MOBILE4 In-House Version 02 to update the array name
C  EFLOSS to EFREFL.
C  Modified for MOBILE5 to add output pollutant and idle pollutant strings.
C
C
      CHARACTER*4 SCNAME
      CHARACTER*6 NAMEF6(3)
      CHARACTER*11 IDLEF(3)
      INTEGER PRTFLG,HCFLAG,COLDFG,OUTFMT
C
      COMMON /CITPAR/ SCNAME(4)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
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
      DATA IDLEF/
     *'"Idle HC ",','"Idle CO ",','"Idle NOx",'/
C
      DATA NAMEF6/
     *'"HC ",','"CO ",','"NOx",'/
      CALL OUTPOL(IP1,IP2)
C
      ICY2=ICY-((ICY/100)*100)
C
      DO 10 IP=IP1,IP2
      IAMBT=IFIX(TEMEXH(IP)+.5)
      IF(OUTFMT.EQ.1) WRITE(IOUREP,200)
     * IREJN,ICY2,SPD,IAMBT,PCCN,PCHC,PCCC,ALT,
     * IP,(EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),(EFFTP(IP,IV),IV=4,MAXVEH),
     * VFTP(IP),VMTMIX,SCNAME
C
  200 FORMAT(
     *' ',I1,1X,I2,1X,7(F4.1,'/'),F4.1,I4,3F6.1,1X,F6.0,
     * 1X,I1,10F8.3,8F6.3,1X,4A4)
C
      IF(OUTFMT.EQ.6) WRITE(IOUREP,210)
     * NAMEF6(IP),ICY,(EFFTP(IP,IV),IV=1,3),EFFTP(IP,9),
     * (EFFTP(IP,IV),IV=4,MAXVEH),VFTP(IP),
     * IAMBT,PCCN,PCHC,PCCC,ALT,IREJN,SPD,
     * VMTMIX,SCNAME
C
  210 FORMAT(
     * A6,I4,',',10(F8.3,','),I4,',',3(F6.1,','),F6.0,
     * ',',I1,',',7(F4.1,','),F4.1,',',8(F6.3,','),4A4)
C
C  See OUTDT2 comment on id versus ef fields.
C
      IF(IP.EQ.1.AND.HCFLAG.GE.2.AND.OUTFMT.EQ.1) WRITE(IOUREP,220)
     *   (EFEXH(IV),IV=1,3),EFEXH(9),(EFEXH(IV),IV=4,MAXVEH),VEXH,
     *   (EFEVAP(IV),IV=1,3),EFEVAP(9),EFEVAP(4),EFEVAP(8),VEVAP,
     *   (EFREFL(IV),IV=1,3),EFREFL(9),EFREFL(4),VLOSS,
     *   (EFRUNL(IV),IV=1,3),EFRUNL(9),EFRUNL(4),VRUNLS,
     *   (EFRSTL(IV),IV=1,3),EFRSTL(9),EFRSTL(4),EFRSTL(8),VRSTLS
C
  220 FORMAT(' ',74X,'X',10F8.3/
     *       ' ',74X,'V',5F8.3,24X,2F8.3/
     *       ' ',74X,'R',5F8.3,32X,F8.3/
     *       ' ',74X,'T',5F8.3,32X,F8.3/
     *       ' ',74X,'S',5F8.3,24X,2F8.3)
C
      IF(IP.EQ.1.AND.HCFLAG.GE.2.AND.OUTFMT.EQ.6) WRITE(IOUREP,230)
     *   ICY,(EFEXH(IV),IV=1,3),EFEXH(9),(EFEXH(IV),IV=4,MAXVEH),VEXH,
     *   ICY,(EFEVAP(IV),IV=1,3),EFEVAP(9),EFEVAP(4),EFEVAP(8),VEVAP,
     *   ICY,(EFREFL(IV),IV=1,3),EFREFL(9),EFREFL(4),VLOSS,
     *   ICY,(EFRUNL(IV),IV=1,3),EFRUNL(9),EFRUNL(4),VRUNLS,
     *   ICY,(EFRSTL(IV),IV=1,3),EFRSTL(9),EFRSTL(4),EFRSTL(8),VRSTLS
C
  230 FORMAT('"Exhaust",',I4,10(',',F8.3),','/
     * '"Evap",',I4,',',5(F8.3,','),3('" ",'),2(F8.3,',')/
     * '"Refuel",',I4,',',5(F8.3,','),4('" ",'),F8.3,','/
     * '"Running",',I4,',',5(F8.3,','),4('" ",'),F8.3,','/
     * '"Resting",',I4,',',5(F8.3,','),3('" ",'),2(F8.3,','))
C
   10 CONTINUE
C
      IF(IDLFLG.EQ.1) RETURN
      DO 20 IP=IP1,IP2
      IPIDLE=IP + 3
      IF(OUTFMT.EQ.1) WRITE(IOUREP,240)
     * IPIDLE,(EFIDLE(IP,IV),IV=1,3),EFIDLE(IP,9),
     * (EFIDLE(IP,IV),IV=4,MAXVEH),VIDLE(IP)
  240 FORMAT(' ',74X,I1,10F8.3)
C
      IF(OUTFMT.EQ.6) WRITE(IOUREP,250)
     * IDLEF(IP),ICY,(EFIDLE(IP,IV),IV=1,3),EFIDLE(IP,9),
     * (EFIDLE(IP,IV),IV=4,MAXVEH),VIDLE(IP)
  250 FORMAT(A11,I4,',',10(F8.3,','))
C
   20 CONTINUE
C
      RETURN
      END
