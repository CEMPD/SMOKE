C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETVMT(INERR,*)
C
C  GETVMT gets the fleet vehicle miles travelled (vmt) mix by vehicle type.
C
C  Called by ONESEC and PARSEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C    /MAXIMA/ MAXVEH
C
C  Output on return:
C
C    parameter list: INERR
C    common block: /VMXCOM/ VMTMIX
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  VMTCHK   R    vmt mix check: sum of the vmt fractions = 1.0 +/- 1.E-6
C                (i.e., a normalized vmt distribution has been entered)
C
C  Notes:
C
C  None.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter VMT split:')
      READ(IOUGEN,100,END=99) VMTMIX
  100 FORMAT(8F4.3)
C
      VMTCHK=0.0
      DO 10 IV=1,MAXVEH
      VMTCHK=VMTMIX(IV)+VMTCHK
      IF(VMTMIX(IV).LT.0.0.OR.VMTMIX(IV).GT.1.)
     *      CALL QUITER(VMTMIX(IV),0,17,INERR)
   10 CONTINUE
      IF(VMTCHK+1.E-6.LT.1..OR.VMTCHK-1.E-6.GT.1.)
     *  CALL QUITER(VMTCHK,0,18,INERR)
C
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,127,INERR)
C
      RETURN1
      END
