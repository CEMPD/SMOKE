C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETBYM(INERR,*)
C
C  GETBYM gets the By Model Year Tables selection flags.
C
C  Called by ONESEC.
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
C    common block: /BYMYC1/: IVPICK,IMPICK
C
C  Notes:
C
C  GETBYM was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT('&','Enter By Model Year Table picks:')
      READ(IOUGEN,100,ERR=99,END=99) IVPICK,IMPICK
  100 FORMAT(8I1,I2)
C
      DO 10 IV=1,MAXVEH
      IF(IVPICK(IV).LT.1.OR.IVPICK(IV).GT.2)
     *      CALL QUITER(0.0,IVPICK(IV),117,INERR)
   10 CONTINUE
      IF(IMPICK.LT.1.OR.IMPICK.GT.2)
     *      CALL QUITER(0.0,IMPICK,117,INERR)
C
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,134,INERR)
C
      RETURN1
      END
