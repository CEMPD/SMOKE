C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETPRG(INERR,*)
C
C  GETPRG gets and verifies purge check program parameters and effectiveness
C  rates.
C  This checks if a Purge Program is covered by an IM240.
C
C  Called by ONESEC.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS2/ IMFLAG
C    /IOUCOM/ IOUASK
C    /IMPAR1/ MODYR
C    /IMPAR2/ ILDT, ITEST
C    /YEARS4/ IY1960,IY2020
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  LVT      I    flag indicating whether or not a vehicle type is
C                covered by the Pressure Check program value: must be
C                either (1) no or (2) yes.
C  MY       I    see top of program.
C
C
C  Notes:
C
C  GETPRG was added in MOBILE5.3 so that it is called when PRGFLG=2
C  If the purge program is not covered by an IM240 the user gets
C  a warning about dynomometer.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER PRGPGM,PRGFQT,PRGFLG
      INTEGER ALHFLG
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG(4)
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C  Identify program parameters: start year, first year covered, last
C  year covered, vehicle types covered, inspection program type,
C  inspection frequency, and compliance rate.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ', 'Enter purge functional check start year,'/
     *   ' first & last years included,'/
     *   ' vehicle types covered, inspection type & frequency,'/
     *   ' and compliance rate (format =',
     *   ' 3(I2,1X),4I1,1X,2I1,1X,F4.0):')
      READ(IOUGEN,100,END=90,ERR=90) LPGSY,LPG1ST,LPGLST,LVPRG,
     *   PRGPGM,PRGFQT,CRPRG
  100 FORMAT(3(I2,1X),4I1,1X,2I1,1X,F4.0)
C
C  Verify parameters:
C
      CALL YRTEST(LPGSY,75,IY1960,IY2020,INERR)
      CALL YRTEST(LPG1ST,76,IY1960,IY2020,INERR)
      CALL YRTEST(LPGLST,77,IY1960,IY2020,INERR)
C
      DO 10 IVTAM=1,4
      LVT=LVPRG(IVTAM)
      IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR)
   10 CONTINUE
C
      IF(PRGPGM.LT.1.OR.PRGPGM.GT.2) CALL QUITER(0.,PRGPGM,103,INERR)
      IF(PRGFQT.LT.1.OR.PRGFQT.GT.2) CALL QUITER(0.,PRGFQT,104,INERR)
      IF(CRPRG.LT.0..OR.CRPRG.GT.100.) CALL QUITER(CRPRG,1,105,INERR)
C
      ICALL=0
      DO 60 IVTAM=1,4
      IF(LVPRG(IVTAM).EQ.1) GOTO 60
      IF(IMFLAG.EQ.1) THEN
      ICALL=1
      GOTO 60
      ENDIF
C
      DO 50 MY=LPG1ST, LPGLST
C
      IF(IMFLAG.EQ.3) GOTO 30
C
C  Assume the one I/M program case
C
      IF(ILDT(IVTAM,1).EQ.1)  ICALL=1
C
      IF(ITEST(1).NE.4)  ICALL=1
C
      IF(MY.LT.MODYR(1,1))  ICALL=1
C
      IF(MY.GT.MODYR(2,1))  ICALL=1
      GOTO 50
C
C  Assume the two I/M program case
C
   30 IF(ILDT(IVTAM,1).EQ.1) GOTO 40
C
      IF(ITEST(1).NE.4) GOTO 40
C
      IF(MY.LT.MODYR(1,1)) GOTO 40
C
C  Assume 1st I/M program is IM240
C
      IF(MY.LE.MODYR(2,1)) GOTO 50
   40 IF(ILDT(IVTAM,2).EQ.1)  ICALL=1
C
      IF(ITEST(2).NE.4)  ICALL=1
C
      IF(MY.LT.MODYR(1,2))  ICALL=1
C
      IF(MY.GT.MODYR(2,2))  ICALL=1
C
   50 CONTINUE
   60 CONTINUE
C
C  Warn users about dynomometer use with Purge Program.
C
      IF(ICALL.EQ.1) CALL QUITER(0.,0,114,INERR)
      RETURN
C
C  "Fatal" error in trying to read ATP program description
C
   90 CALL QUITER(0.,0,107,INERR)
C
   99 RETURN1
      END
