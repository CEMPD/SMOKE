C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE USDPTR
C
C  USDPTR finds the pointers into the user supplied data name array.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ TAMFLG
C    /FLAGS2/ MYMRFG,NEWFLG
C
C  Output on return:
C
C    common block: /USDATA/ NUSD,IUSD
C
C  Notes:
C
C  USDPTR was modified in Version 01 to reinitialize the "name" parameters.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG
      CHARACTER*8 USNAME
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /USDATA/ USNAME(4,4),NUSD,IUSD(4)
C
C  Initially set the "name" parameters to the "no user supplied one-time data"
C  setting.
C
      DO 10 NUSD=1,4
      IUSD(NUSD)=0
   10 CONTINUE
      NUSD=0
C
C  Setup scenario header user supplied data pointer vector.
C
      IF(TAMFLG.LE.1.AND.NEWFLG.EQ.1.AND.MYMRFG.EQ.1) GOTO 99
      IF(TAMFLG.LE.1.AND.NEWFLG.EQ.3.AND.MYMRFG.EQ.1) GOTO 99
      IF(TAMFLG.LE.1.AND.NEWFLG.EQ.5.AND.MYMRFG.EQ.1) GOTO 99
      IF(TAMFLG.LE.1) GOTO 20
      NUSD=1
      IUSD(1)=1
   20 IF(NEWFLG.EQ.1) GOTO 30
      IF(NEWFLG.EQ.3) GOTO 30
      IF(NEWFLG.EQ.5) GOTO 30
      NUSD=NUSD+1
      IUSD(NUSD)=2
   30 IF(MYMRFG.EQ.1) GOTO 99
      IF(MYMRFG.EQ.4) GOTO 40
      NUSD=NUSD+1
      IF(MYMRFG.EQ.2) IUSD(NUSD)=3
      IF(MYMRFG.EQ.3) IUSD(NUSD)=4
      GOTO 99
   40 NUSD=NUSD+1
      IUSD(NUSD)=3
      NUSD=NUSD+1
      IUSD(NUSD)=4
C
   99 RETURN
      END
