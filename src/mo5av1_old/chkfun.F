C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION CHKFUN(MYR,IVTAM,ICY)
C
C  CHKFUN discovers the status of the Pressure/Purge Check
C  for a given MYR, Vehicle, and CY for functionality.
C
C  Called by FAIL.
C
C  Input on call:
C
C    Parameter list: MYR,IVTAM,ICY
C    common blocks:
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG
C    /PRGCH2/ PRGFLG
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS
C    /PRSCH2/ PRSFLG
C
C  Output on return:
C
C    Function: CHKFUN
C    1=none
C    2=pressure only
C    3=purge only
C    4=both
C
C
C  Notes:
C
C  This function was developed for FAIL re-write ordered as part of
C  a March 31, 1992 assignment to allow multiple functionality
C  checking.  It was added to job-out some of the "thought" in
C  the new MOBILE5.3 FAIL.  This was added in MOBILE5.3.
C
C
      INTEGER PRGPGM,PRGFQT,PRSPGM,PRSFQT
      INTEGER PRGFLG,PRSFLG
C
      COMMON /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG(4)
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS(4)
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
C
      CHKFUN = 1
C
C  Here is one possible exclusion if the vehicle is not available.
C
      IF(IVTAM.GT.4) RETURN
C
      IF(PRSFLG.EQ.1) GOTO 10
C
C  The following is kept as four separate IFs for clarity only
C
      IF(ICY.LT.LPRSY) GOTO 10
      IF(MYR.LT.LPR1ST) GOTO 10
      IF(MYR.GT.LPRLST) GOTO 10
      IF(LVPRS(IVTAM).EQ.1) GOTO 10
C
C  At this point Pressure is yes.
C
      CHKFUN = 2
C
C  Jumping here indicates Pressure is no.
C
   10 IF(PRGFLG.EQ.1) RETURN
C
C  The following is kept as four separate IFs for clarity only.
C
      IF(ICY.LT.LPGSY) RETURN
      IF(MYR.LT.LPG1ST) RETURN
      IF(MYR.GT.LPGLST) RETURN
      IF(LVPRG(IVTAM).EQ.1) RETURN
C
C  At this point Purge is yes.
C
      IF(CHKFUN.EQ.2) THEN
         CHKFUN = 4
         RETURN
      ENDIF
C
C  At this point Purge is the only function in effect.
C
      CHKFUN = 3
      RETURN
      END
