C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETPRS(INERR,*)
C
C  GETPRS gets and verifies pressure check program parameters and effectiveness
C  rates.
C
C
C  Called by ONESEC.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUASK
C    /YEARS4/ IY1960,IYEND
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  LVT      I    flag indicating whether or not a vehicle type is
C                covered by the Pressure Check program: value must be
C                either 1 (no) or 2 (yes).
C
C
C  Notes:
C
C  it was added in MOBILE5.3 so that GETPRS is called when PRSFLG=2
C
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C       Change PRSPGM range checking in the IF statement
C
      INCLUDE 'PRSCH1.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'YEARS4.I'
C
C  Identify program parameters: start year, first year covered, last
C  year covered, vehicle types covered, inspection program type,
C  inspection frequency, and compliance rate.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ', 'Enter pressure functional check start year,'/
     *   ' first & last years included,'/
     *   ' vehicle types covered, inspection type & frequency,'/
     *   ' and compliance rate (format =',
     *   ' 3(I2,1X),4I1,1X,2I1,1X,F4.0):')
      READ(IOUGEN,100,END=90,ERR=90) LPRSY,LPR1ST,LPRLST,LVPRS,
     *   PRSPGM,PRSFQT,CRPRS
  100 FORMAT(3(I2,1X),4I1,1X,2I1,1X,F4.0)
C
C  Verify parameters:
C
      CALL YRTEST(LPRSY,75,IY1960,IYEND,INERR)
      CALL YRTEST(LPR1ST,76,IY1960,IYEND,INERR)
      CALL YRTEST(LPRLST,77,IY1960,IYEND,INERR)
C
      DO 10 IVTAM=1,4
      LVT=LVPRS(IVTAM)
      IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR)
   10 CONTINUE
C
      IF(PRSPGM.LT.1.OR.PRSPGM.GT.3) CALL QUITER(0.,PRSPGM,103,INERR)
      IF(PRSFQT.LT.1.OR.PRSFQT.GT.2) CALL QUITER(0.,PRSFQT,104,INERR)
      IF(CRPRS.LT.0..OR.CRPRS.GT.100.) CALL QUITER(CRPRS,1,105,INERR)
C
      RETURN
C
C  "Fatal" error in trying to read ATP program description
C
   90 CALL QUITER(0.,0,107,INERR)
C
   99 RETURN1
      END
