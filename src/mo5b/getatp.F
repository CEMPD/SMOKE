C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETATP(INERR,*)
C
C  GETATP gets and verifies anti-tampering program parameters and effectiveness
C  rates.
C
C  Called by ONESEC.
C
C  Calls ATPAER and CALAER.
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
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /EGSCAL/ AER
C    /EVPDAT/ EVPAER
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COLSUM   R    column sum values from AER: each column's sum must be
C                within the closed interval from 0.0 to 1.0.
C  EVP      R    evaporative A-T effectiveness rate: each rate must be
C                within the closed inetrval from 0.0 to 1.0.
C  LVT      I    flag indicating whether or not a vehicle type is
C                covered by the A-T program: value must be either 1 (no)
C                or 2 (yes).
C  RATE     R    A-T effectiveness rate against exhaust emissions: each
C                rate must be within the closed interval from 0.0 to 1.0
C
C  Notes:
C
C  GETATP was modified in MOBILE5.
C  ATP vs I/M error checking is not done it GETATP.  This is now done in
C  three places, ONESEC, QUITER and CHCKIM. (In increasing order of
C  responsibility and importance).
C
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'EGSCAL.I'
      INCLUDE 'EVPDAT.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'YEARS4.I'
C
C  Identify program parameters: start year, first year covered, last
C  year covered, vehicle types covered, inspection program type,
C  inspection frequency, compliance rate, and disablement inspections
C  conducted.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ', 'Enter ATP start year, first & last years included,'/
     *   ' vehicle types covered, inspection type & frequency,'/
     *   ' compliance rate, and inspections conducted (format ='/
     *   ' 3(I2,1X),4I1,1X,2I1,1X,F4.0,1X,8I1):')
      READ(IOUGEN,100,END=90,ERR=90) LAPSY,LAP1ST,LAPLST,LVTFLG,
     *   ATPPGM,ATPFQT,CRATP,DISTYP
  100 FORMAT(3(I2,1X),4I1,1X,2I1,1X,F4.0,1X,8I1)
C
C  Verify parameters:
C
      CALL YRTEST(LAPSY,75,IY1960,IYEND,INERR)
      CALL YRTEST(LAP1ST,76,IY1960,IYEND,INERR)
      CALL YRTEST(LAPLST,77,IY1960,IYEND,INERR)
C
      DO 10 IVTAM=1,4
      LVT=LVTFLG(IVTAM)
      IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR)
   10 CONTINUE
C
      IF(ATPPGM.LT.1.OR.ATPPGM.GT.2) CALL QUITER(0.,ATPPGM,103,INERR)
      IF(ATPFQT.LT.1.OR.ATPFQT.GT.2) CALL QUITER(0.,ATPFQT,104,INERR)
      IF(CRATP.LT.0..OR.CRATP.GT.100.) CALL QUITER(CRATP,1,105,INERR)
C
      DO 20 I=1,8
      IF(DISTYP(I).LT.1.OR.DISTYP(I).GT.2)
     *   CALL QUITER(0.,DISTYP(I),106,INERR)
   20 CONTINUE
      IF(DISTYP(6).EQ.2.AND.DISTYP(8).NE.2)
     *   CALL QUITER(0.,0,121,INERR)
C
C  Assign base ATP effectiveness rates.
C
      CALL ATPAER(1,ATPPGM,ATPFQT,CRATP,DISTYP)
C
C  Calculate ATP effectiveness rates.
C
      CALL CALAER
C
C  EVPAER and AER checks.
C
      DO 40 IAY=1,2
      DO 30 IHG=1,3
      EVP=EVPAER(IHG,IAY)
      IF(EVP.LT.0.0.OR.EVP.GT.1.0) CALL QUITER(EVP,0,79,INERR)
   30 CONTINUE
   40 CONTINUE
C
C  Note: ICOL / IROW here refers to the input file layout of the arrays.
C  In subroutine ATPAER the labels are switched, to indicate the actual
C  orientation of the matrices in the effectiveness computation process.
C
      DO 80 IAY=1,2
      DO 70 IPG=1,2
      DO 60 ICOL=1,11
      COLSUM=0.0
      DO 50 IROW=1,11
      RATE=AER(ICOL,IROW,IPG,IAY)
      IF(RATE.LT.0.0.OR.RATE.GT.1.0) CALL QUITER(RATE,0,80,INERR)
      COLSUM=COLSUM+RATE
   50 CONTINUE
      IF(COLSUM.LT.0.0.OR.COLSUM.GT.1.0) CALL QUITER(COLSUM,0,81,INERR)
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
      RETURN
C
C  "Fatal" error in trying to read ATP program description
C
   90 CALL QUITER(0.,0,107,INERR)
C
   99 RETURN1
      END
