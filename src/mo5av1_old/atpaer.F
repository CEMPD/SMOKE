C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE ATPAER(IATP,APGM,BFQT,CRATE,DDISTP)
C
C  ATPAER is designed to simplify the anti-tampering effectiveness rates
C  input by using a set of its parameters.
C
C  Called by GETATP and TAMPER
C
C  Calls ENFORC.
C
C  Input on call:
C
C    Parameter list: IATP,APGM,BFQT,CRATE,DDISTP
C    common blocks:
C    /ATPAR4/ EFF
C    /FLAGS2/ IMFLAG
C
C  Output on return:
C
C    Parameter list on return
C    common blocks:
C    /ATPAR3/ AIR,CAT,ECK,TNK,EGR
C    /ATPAR4/ CAN,PCV,CAP
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------------
C  NCRATE   R    New calculated compliance rate as a function
C                of DDISTP(ID).
C  NAPGM    I    New calculated inspection program type as a
C                function of DDISTP(ID).
C  NBFQT    I    New calculated frequency of program as a
C                function of DDISTP(ID).
C
C
C  Notes:
C
C  ATPAER was modified in Version 01 to reinitialize AIR/CAT/ECK/TNK/EGR/CAN/
C  PCV/CAP.
C  and again in MOBILE5.3 to handle IM240 catalyst removal program.
C
C
      REAL NCRATE
      INTEGER APGM,BFQT
      INTEGER ALHFLG,DDISTP(8)
C
      COMMON /ATPAR3/ AIR(2),CAT(2),ECK(3),TNK(2),EGR(2)
      COMMON /ATPAR4/ CAN(2),PCV(2),CAP(2),EFF(2,3,2,15)
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
C
C  Initialize
C
      DO 5 IAY=1,2
      AIR(IAY)=0.0
      CAT(IAY)=0.0
      ECK(IAY)=0.0
      TNK(IAY)=0.0
      EGR(IAY)=0.0
      CAN(IAY)=0.0
      PCV(IAY)=0.0
      CAP(IAY)=0.0
   5  CONTINUE
      ECK(3)=0.0
C
C  Loop through for Previous/Subsequent (IAY)
C
      DO 80 IAY=1,2
C
      ICOL=IAY+1
C
C  Assign Air Pump (ID=1)
C  (Annual: IROW=1, Biennial: IROW=2)
C
      IF(DDISTP(1).EQ.1) GOTO 10
      CALL NPARAM(DDISTP,1,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=1
      IF(NBFQT.EQ.2) IROW=2
C
      AIR(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign Catalyst (ID=2)
C
   10 IF(DDISTP(2).EQ.1) GOTO 20
      CALL NPARAM(DDISTP,2,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=3
C
      CAT(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign Catalyst/Inlet (Neck) (ID=3)
C  (No CAT inspection: IROW=4,  With I/M: IROW=5,
C  Without I/M program(s): IROW=6,
C  With Plumbtesmo Check: IROW=7; also assign NECK Previous Cat Removed)
C
   20 IF(DDISTP(3).EQ.1) GOTO 30
      CALL NPARAM(DDISTP,3,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IF(DDISTP(2).EQ.1) IROW=4
      IF(DDISTP(2).NE.1.AND.IMFLAG.NE.1.AND.DDISTP(4).EQ.1) IROW=5
      IF(DDISTP(2).NE.1.AND.IMFLAG.EQ.1.AND.DDISTP(4).EQ.1) IROW=6
      IF(DDISTP(2).NE.1.AND.DDISTP(4).NE.1) IROW=7
C
      NCKIAY=IAY+1
      ECK(NCKIAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
      ECK(1)=EFF(NAPGM,1,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign Tailpipe Lead Deposit (Tank) (ID=4)
C  (With CAT inspection: IROW=8, Without CAT inspection: IROW=9)
C
   30 IF(DDISTP(4).EQ.1) GOTO 40
      CALL NPARAM(DDISTP,4,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=8
      IF(DDISTP(2).NE.1) IROW=9
C
      TNK(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign EGR (ID=5)
C
   40 IF(DDISTP(5).EQ.1) GOTO 50
      CALL NPARAM(DDISTP,5,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=10
C
      EGR(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign Evap System (Canister) (ID=6)
C  (Annual: IROW=11, Biennial: IROW=12)
C  NO credit given if gas cap inspection not also included
C
   50 IF(DDISTP(6).EQ.1.OR.DDISTP(8).EQ.1) GOTO 60
      CALL NPARAM(DDISTP,6,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=11
      IF(BFQT.EQ.2) IROW=12
C
      CAN(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign PCV (ID=7)
C  (Annual: IROW=13, Biennial: IROW=14)
C
   60 IF(DDISTP(7).EQ.1) GOTO 70
      CALL NPARAM(DDISTP,7,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=13
      IF(NBFQT.EQ.2) IROW=14
C
      PCV(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
C  Assign Gas Cap (ID=8)
C
   70 IF(DDISTP(8).EQ.1) GOTO 80
      CALL NPARAM(DDISTP,8,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      IROW=15
C
      CAP(IAY)=EFF(NAPGM,ICOL,IATP,IROW)*ENFORC(NCRATE,2)
C
   80 CONTINUE
      RETURN
      END
