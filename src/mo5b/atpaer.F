C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE ATPAER(IATP,APGM,BFQT,CRATE,DDISTP)
C
C  ATPAER is designed to simplify the anti-tampering effectiveness
C  rates input by using a set of its parameters.
C
C  Called by GETATP and TAMPER
C
C  Calls ENFORC.
C
C  Input on call:
C
C    Parameter list: IATP,APGM,BFQT,CRATE,DDISTP
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /ATPAR2/ ATPPGM,ATPFQT
C    /ATPAR4/ EFF
C    /FLAGS3/ ATPFLG
C    /IMPAR1/ ICYIM
C    /IMPAR2/ ITEST
C    /IMPAR6/ IFREQ
C    /LOOKUP/ IVTAM
C    /MAXIMA/ MAXYRS
C    /MYCODE/ MY,IDX
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
C  ATPAER was modified in Version 01 to reinitialize
C  AIR/CAT/ECK/TNK/EGR/CAN/PCV/CAP.
C  and again in MOBILE5.3 to handle IM240 catalyst removal program.
C  12-April-1993 @ ARC-bsg Subtask 2-212 to correctly apply one year
C  biennial effects of IM240 and ATP
C  21-September-1993 @ ARC-bsg Request 2-282 fix the I/M effect on
C  tampering by adding more checking to determine I/M status
C
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'ATPAR3.I'
      INCLUDE 'ATPAR4.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYCODE.I'
C
      REAL NCRATE
      INTEGER APGM,BFQT
      INTEGER DDISTP(8)
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
C  This model year (MY) is an example model year in the range of
C  interest
C
      ICY=MY-IDX+MAXYRS
C
C  Loop through for Previous/Subsequent (IAY)
C
      DO 90 IAY=1,2
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
C  With Plumbtesmo Check: IROW=7; also assign NECK Previous CatRemoved)
C
   20 IF(DDISTP(3).EQ.1) GOTO 30
      CALL NPARAM(DDISTP,3,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
      LIM240=MYRIMT(MY,IVTAM,ICY)
      IF(DDISTP(2).EQ.1) IROW=4
      IF(DDISTP(2).NE.1.AND.LIM240.GT.0.AND.DDISTP(4).EQ.1) IROW=5
      IF(DDISTP(2).NE.1.AND.LIM240.EQ.0.AND.DDISTP(4).EQ.1) IROW=6
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
C  One year-old biennial ATP have only half credit
C
   80 IF(ATPFLG.EQ.2.AND.ICY.EQ.LAPSY+1) THEN
        IF(LVTFLG(IVTAM).EQ.2.AND.ATPFQT.EQ.2
     *  .AND.MY.GE.LAP1ST.AND.MY.LE.LAPLST) THEN
          AIR(IAY)=AIR(IAY)*0.5
          CAT(IAY)=CAT(IAY)*0.5
          ECK(IAY)=ECK(IAY)*0.5
          TNK(IAY)=TNK(IAY)*0.5
          EGR(IAY)=EGR(IAY)*0.5
          CAN(IAY)=CAN(IAY)*0.5
          PCV(IAY)=PCV(IAY)*0.5
          CAP(IAY)=CAP(IAY)*0.5
C
C  Skipping IM240 check prevents double decreasing of benefit
C  due to one year old biennial programs
C
          GOTO 90
        ENDIF
      ENDIF
C
C  When the IM240 is biennial the program is one year old
C  reduce benefits to one-half
C
      LIM240=MYRIMT(MY,IVTAM,ICY)
      IF(LIM240.NE.0) THEN
        IF(ITEST(LIM240).EQ.4.AND.ICY.EQ.ICYIM(LIM240)+1
     *  .AND.IFREQ(LIM240).EQ.2) THEN
          AIR(IAY)=AIR(IAY)*0.5
          CAT(IAY)=CAT(IAY)*0.5
          ECK(IAY)=ECK(IAY)*0.5
          TNK(IAY)=TNK(IAY)*0.5
          EGR(IAY)=EGR(IAY)*0.5
          CAN(IAY)=CAN(IAY)*0.5
          PCV(IAY)=PCV(IAY)*0.5
          CAP(IAY)=CAP(IAY)*0.5
        ENDIF
      ENDIF
   90 CONTINUE
      RETURN
      END
