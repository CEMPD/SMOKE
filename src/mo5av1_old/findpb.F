C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE FINDPB(ICY)
C
C  FINDPB looks up the misfueling categories tampering rates for the given
C  IVTAM x IDX case, when ICY exceeds LPOD, and assigns them to BTR.
C
C  Called by TAMPER.
C
C  Input on call:
C
C    common blocks:
C    /MAXIMA/ MAXYRS
C    /LOOKUP/ IVTAM
C    /MYCODE/ MY,IDX
C    /TAMPB1/ LPOD,PBBTR
C
C  Output on return:
C
C    common blocks:
C    /TAMEQ4/ BTR (misfueling categories only)
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  LEADX    I    ?DX order pointer to the LPOD misfueling rates corresponding to
C                ICY's IVTAM x IDX case passed into the routine
C
C  Notes:
C
C  PBSAVE was added for MOBILE4.1v6.
C
C  PBSAVE is identical to TAMPER, except that LPOD is used instead of ICY and
C  the save to PBBTR replaces the calls to EEFGRP, EMIRAT, BAGEME, and EVPEME.
C
C
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMEQ4/ BTR(9,2)
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
C
C  The correct misfueling rates are those cumulated by the vehicle through
C  LOPD.  Afterward, no further misfueling related tampering occurs.  If
C  the vehicle was produced after LPOD, of course, misfueling tampering is 0.
C
      IF(MY.GT.LPOD) GOTO 20
      LEADX=IDX+(ICY-LPOD)
C
C  Retrieve the corresponding set of LPOD misfueling categories tampering rates.
C
      DO 10 IAY=1,2
      BTR(3,IAY)=PBBTR(1,IAY,LEADX,IVTAM)
      BTR(4,IAY)=PBBTR(2,IAY,LEADX,IVTAM)
      BTR(9,IAY)=PBBTR(3,IAY,LEADX,IVTAM)
  10  CONTINUE
      RETURN
C
  20  DO 30 IAY=1,2
      BTR(3,IAY)=0.0
      BTR(4,IAY)=0.0
      BTR(9,IAY)=0.0
  30  CONTINUE
C
      RETURN
      END
