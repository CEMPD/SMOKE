      SUBROUTINE FINDPB(ICY,ACTUAL_MY,IDX) 
C 
C  FINDPB looks up the misfueling categories tampering rates for the given 
C  IVTAM x IDX case, when ICY exceeds LPOD, and assigns them to BTR. 
C 
C  Called by TAMPER. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I 
C  01 Oct 1998 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    common blocks: 
C    /MAXIMA/ MAXYRS 
C    /LOOKUP/ IVTAM 
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
      IMPLICIT NONE 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'TAMEQ4.I' 
      INCLUDE 'TAMPB1.I' 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER LEADX 
      INTEGER IAY 
 
C 
C  The correct misfueling rates are those cumulated by the vehicle through 
C  LOPD.  Afterward, no further misfueling related tampering occurs.  If 
C  the vehicle was produced after LPOD, of course, misfueling tampering is 0. 
C 
      IF(ACTUAL_MY.GT.LPOD) GOTO 20 
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
