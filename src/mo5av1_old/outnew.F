C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTNEW
c
c  LMOS Modification 04/08/93
c  Comment out OUTNEW which writes to M5OUT.DAT
c
cC
cC  OUTNEW writes out the new (user supplied replacement) basic emission
cC  rate parameters and data.
cC
cC  Called by OUTHD3 and OUTHD4.
cC
cC  Input on call:
cC
cC    common blocks:
cC    /BASEQ3/ KINKS
cC    /BASEQ4/ NEWPAR,BERNEW,NEWFIT,NEWCT
cC    /IOUCOM/ IOUREP
cC    /STRING/ NOYES
cC
cC  Notes:
cC
cC  None.
cC
cC
c      CHARACTER*4 NOYES,COMMA,PERIOD
c      CHARACTER*1 COLON
cC
c      COMMON /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW,KINKS
c      COMMON /BASEQ4/NEWPAR(5,100),BERNEW(3,100),NEWFIT(100),NEWCT,MAXCT
c      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
c      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
cC
c      WRITE(IOUREP,210)
c  210 FORMAT(
c     *'0',T16,'Emission Factor Modification Profile'/
c     *'+',T16,'____________________________________')
cC
c      IF(KINKS.EQ.1) WRITE(IOUREP,220)
c  220 FORMAT(
c     *'0','Equation  Reg  Veh  Pol  First MY  Last MY    Base     ',
c     *    'DR           Altered'/
c     *'+','________  ___  ___  ___  ________  _______    ____    _',
c     *    '___          _______')
cC
c      IF(KINKS.EQ.2) WRITE(IOUREP,230)
c  230 FORMAT(
c     *'0','Equation  Reg  Veh  Pol  First MY  Last MY    Base     ',
c     *    'DR   50K DR  Altered'/
c     *'+','________  ___  ___  ___  ________  _______    ____    _',
c     *    '___  ______  _______')
cC
c      DO 10 IFAC=1,NEWCT
c      INO=1
c      IF(NEWFIT(IFAC).EQ.1) INO=2
cC
c      IF(BERNEW(3,IFAC).LT.0.0) WRITE(IOUREP,240) IFAC,
c     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,2),NOYES(INO)
c  240 FORMAT(' ',T4,I3,T13,I1,T18,I1,T23,I1,T29,I4,T38,I4,T46,F6.2,
c     *           T54,F6.2,T72,A4)
cC
c      IF(BERNEW(3,IFAC).GE.0.0) WRITE(IOUREP,250) IFAC,
c     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,3),NOYES(INO)
c  250 FORMAT(' ',T4,I3,T13,I1,T18,I1,T23,I1,T29,I4,T38,I4,T46,F6.2,
c     *           T54,F6.2,T62,F6.2,T72,A4)
c   10 CONTINUE
cC
c
c  LMOS Modification 04/08/93
c  End of commented out code
c
      RETURN
      END
