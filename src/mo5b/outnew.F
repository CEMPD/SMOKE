C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTNEW
c
cmh  LMOS Modification from Mobil5a
cmh  Comment out OUTNEW which writes to M5OUT.DAT
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
c      INCLUDE 'BASEQ3.I'
c      INCLUDE 'BASEQ4.I'
c      INCLUDE 'IOUCOM.I'
c      INCLUDE 'STRING.I'
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
cmh  LMOS Modification from Mobile5a
cmh  End of commented out code
c
      RETURN
      END
