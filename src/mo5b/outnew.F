C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTNEW
C
C  OUTNEW writes out the new (user supplied replacement) basic emission
C  rate parameters and data.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /BASEQ3/ KINKS
C    /BASEQ4/ NEWPAR,BERNEW,NEWFIT,NEWCT
C    /IOUCOM/ IOUREP
C    /STRING/ NOYES
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'BASEQ3.I'
      INCLUDE 'BASEQ4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'STRING.I'
C
      WRITE(IOUREP,210)
  210 FORMAT(
     *'0',T16,'Emission Factor Modification Profile'/
     *'+',T16,'____________________________________')
C
      IF(KINKS.EQ.1) WRITE(IOUREP,220)
  220 FORMAT(
     *'0','Equation  Reg  Veh  Pol  First MY  Last MY    Base     ',
     *    'DR           Altered'/
     *'+','________  ___  ___  ___  ________  _______    ____    _',
     *    '___          _______')
C
      IF(KINKS.EQ.2) WRITE(IOUREP,230)
  230 FORMAT(
     *'0','Equation  Reg  Veh  Pol  First MY  Last MY    Base     ',
     *    'DR   50K DR  Altered'/
     *'+','________  ___  ___  ___  ________  _______    ____    _',
     *    '___  ______  _______')
C
      DO 10 IFAC=1,NEWCT
      INO=1
      IF(NEWFIT(IFAC).EQ.1) INO=2
C
      IF(BERNEW(3,IFAC).LT.0.0) WRITE(IOUREP,240) IFAC,
     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,2),NOYES(INO)
  240 FORMAT(' ',T4,I3,T13,I1,T18,I1,T23,I1,T29,I4,T38,I4,T46,F6.2,
     *           T54,F6.2,T72,A4)
C
      IF(BERNEW(3,IFAC).GE.0.0) WRITE(IOUREP,250) IFAC,
     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,3),NOYES(INO)
  250 FORMAT(' ',T4,I3,T13,I1,T18,I1,T23,I1,T29,I4,T38,I4,T46,F6.2,
     *           T54,F6.2,T62,F6.2,T72,A4)
   10 CONTINUE
C
      RETURN
      END
