      SUBROUTINE OUTNEW 
C 
C  OUTNEW writes out the new (user supplied replacement) basic emission 
C  rate parameters and data. 
C 
C  Called by OUTHD3. 
C 
C  Changes: (Last change first) 
C 
C 11 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /BASEQ3/ KINKS 
C    /BASEQ4/ NEWOML,NEWPAR,BERNEW,NEWFIT,NEWCT 
C    /IOUCOM/ IOUREP 
C    /STRING/ NOYES 
C 
C  Notes: 
C 
C  11 Jun 97 DynTel@MLA 2-622  Altered code to print outputs for the new 
C            vehicle types that were defined in Mobile6. 
C 
C  13 Apr 98 @DynTel-RJD 2-654 Added field to the output to indicate if 
C            the user input basic exhaust emissions apply to vehicles 
C            in starting or running mode. 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ3.I' 
      INCLUDE 'BASEQ4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'STRING.I' 
      INTEGER IFAC, INO, I 
C       
      WRITE(IOUREP,210) 
  210 FORMAT( 
     *'0',T16,'Emission Factor Modification Profile'/ 
     *'+',T16,'____________________________________') 
C 
      IF(KINKS.EQ.1) WRITE(IOUREP,220) 
  220 FORMAT( 
     *'0','Equation  Mode Reg  Veh  Pol  First MY  Last MY    Base ', 
     *    '    DR           Altered'/ 
     *'+','________  ____ ___  ___  ___  ________  _______    ____ ', 
     *    '   ____          _______') 
C 
      IF(KINKS.EQ.2) WRITE(IOUREP,230) 
  230 FORMAT( 
     *'0','Equation  Mode Reg  Veh  Pol  First MY  Last MY    Base ', 
     *    '    DR   50K DR  Altered'/ 
     *'+','________  ____ ___  ___  ___  ________  _______    ____ ', 
     *    '   ____  ______  _______') 
C 
      DO 10 IFAC=1,NEWCT 
      INO=1 
      IF(NEWFIT(IFAC).EQ.1) INO=2 
C 
      IF(BERNEW(3,IFAC).LT.0.0) WRITE(IOUREP,240) IFAC, NEWOML(IFAC), 
     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,2),NOYES(INO) 
  240 FORMAT(' ',T4,I3,T13,A1,T17,I2,T23,I2,T28,I1,T34,I4,T43,I4,T51, 
     *           F6.2,T59,F6.2,T77,A4) 
C 
      IF(BERNEW(3,IFAC).GE.0.0) WRITE(IOUREP,250) IFAC, NEWOML(IFAC), 
     *   (NEWPAR(I,IFAC),I=1,5),(BERNEW(I,IFAC),I=1,3),NOYES(INO) 
  250 FORMAT(' ',T4,I3,T13,A1,T17,I2,T23,I1,T28,I1,T34,I4,T43,I4,T51, 
     *           F6.2,T59,F6.2,T67,F6.2,T77,A4) 
   10 CONTINUE 
C 
      RETURN 
      END 
