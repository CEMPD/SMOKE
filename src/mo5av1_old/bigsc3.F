C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSC3(IP,MY,IDX,ISP)
C
C  BIGSC3 calculates the SALHCF for the Low, Mid and High
C  Speeds using the MOBILE4 algorithm.
C
C  Called by BIGSAL.
C
C  Calls ISPPTR and EXP
C
C  Input on call:
C
C    parameter list: IP,MY,IDX,ISP
C    common blocks:
C    /SCENE1/ SPD
C    /SPEED1/ SK2C,SK4C
C    /SPEED3/ SADJ
C    /SPEED4/ IVA
C    /SPEED5/ MYHSPC
C    /SPEED8/ SIXT5
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCIADJ
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  A        R    coefficient "A" of new scf equation: A/X+B or EXP(A+BX+CX**2)
C  A16      R    16.0 mph based coefficient "A" of new scf equation: A/X+B
C  A2       R    SADJ based coefficient "A" of new scf equation: A/X+B
C  B        R    coefficient "B" of new scf equation: A/X+B or EXP(A+BX+CX**2)
C  B16      R    16.0 mph based coefficient "B" of new scf equation: A/X+B
C  B2       R    SADJ based coefficient "B" of new scf equation: A/X+B
C  BB       R    slope coefficient "BB" of high speed scf equ.
C  MYR      I    model year pointer for SIXT5 array
C  S        R    Vehicle speed, held between 2.5 and 48 mph, used in equations.
C  SP       R    Vehicle speed for high speed equation
C  SCADJ2   R    new scf equation computed using opmode adjusted speed value
C  SCUNA2   R    new scf equation computed using unadjusted user inputted speed
C  SCUNAi   R    idle scf adjustment computed using 16.0 mph speed
C
C  Notes:
C
C  BIGSC3 added for MOBILE4.1 v5.
C  BIGSC3 was modified in MOBILE4.1.
C
C
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,3),SK4C(2,14,3,3)
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /SPEED8/ SIXT5(11,3,3)
C
C  Set the speed.
C
      S=SPD(IVA)
      IF(S.LT.2.5) S=2.5
      IF(S.GT.48.0) S=48.0
C
C  For NOx MY 1977 - 1979, IV 1, 2, 3 ( LDGV,
C  LDGT1, and LDGT2 ); A, B, and C are constants.
C  The equation for NOx is EXP(A+B*X+C*X**2) / EXP(A+B*Y+C*Y**2),
C  where X & Y are defined as in the HC & CO case.
C  Hard Code the 1977 - 1979 NOx Coeff's.
C
      IF(MY.LT.1980.AND.IP.EQ.3) THEN
      A = 0.34665
      B =-0.026082
      C = 4.2835E-04
C
      SCUNA2=EXP(A+B*S+C*S**2)
      SCADJ2=EXP(A+B*SADJ+C*SADJ**2)
      SCUNAI=EXP(A+B*16.0+C*16.0**2)
      GOTO 10
      ENDIF
C
C  New equation form is (A/X+B) / (A/Y+B), where X = inputted speed
C  and Y = operating mode adjusted FTP speed.
C
      ISP=ISPPTR(MY,IP,IVA,3)
C
C  Adjust ISP to work with array SK2C and SK4C
C
      IF(IVA.GT.1) ISP=ISP+1
      IF(MY.LT.1980) ISP=1
      IF(ISP.GT.14) ISP=14
C
C  Set low speed coefficients (<=19.6 mph).
C
      IF(S.LE.19.6) THEN
      A=SK2C(1,ISP,IP,IVA)
      B=SK2C(2,ISP,IP,IVA)
      ENDIF
C
C     speed adjust
C
      IF(SADJ.LE.19.6) THEN
      A2=SK2C(1,ISP,IP,IVA)
      B2=SK2C(2,ISP,IP,IVA)
      ENDIF
C
      A16=SK2C(1,ISP,IP,IVA)
      B16=SK2C(2,ISP,IP,IVA)
C
C  Set mid to high speed coefficients (>19.6 mph).
C
      IF(S.GT.19.6) THEN
      A=SK4C(1,ISP,IP,IVA)
      B=SK4C(2,ISP,IP,IVA)
      ENDIF
C
C     speed adjust
C
      IF(SADJ.GT.19.6) THEN
      A2=SK4C(1,ISP,IP,IVA)
      B2=SK4C(2,ISP,IP,IVA)
      ENDIF
C
C  Determine the SCFs and apply.
C
      SCUNA2=(A/S)+B
      SCADJ2=(A2/SADJ)+B2
      SCUNAI=(A16/16.0)+B16
C
 10   SALHCF(IDX,IP,IVA)=SALHCF(IDX,IP,IVA)*(SCUNA2/SCADJ2)
C
C  Determine the Idle Speed Correction Factor Adjustment.
C
      SCIADJ(IDX,IP,IVA)=SCUNAI/SCADJ2
C
C  Apply high speed correction when necessary.  For HC/CO and speeds less
C  than or equal to 55 mph, use 48 mph value.  SCF is constant in this case.
C
      IF(IP.EQ.3.AND.SPD(IVA).GT.48.0.OR.IP.NE.3.AND.SPD(IVA).GT.55.0)
     *  THEN
C
C  New equation form for High Speed(.GT.48MPH) is:
C  Speed Correction Factor * ( 1.0 + ( BB * speed-(set spd ))
C  SP is either 48 for NOx or 55 for HC/CO
C
        SP=55.0
        IF(IP.EQ.3) SP=48.0
C
        MYR=MY-1981+1
        IF(MY.LT.1981) MYR=1
        IF(MY.GT.1991) MYR=11
        BB=(SIXT5(MYR,IP,IVA)-SALHCF(IDX,IP,IVA))/(65.0-SP)
C
        SALHCF(IDX,IP,IVA)=SALHCF(IDX,IP,IVA)*(1.0+(BB*(SPD(IVA)-SP) ))
C
        ENDIF
C
      RETURN
      END
