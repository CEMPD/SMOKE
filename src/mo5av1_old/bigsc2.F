C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSC2(IP,IV,MY,IDX,ISP)
C
C  BIGSC2 computes SALHCF for all Motorcycles and for older (MY prior to
C  LB1STS(IV) ) LDGV/T, and for the High Speed ( .GT. 48 MPH ) CF.
C  It also computes the Idle Speed Correction Factor Adjustment for this group.
C
C  Called by BIGSAL.
C
C  Calls ISPPTR and EXP.
C
C  Input on call:
C
C    parameter list: IP,IV,MY,IDX,ISP
C    common blocks:
C    /SCENE1/ SPD
C    /SPEED1/ SK1
C    /SPEED3/ SADJ,SCUNA1
C    /SPEED4/ IVA,SCUNA0
C    /SPEED8/ SIXT5
C
C  Output on return:
C
C    common blocks:
C    /SPEED3/ SCADJ1,ISPGRP
C    /SPEED6/ SALHCF,SCIADJ
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BB       R    slope coefficient "BB" of high speed scf equ.
C  IVA2     I    used to point to the LDGV for motorcycles in SIXT5.
C  MYR      I    model year pointer for SIXT5 array
C  SCF1     R    old scf equation computed using opmode adjusted speed value
C  SP       R    Vehicle speed for high speed equation
C
C  Notes:
C
C  BIGSC2 added for MOBILE4.1 v4.
C  BIGSC2 was modified for MOBILE4.1.
C
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,3),SK4C(2,14,3,3)
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /SPEED8/ SIXT5(11,3,3)
C
      IF(IV.EQ.8) ISP=ISPPTR(MY,IP,0,2)
C
      IF(IDX.GT.1.AND.ISPGRP(IP,IVA).EQ.ISP) GOTO 20
C
C  Use MOBILE3's 5th degree polynomial equation.
      SCF1=SK1(6,ISP,IP)
      DO 10 ICOEF5=1,5
      SCF1=SCF1*SADJ+SK1(6-ICOEF5,ISP,IP)
   10 CONTINUE
C
      IF(IP.NE.3) SCF1=EXP(SCF1)
      SCADJ1(IP,IVA)=SCF1
      ISPGRP(IP,IVA)=ISP
C
 20   CONTINUE
C
C  Original SALHCF now used for low/mid speed CF only
      SALHCF(IDX,IP,IV)=SALHCF(IDX,IP,IV)*
     1                 (SCUNA1(ISP,IP,IV)/SCADJ1(IP,IVA))
C
C  Compute the Idle Speed Correction Factor Adjustment.
C  Use the 16mph SCUNA0 value.
C
      SCIADJ(IDX,IP,IV)=SCUNA0(ISP,IP,IV)/SCADJ1(IP,IVA)
C
C  Apply high speed correction when necessary.  For HC/CO and speeds less
C  than or equal to 55 mph, use 48 mph value.  SCF is constant in this case.
C
      IF(IP.EQ.3.AND.SPD(IV).GT.48.0.OR.IP.NE.3.AND.SPD(IV).GT.55.0)
     *  THEN
C
C  New equation form for High Speed(.GT.48MPH) is:
C  Speed Correction Factor * ( 1.0 + ( BB * speed-(set spd ))
C  SP is either 48 for NOx or 55 for HC/CO
C  Use 1981 model year LDGV SIXT5 value for motorcyles.
C
        SP=55.0
        IF(IP.EQ.3) SP=48.0
C
        MYR=MY-1981+1
        IF(MY.LT.1981) MYR=1
        IF(MY.GT.1991) MYR=11
        IVA2=IVA
        IF(IV.EQ.8) THEN
          IVA2=1
          MYR=1
        ENDIF
C
        BB=(SIXT5(MYR,IP,IVA2)-SALHCF(IDX,IP,IV))/(65.0-SP)
C
        SALHCF(IDX,IP,IV)=SALHCF(IDX,IP,IV)*(1.0+(BB*(SPD(IV)-SP) ))
C
      ENDIF
C
      RETURN
      END
