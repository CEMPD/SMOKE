C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /SPEED4/ IVA
C    /SPEED6/ SALHCF,SCFIDL,IDLYES
C    /SPEED8/ SIXT5
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCFIDL
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
C  BIGSC2 was modified for MOBILE5a in house version for new idle algorithm.
C  15 July 1993 @ ARC-rwk Task 2-277
C
C
      INCLUDE 'SCENE1.I'
      INCLUDE 'SPEED1.I'
      INCLUDE 'SPEED3.I'
      INCLUDE 'SPEED4.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'SPEED8.I'
C
      IF(IV.EQ.8) ISP=ISPPTR(MY,IP,0,2)
C
C  Use MOBILE3's 5th degree polynomial equation.
C
      SCF1=SK1(6,ISP,IP)
      DO 10 ICOEF5=1,5
      SCF1=SCF1*SADJ+SK1(6-ICOEF5,ISP,IP)
   10 CONTINUE
C
      IF(IP.NE.3) SCF1=EXP(SCF1)
C
C  Original SALHCF now used for low/mid speed CF only
C
      IF(IDLYES) THEN
        SCFIDL=SCFIDL*(SCUNA1(ISP,IP,IV)/SCF1)
      ELSE
        SALHCF(IDX,IP,IV)=SALHCF(IDX,IP,IV)*(SCUNA1(ISP,IP,IV)/SCF1)
      ENDIF
C
C  Apply high speed correction when necessary.  For HC/CO and speeds less
C  than or equal to 55 mph, use 48 mph value.  SCF is constant in this case.
C
      IF((IP.EQ.3.AND.SPD(IV).GT.48.0.OR.IP.NE.3.AND.SPD(IV).GT.55.0)
     *  .AND..NOT.IDLYES) THEN
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
