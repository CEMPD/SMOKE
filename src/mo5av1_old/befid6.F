C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID6(MY,IDX,VMTAGE)
C
C  BEFID6 returns the basic idle emission factor adjusted for I/M,
C  methane, tampering and oxy fuel, depending on the case.
C
C  Called by IDLCAL.
C
C  Calls BEFUN, IDLPTR, PCLEFT and FUEL.
C
C  Input on call:
C
C    parameter list: MY,IDX,VMTAGE
C    common blocks:
C    /IDLEQ1/ ZMLIDL,DRIDL
C    /MAXIMA/ MAXYRS
C    /OMTCOM/ OMTCF,OMTCFF
C    /REGION/ IREJN
C    /RVPEX2/ RVPCF
C    /SPEED6/ SCIADJ
C    /TAMID3/ OFFIDL
C
C  Output on return:
C
C    function: BEFID6
C
C  Local array subscripts:
C
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BER      R    untampered basic emission rate
C  CREDIM   R    I/M adjustment to idle emission factor
C
C
C  Notes:
C
C     BEFID6 is called by IDLCAL and calculates idle emissions
C     for NOx HDGV.  BEFID6 was created for MOBILE5.
C
C
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      BER=0.0
      BEFID6=0.0
C
C  Assign ICY and determine the I/M credit.
C
      ICY=MY-IDX+MAXYRS
      CREDIM=PCLEFT(MY,ICY,3,4)
C
      IF(MY.LT.1984) GOTO 10
C
C  NOx numbers are hard coded
C
      IF(MY.LE.1985) BEFID6=7.11
      IF(MY.GE.1986) BEFID6=2.35
      GOTO 20
C
C  In the pre-1984 model year case use the MOBILE4 revised idle EF
C  calculation algorithm.
C
   10 IGIDL=IDLPTR(MY,3,4)
      BEFID6=ZMLIDL(IGIDL,3,4,IREJN)+DRIDL(IGIDL,3,4,IREJN)*VMTAGE
C
C  Apply the I/M credit.
C
   20 BEFID6=BEFID6*CREDIM
C
C  Add in the Idle offset for Tampering
C
      BEFID6=BEFID6+OFFIDL(3,IDX,4)
C
C     IF(BEFID6.LT.0.0) BEFID6=0.0
C
C  Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID6=BEFID6*(OMTCF(IDX,3,4)/OMTCFF(IDX,3,4))
C
C  Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID6=BEFID6*SCIADJ(IDX,3,4)
C
C  Add RVP Correction Factor
C
      BEFID6=BEFID6*RVPCF(IDX,3,4)
C
C  Determine and apply the oxy fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
      BER=BER*CREDIM
      FUELCF=FUEL(ICY,3,4,MY,BER)
      BEFID6=BEFID6*FUELCF
      RETURN
      END
