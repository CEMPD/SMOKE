C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID5(MY,IDX,VMTAGE)
C
C  BEFID5 returns the basic idle emission factor adjusted for I/M,
C  methane, tampering and oxy fuel.
C
C  Called by IDLCAL
C
C  Calls IDLPTR, IERPTR, PCLEFT and FUEL.
C
C  Input on call:
C
C    parameter list: MY,IDX,VMTAGE
C    common blocks:
C    /BASEQ8/ IHDVMY,HDVZML,HDVBDR
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
C    function: BEFID5
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  BER      R    untampered basic emission rate
C  CREDIM   R    I/M adjustment to idle emission factor
C  IGIDL    I    pointer
C  IDLGER   I    pointer
C
C
C  Notes:
C
C     BEFID5 is called by IDLCAL and calculates idle emissions
C     for CO HDGV.  BEFID5 was created for MOBILE5.
C
C
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASEQ8/ IHDVMY(23,2),HDVZML(23,2),HDVBDR(23,2)
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      BER=0.0
      BEFID5=0.0
C
C  Assign ICY and determine the I/M credit.
C
      ICY=MY-IDX+MAXYRS
      CREDIM=PCLEFT(MY,ICY,2,4)
C
      IF(MY.LT.1984) GOTO 10
C
C  Obtain an untampered Basic Emission Rate
C
      CALL IERPTR(MY,2,4,IDLGER)
C
C  HDV CO Idle correction see Block Data 04
C
      BER=HDVZML(IDLGER,IREJN)+HDVBDR(IDLGER,IREJN)*VMTAGE
C
C  Convert the uncorrected (FTP conditions) exhaust from gpm to
C  gph and then adjust for no cold start component (idle occurs
C  under Hot FTP conditions).  CHGHRS=20.0
C
      BEFID5=BER*20.0
      IF(MY.LE.1985) BEFID5=BEFID5-22.88
      GOTO 20
C
C  In the pre-1984 model year case use the MOBILE4 revised idle EF
C  calculation algorithm.
C
   10 IGIDL=IDLPTR(MY,2,4)
      BEFID5=ZMLIDL(IGIDL,2,4,IREJN)+DRIDL(IGIDL,2,4,IREJN)*VMTAGE
C
C  Apply the I/M credit.
C
   20 BEFID5=BEFID5*CREDIM
C
C  Add in the Idle offset for Tampering
C
      BEFID5=BEFID5+OFFIDL(2,IDX,4)
C     IF(BEFID5.LT.0.0) BEFID5=0.0
C
C  Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID5=BEFID5*(OMTCF(IDX,2,4)/OMTCFF(IDX,2,4))
C
C  Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID5=BEFID5*SCIADJ(IDX,2,4)
C
C  Add RVP Correction Factor
C
      BEFID5=BEFID5*RVPCF(IDX,2,4)
C
C  Determine and apply the oxy fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
      BER=BER*CREDIM
      FUELCF=FUEL(ICY,2,4,MY,BER)
      BEFID5=BEFID5*FUELCF
      RETURN
      END
