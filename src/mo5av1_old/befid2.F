C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID2(MY,IDX,IP,IV,VMTAGE)
C
C  BEFID2 returns the basic idle emission factor adjusted for I/M,
C  methane, tampering and oxy fuel, depending on the case.
C
C  Called by IDLCAL.
C
C  Calls IDLPTR, PCLEFT, TOGFID, VOCFID and FUEL.
C
C  Input on call:
C
C    parameter list: MY,IDX,IP,IV,VMTAGE
C    common blocks:
C    /FLAGS1/ OXYFLG
C    /FLAGS4/ NMHFLG
C    /IDLEQ1/ ZMLIDL,DRIDL
C    /IDLEQ3/ UTIDOF
C    /MAXIMA/ MAXYRS
C    /OMTCOM/ OMTCF,OMTCFF
C    /REGION/ IREJN
C    /RVPEX2/ RVPCF
C    /SPEED6/ SCIADJ
C    /TAMID3/ OFFIDL
C
C  Output on return:
C
C    function: BEFID2
C
C  Local array subscripts:
C
C  MTHIDL(2)  -  MTHIDL ( IREJN )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  CREDIM   R    I/M adjustment to idle emission factor
C  MTHIDL   R    methane offset to idle EF in gph (assumes methane is
C                approx. 6% of total HC over the affected MY range for
C                Bags 2 & 3)
C
C
C  Notes:
C
C     BEFID2 is called by IDLCAL for pre 1977 cars and for pre
C     1979 trucks.
C     BEFID2 calculates idle emissions for IV = 1-3
C     BEFID2 was created for MOBILE5.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      REAL MTHIDL(2)
C
      DATA MTHIDL/0.9,1.5/
C
      BEFID2=0.0
C
C  Assign ICY and determine the I/M credit.
C
      ICY=MY-IDX+MAXYRS
      CREDIM=PCLEFT(MY,ICY,IP,IV)
C
C  Use the MOBILE4 revised idle EF calculation algorithm and then
C  make corrections for, tampering offset and methane offset.
C
      IGIDL=IDLPTR(MY,IP,IV)
      BEFID2=ZMLIDL(IGIDL,IP,IV,IREJN)+DRIDL(IGIDL,IP,IV,IREJN)*VMTAGE
C
C    Apply the I/M Credit
C
       BEFID2 = BEFID2*CREDIM
C
C   Add in the Idle offset for Tampering
C
       BEFID2=BEFID2+OFFIDL(IP,IDX,IV)
C
C   HC Only Subtract off the Methane( LD gas and MC )
C   Compute idle VOC, TOG, or NMOG, when necessary.
C
      IF(IP.EQ.1) THEN
          BEFID2 = BEFID2 - MTHIDL(IREJN)
C
          IF(NMHFLG.EQ.3) BEFID2 = BEFID2 * VOCFID(MY,IV)
          IF(NMHFLG.GE.4) BEFID2 = BEFID2 * TOGFID(MY,IV,OXYFLG)
C
C         Add back in the methane for THC and NMOG.
C
          IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.5)
     *       BEFID2 = BEFID2 + MTHIDL(IREJN)
          IF(BEFID2.LT.0.0) BEFID2=0.0
      ENDIF
C
C    Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID2=BEFID2*(OMTCF(IDX,IP,IV)/OMTCFF(IDX,IP,IV))
C
C    Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID2=BEFID2*SCIADJ(IDX,IP,IV)
C
C    Add RVP Correction Factor to gasoline vehicles and trucks
C
      BEFID2 = BEFID2 * RVPCF(IDX,IP,IV)
C
C  Determine and apply the oxy fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
      BER=BER*CREDIM
      FUELCF=FUEL(ICY,IP,IV,MY,BER)
      BEFID2=BEFID2*FUELCF
      RETURN
      END
