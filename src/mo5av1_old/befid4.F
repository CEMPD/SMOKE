C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID4(MY,IDX,IP,IV,VMTAGE)
C
C  BEFID4 returns the basic idle emission factor adjusted for
C  I/M, methane and oxy fuel, depending on the case.
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
C    /MAXIMA/ MAXYRS
C    /OMTCOM/ OMTCF,OMTCFF
C    /REGION/ IREJN
C    /SPEED6/ SCIADJ
C
C  Output on return:
C
C    function: BEFID4
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
C                approx.
C
C
C  Notes:
C
C     BEFID4 is called by IDLCAL for IV = 5, 6, and 8.
C     BEFID4 calculates their idle emissions.
C     BEFID4 was created for MOBILE5
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
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
C
      REAL MTHIDL(2)
C
      DATA MTHIDL/0.9,1.5/
C
      BEFID4=0.0
      CREDIM=PCLEFT(MY,ICY,IP,IV)
C
C
C  In some MY cases use the MOBILE4 revised idle EF calculation
C  algorithm.
C
C  Compute base idle emission factor and then make corrections
C  for the methane offset.
C
      IGIDL=IDLPTR(MY,IP,IV)
      BEFID4=ZMLIDL(IGIDL,IP,IV,IREJN)+DRIDL(IGIDL,IP,IV,IREJN)*VMTAGE
C
C   Subtract off the Methane ( MC and HC only )
C   First check if not right IV or IP and skip out
C   Compute idle VOC, TOG, or NMOG, when necessary.
C
      IF(IV.EQ.8.AND.IP.EQ.1) THEN
          BEFID4 = BEFID4 - MTHIDL(IREJN)
C
          IF(NMHFLG.EQ.3) BEFID4 = BEFID4 * VOCFID(MY,IV)
          IF(NMHFLG.GE.4) BEFID4 = BEFID4 * TOGFID(MY,IV,OXYFLG)
C
C         Add back in the methane for THC and NMOG.
C
          IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.5)
     *       BEFID4 = BEFID4 + MTHIDL(IREJN)
          IF(BEFID4.LT.0.0) BEFID4=0.0
      ENDIF
C
C    Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID4=BEFID4*(OMTCF(IDX,IP,IV)/OMTCFF(IDX,IP,IV))
C
C    Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID4=BEFID4*SCIADJ(IDX,IP,IV)
C
C  No RVP adjustment is made to exhaust emissions for MC.
C  Determine and apply the oxy fuel correction factor,
C  as necessary (MC only).  Use BEFUN to get full ( non-idle )
C  emission factor based OXYCF.
C
      FUELCF=1.0
      IF(IV.EQ.8) THEN
         BER=BER*CREDIM
         FUELCF=FUEL(ICY,IP,IV,MY,BER)
      ENDIF
      BEFID4=BEFID4*FUELCF
      RETURN
      END
