C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEF(MY,IDX,ICY,IP,IV,VMTAGE)
C
C  BEF returns the basic emission factor adjustment for operating mode,
C  CO temperature offset, methane offset, temperature,
C  inspection / maintenance tampering offset and reformulated fuels.
C
C  Called by BEFIDL, EFCALX and HCCALX.
C
C  Calls PCLEFT, BEFUN and FUEL.
C
C  Input on call or from calls:
C
C    parameter list: MY,IDX,ICY,IP,IV,VMTAGE
C    common blocks:
C    /BASE12/ CAPIV1, CAPIV2
C    /BYMYC1/ IMPICK
C    /FLAGS1/ OXYFLG
C    /FLAGS2/ IMFLAG,IMTIER
C    /FLAGS4/ NMHFLG
C    /OFFSET/ OFFCO,OFFMTH
C    /OMTCOM/ OMTCF,OMTTAM,OMTCFF
C    /LEVBLK/ LEVFLG,LEVYRS
C    /RFORM1/ RFGFLG,IASTM,IPHASE
C    /RFORM2/ RFGSYR
C    /YEARS4/ IY1994
C
C  Output on return:
C
C    function: BEF
C    common blocks:
C    /BYMYC2/ BYBEF4,BYTAM
C    /BYMYC5/ BYIMCR,BYIMRE
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  AEF      R    Comes from GETLEV, is like BEF.
C  BEEF     R    Value holder for BEF for passing it in subroutines
C  ICAP     I    Flag to indicate capped methane offset; 1=no 2=yes
C  KINK50   I    key to whether or not to use the 50K+ deterioration rates:
C                1 = no:  no 50K kink in the curve or do not reach it
C                2 = yes: have 50K kink & do reach it (mileage exceeds 50K)
C  PCTREM   R    Temporary variable used to create complement.
C
C
C  Notes:
C
C  BEF was modified for MOBILE4 In-House Version 02 to save by model year
C  computations for the BYMY Tables.
C  BEF was changed in MOBILE4.1 to calculate VOC, NMOG and TOG
C  and adjust for oxy fuels.
C  BEF was updated to include LDGT ef kinks in MOBILE5.
C  BEF was split into BEF, BEFUN, BEFHDV for MOBILE5 v2.
C  BEF was modified in MOBILE5v13 to include a capped OFFMTH.
C
C
      CHARACTER*1 ASTMCL
      INTEGER PRTFLG,HCFLAG,COLDFG
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER RFGSYR,ALHFLG
      INTEGER RFGFLG,WINFLG
      LOGICAL RFGON
C
      COMMON /BASE11/ SLOIM1, SLOIM2
      COMMON /BASE12/ CAPIV1, CAPIV2
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C  ICAP is initialized to no cap then checked, when Tier 1 with
C  special I/M for IV=1,2 and model is greater than IY1994
C
      ICAP=1
      IF(MY.GE.IY1994.AND.IMTIER.EQ.1.AND.IV.LE.2) THEN
         IF(IV.EQ.1) THEN
            IF(VMTAGE.GE.CAPIV1) ICAP=2
         ELSE
            IF(VMTAGE.GE.CAPIV2) ICAP=2
         ENDIF
      ENDIF
C
C  Initialize BEF to the UN-Corrected value.
C  and the SLOIM1 and SLOIM2 from the S1 and S2
C
      SLOIM1 = 0.0
      SLOIM2 = 0.0
C
      CALL BEFUN(MY,IP,IV,VMTAGE,S1,S2,BEF)
      BEFADJ=BEF
C
      SLOIM1 = S1
      SLOIM2 = S2
C
C  Obtain % of emission factor LEFT after I/M program effect is
C  computed.
C
      PCTREM=1.0
C
C  Applys to IV=1,2 LDGV and LDGT1
C
C  New flag IMTIER = 1, if Tier 1 I/M Special Credit don't use PCLEFT
C  or       IMTIER = 0, if No Tier 1 I/M Special Credit use PCLEFT
C  if IMFLAG is 4 or 5 sets IMTIER to 1
C  if IMTIER is 1 and MY gt 1994 that's: Tier 1 Special I/M Credit
C  and no partical [I/M] credit remains
C
      IF(IMTIER.EQ.0.OR.MY.LT.IY1994) THEN
         IF(IV.LE.2)PCTREM=PCLEFT(MY,ICY,IP,IV)
      ENDIF
C
C  apply to the IV=3 LDGT2's or IV=4 HDGV's
C
      IF(IV.EQ.3.OR.IV.EQ.4) PCTREM=PCLEFT(MY,ICY,IP,IV)
C
C  Determine the fuel correction factor, as necessary.
C  Use PCLEFT adjusted BEFUN value as input to FUEL.
C
      FUELCF=1.0
C
C  Tier 1 vehicles use current standard vehicle ef's to
C  calculate fuel adjustments.
C
      IF(IV.LE.4.OR.IV.EQ.8) THEN
         BEFADJ=BEFADJ*PCTREM
         FUELCF=FUEL(MY,IDX,IP,IV,BEFADJ,ICAP)
      ENDIF
C
C  Adjust for operating mode, temperature, oxy fuel, and, if HC,
C  methane offset and VOC, TOG, or NMOG FID, as necessary.
C
      BEF=BEF*OMTCF(IDX,IP,IV)*FUELCF
C
      IF(IP.EQ.1) THEN
         BEF = BEF - OFFMTH(IDX,IV,ICAP)
C
C        Adjust NMHC part first
C
         IF(NMHFLG.EQ.3) BEF=BEF*VOCFID(MY,IV)
         IF(NMHFLG.EQ.4.OR.NMHFLG.EQ.5) BEF=BEF*TOGFID(MY,IV,OXYFLG)
C
C        Add back in the methane for THC and TOG.
C
         IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.4) BEF=BEF+OFFMTH(IDX,IV,ICAP)
      ENDIF
      JDX=MAXYRS-IDX+1
C
C  Use the PerCentage of emission factor LEFT after I/M program effect is
C  computed.  Save the complement as the I/M Credit and apply the Credit to
C  BEF to get the I/M Reduction, if the I/M By Model Year Table option has
C  been selected.
C
      IF(IMPICK.EQ.2.AND.IV.LE.4) THEN
         BYIMCR(IP,JDX,IV)=1.0-PCTREM
         BYIMRE(IP,JDX,IV)=BEF*BYIMCR(IP,JDX,IV)
      ENDIF
C
      BEF=BEF*PCTREM
C
      IF(IV.GT.4) THEN
         BYBEF4(IP,JDX,IV)=BEF
         RETURN
      ENDIF
C
C  Adjust for offset model (CO, LDGV & LDGT only) and save as untampered
C  by model year base emission factor. Adjust offset for oxy fuel, too.
C
      IF(IV.LE.3.AND.IP.EQ.2) BEF=BEF+OFFCO(IDX,IV)*FUELCF
      BYBEF4(IP,JDX,IV)=BEF
C
C  No tampering effect adjustment needed for non-Tier 1
C  (Pre-1994) vehicles LDGT2'S (IV=3) and HDGV's (IV=4),
C  ( they recieve Full tampering effect )
C
      IF(MY.LT.IY1994.OR.IV.EQ.3.OR.IV.EQ.4) THEN
         BYTAM(IP,JDX,IV)=OMTTAM(IDX,IP,IV)*FUELCF
         BEF = BEF+BYTAM(IP,JDX,IV)
         RETURN
      ENDIF
C
C  load BEF into BEEF for LEV Vehicle Type processing
C
      BEEF=BEF
C
C  Compute the LEV Low Emission Vehicles
C  Corrected basic emission factor AEF
C
C  Initialize sizes for Tier 1 (non-LEV) case
C
      ATIER = 1.0
      ALEV  = 0.0
C
C  If LEV program, determine new BEF and Tier 1/LEV sizes.
C  LEV does not apply to LDGT's (IV=3)
C
      IF(LEVFLG.EQ.2.AND.IV.LE.2.AND.MY.GE.LEVYRS) THEN
         CALL GETLEV(ICY,MY,IP,IV,VMTAGE,AEF,IDX,BEEF)
         BEF=AEF
C
         ATIER=PCTLEV(MY,IV,1)
         AZEV =PCTLEV(MY,IV,8)
         ALEV = 1.0 - (ATIER + AZEV)
C
      ENDIF
C
C  When special I/M is selected, zero out the
C  corresponding, Tier 1/LEV size.
C
      IF(IMTIER.EQ.1) ATIER = 0.0
      IF(IMFLG.EQ.2)  ALEV  = 0.0
C
C  Determine adjusted tampering offset
C
      BYTAM(IP,JDX,IV)=
     * OMTTAM(IDX,IP,IV) * (ATIER + ALEV) * FUELCF
C
      BEF=BEF+BYTAM(IP,JDX,IV)
C
      RETURN
      END
