C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION BEF(MY,IDX,ICY,IP,IV,VMTAGE,IY)
C
C  BEF returns the basic emission factor adjustment for operating mode,
C  CO temperature offset, methane offset, temperature,
C  inspection / maintenance tampering offset and reformulated fuels.
C
C  Called by EFCALX and HCCALX.
C
C  Calls BEFUN, PCLEFT, FUEL, VOCFID, TOGFID, GETLEV and PCTLEV.
C
C  Input on call or from calls:
C
C    parameter list: MY,IDX,ICY,IP,IV,VMTAGE,IY
C    common blocks:
C    /BYMYC1/ IMPICK
C    /FLAGS1/ OXYFLG
C    /FLAGS2/ IMFLAG,IMTIER
C    /FLAGS4/ NMHFLG
C    /LEVBLK/ LEVFLG,IMFLG,LEVYRS,APPLEV
C    /OFFSET/ OFFCO,OFFMTH
C    /OMTCOM/ OMTCF,OMTTAM
C    /RFORM1/ RFGON
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
C  20-April-1994 @CSC ljn Request 2-400 added IMFLG test to correct
C  the LDGV's BYMY Credit/Reduction Print-out (OUTFMT=5)
C  22-April-1994 @CSC-trm Request #2-421 Added code for detergent gasoline
C                  adjustment. Added ICY to parameter list of call to FUEL
C  22-April-1994 @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEVs
C  24-August-1994 @CSC-pme request 446. Include file BASE12.I, 
C                 array OFFMTH's third dimension, and ICAP comment and
C                 definition were removed. 
C  October-3-1994 @CSC-pme request 2-472. A check for APPLEV variable to 
C  see whether or not to apply LEV stds to LDGT2 was included.
C  November-4-1994 @ CSC-bsg Request 0-472  Discovered that tampering should
C  be zero for guidance memo tier 1 emissions for LDGT2.
C  5-February-1996 @DynTel-yc request 2-613. Warning messages appeared during
C                  the compilation of BEF.FOR was fixed by setting flag IBEFSW
C                  on when BEFUN is called by BEF and off when BEFUN is  
C                  called by GETLEV.  SLOIM1 and SLOIM2 are set in BEFUN only
C                  when flag IBEFSW is on.
C  5-March-1996 @DynTel-yc request 2-616. Effect of RFG on total organic gases 
C                was fixed by calling TOGFID with OXYFLG=2 when RFGON is 
C                true and OXYFLG is 1.
C
      INCLUDE 'BYMYC1.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC5.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'OMTCOM.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'RFORM2.I'
      INCLUDE 'YEARS4.I'
C
      JDX=MAXYRS-IDX+1
C
C  Initialize BEF to the UN-Corrected value.
C  (5-February-1996) @ DynTel-yc Request 2-613 Set the flag IBEFSW=2(on) to
C                      set the SLOIM1 and SLOIM2 in subroutine BEFUN
C
      IBEFSW=2
      CALL BEFUN(MY,IP,IV,VMTAGE,BEF,IBEFSW)
      BEFADJ=BEF
C
C  Obtain % of emission factor LEFT after I/M program effect is
C  computed.
C
      PCTREM=1.0
C
C  New flag IMTIER = 1, if Tier 1 I/M Special Credit don't use PCLEFT
C  or       IMTIER = 0, if No Tier 1 I/M Special Credit use PCLEFT
C  if IMFLAG is 4 or 5 sets IMTIER to 1
C  if IMTIER is 1 and MY gt 1994 that's: Tier 1 Special I/M Credit
C  and no partial [I/M] credit remains
C
C  (20-April-1994) @CSC ljn Request 2-400 added IMFLG test to correct
C  the LDGV's BYMY Credit/Reduction Print-out (OUTFMT=5)
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEVs
C
C Non Max TIER1 I/M or Non Max LEV I/M or Before 1994
C
      IF((IV.LE.3).AND.
     *  (IMTIER.EQ.0.OR.IMFLG.EQ.1.OR.MY.LT.IY1994).OR.
C
C  also applies to the IV=4 HDGV's
C
     *  IV.EQ.4) THEN
        PCTREM=PCLEFT(MY,ICY,IP,IV)
      ENDIF
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
         FUELCF=FUEL(MY,IDX,IP,IV,ICY,BEFADJ,IY)
      ENDIF
C
C  Adjust for operating mode, temperature, oxy fuel, and, if HC,
C  methane offset and VOC, TOG, or NMOG FID, as necessary.
C
      BEF=BEF*OMTCF(IDX,IP,IV)*FUELCF
C
      IF(IP.EQ.1) THEN
         BEF = BEF - OFFMTH(IDX,IV)
C
C        Adjust NMHC part first
C  (5-March-1996) @DynTel-yc Request 2-616  Effect of RFG on TOG is controlled
C                   before calling TOGFID when NMHFLG=4 or 5.
C
         IF(NMHFLG.EQ.3) BEF=BEF*VOCFID(MY,IV)
         IF(NMHFLG.EQ.4.OR.NMHFLG.EQ.5) THEN
           IF(RFGON .AND. OXYFLG.EQ.1) THEN
             BEF=BEF*TOGFID(MY,IV,2)
           ELSE
             BEF=BEF*TOGFID(MY,IV,OXYFLG)
           ENDIF
         ENDIF
C
C        Add back in the methane for THC and TOG.
C
         IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.4) BEF=BEF+OFFMTH(IDX,IV)
      ENDIF
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
C  No tampering effect adjustment needed for non-Tier 1 (Pre-1994) vehicles
C  and HDGV's (IV=4), ( they receive full tampering effect )
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEVs
C
      IF(MY.LT.IY1994.OR.IV.EQ.4) THEN
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
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEVs
C
      IF(LEVFLG.EQ.2.AND.MY.GE.LEVYRS) THEN
         IF(IV.LE.2.OR.(IV.EQ.3.AND.APPLEV.EQ.2))THEN
           CALL GETLEV(ICY,MY,IP,IV,VMTAGE,AEF,IDX,BEEF,IY)
           BEF=AEF
C
           ATIER=PCTLEV(MY,IV,1)
           AZEV =PCTLEV(MY,IV,8)
           ALEV = 1.0 - (ATIER + AZEV)
C
         ENDIF
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
