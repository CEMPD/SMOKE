C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETLEV(ICY,MY,IP,IV,VMTAGE,AEF,IDX,BEEF,IY)
C
C  Determine AEF, emissions using 'California' LEV standards
C
C  Called by BEF
C
C  Calls LEVEF, BEFUN, PCLEFT, FUEL, VOCFID and TOGFID.
C
C  Input on call:
C
C    parameter list: ICY, MY, IP, IV, VMTAGE, IDX, BEEF, IY
C    common blocks:
C
C    /FLAGS1/ OXYFLG
C    /FLAGS2/ IMTIER
C    /FLAGS4/ NMHFLG
C    /OFFSET/ OFFCO,OFFMTH
C    /OMTCOM/ OMTCF,OMTTAM
C    /MAXIMA/ MAXYRS
C    /LEVFLG/ LEVFLG,LEVYRS,IMFLG
C    /BASEQ9/ LEVIMP
C
C  Output on return:
C
C    parameter list: AEF
C    common blocks:
C
C    /BYMYC2/ BYBEF4
C
C  Local array subscripts:
C
C  LEF(2)       -  LEF  ( IX )
C  STRLEV(8,2)  -  STRLEV  ( ILEV, IX )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  BEFFCF  R     Tier1 Basic emission factor adjusted for I/M effects.
C                Used for fuel correction factor calculation.
C  BER     R     Basic Emission Rate BER, Variable passing into LEVEF
C  ICAR    I     Index into LEV Classes
C  ILEV    I     Index into LEVIMP process of LEVs into MOBILE
C  IMTEMP  R     Temporary hold for IMTIER flag for call to BEFUN, this
C                allows the BEF returned to be Non MAX I/M for the LEV's
C  IX      I     Index "A" or "B" for LDGT1
C  LEF     R     LEV emission rates
C  LEVMTH  R     LEV Methane
C  LMY     I     Index into LEVIMP model year process of LEVs
C  PCTNLV  R     I/M credit for EF that are not LEVs
C  PCTREM  R     LEV I/M credit
C  STRLEV  R     Stores (local) LEV EF's
C  TEMLEV  R     Temporary storage(local) for LEV EF's
C  WGT     R     Temporary storage(local) for Weighting of
C                "A" & "B"
C  Notes:
C
C  GETLEV added for MOBILE5, now it processes Low Emitting Vehicles (LEV's).
C  GETLEV was modified in MOBILE5v13 to include a tier 1 with
C  special I/M cap on methane emissions.
C  GETLEV was modified in MOBILE5v15.1 PLUS to change the order
C  of the corrections factors to become more consistent with
C  the order in BEF.
C  20-April-94 @ CSC ljn Request 2-400 TIER1 I/M. save and set flag IMTIER to
C              FALSE before calling BEFUN. This forces BEFUN to return the 
C              Non Max TIER1 I/M Basic Emission Rate to variable BEFADJ.
C
C  22-April-1994 @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C  6-July-1994 @CSC-rwk Request #0-448 Added code for detergent gasoline
C              adjustment. Added ICY to parameter list of call to FUEL.  This was
C              missed on the original task, Request #2-421.
C
C  8-July-1994 @CSC pme Request 2-446
C              LEV's mileage cap removed. 
C              Include file BASE12.I and comment were deleted.
C
C  12-August-1994 @CSC-tmm Request 2-460
C                 Added operating mode, temperature and fuel correction factors
C                 that were inadvertently left out in the formula used to 
C                 determine the methane fraction for LEVs. This code was 
C                 supplied directly by the EPA.
C
C  24-August-1994 @CSC pme Request 2-446. The third dimension of OFFMTH array
C                 and ICAP variable were removed. 
C  5-February-1996 @DynTel request 2-613. Warning messages appeared during the
C                  compilation of BEF.FOR was fixed by setting flag IBEFSW
C                  on when BEFUN is called by BEF and off when BEFUN is  
C                  called by GETLEV.  SLOIM1 and SLOIM2 are set in BEFUN only
C                  when flag IBEFSW is on.
C
      INCLUDE 'BASEQ9.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'OMTCOM.I'
C
      REAL BER, LEF, LEVMTH, STRLEV
      DIMENSION LEF(2), WGT(2), STRLEV(8,2)
C
C  loop through the seven LEV types:
C  ( TIER 1, TLEV-Inter,TLEV, LEV-Inter,LEV, ULEV-Inter, and ULEV )
C  note: ILEV = 1 are "Current" vehicles and "8th" are ZEV's
C
      DO 30 ILEV = 2, 7
        JLEV = ILEV - 1
C
C  loop to fill compute the LDGT1 A's & B's Weighting factors
C
        DO 20 IX = 1, 2
C
          LEF(IX) = 0.0
          STRLEV(ILEV,IX) = 0.0
C
C  There are No LDGV's "B" Weighting factors
C
          IF (IV.EQ.1.AND.IX.EQ.2) GOTO 20
C
          ICAR=IV+(IX-1)
          IF (IV.EQ.3) ICAR=IV+IX
C
C  Get Exhaust LEV EF - Low Emission Vehicles Emission Factor
C
          CALL LEVEF(JLEV,IP,ICAR,BER,VMTAGE)
C
C  Determine TIER1 vehicles (non LEV) EF
C
          IMTEMP=IMTIER
C
C  (20-April-94) @ CSC ljn Request 2-400 TIER1 I/M ......
C  Set IMTIER flag to a Non TIER1 I/M, so the call to BEFUN returns
C  the Non MAX TIER1 I/M value in BEFADJ
C  (5-February-1996) @ DynTel Request 2-613 The flag IBEFSW is set as 1 and 
C                      SLOIM1 and SLOIM2 are not set in BEFUN.
C
          IMTIER=0
          IBEFSW=1
          CALL BEFUN(MY,IP,IV,VMTAGE,BEFADJ,IBEFSW)
C
C  reset IMTIER
C
          IMTIER=IMTEMP
C
C  Calculate LEV methane
C
C  Calculate fuel adjustment for IP,IV,MY and non-LEV EF.
C  Assume LDGV fuel adjustment (100% three way catalyst)
C
          BEFFCF = BEFADJ * PCLEFT(MY,ICY,IP,IV)
          FUELCF = FUEL(MY,IDX,IP,1,ICY,BEFFCF,IY)
C
C  Adjust for operating mode and fuel
C
          BEFADJ = BEFADJ * OMTCF(IDX,IP,IV) * FUELCF
          LEF(IX) = BER *  OMTCF(IDX,IP,IV) * FUELCF
C 
C  Calculate and add the LEV methane adjustment
C
          IF (IP.EQ.1) THEN
            LEVMTH = LEF(IX) *
     *               ((1./(1.-(OFFMTH(IDX,IV)/BEFADJ)))-1.)
            LEF(IX) = LEF(IX) + LEVMTH
          END IF
C
C  Initialize I/M credit to 1 for special I/M Case.
          PCTREM = 1.0
          PCTNLV = 1.0
C  Calculate I/M Credit for not Special I/M Case
C
          IF (IMFLG.EQ.1) THEN
            PCTNLV = PCLEFT(MY,ICY,IP,IV)
            PCTREM = 1.0 - ((1.0-PCTNLV)*(BEFADJ/LEF(IX)))
          END IF
C
C  Determine appropriate HC species
C
          IF (IP.EQ.1) THEN
            LEF(IX)=LEF(IX)-LEVMTH
C
C  Adjust NMHC part first
C
            IF (NMHFLG.EQ.3) LEF(IX)=LEF(IX)*VOCFID(MY,IV)
            IF (NMHFLG.GE.4) LEF(IX)=LEF(IX)*TOGFID(MY,IV,OXYFLG)
C
C  Add back in the methane for THC and TOG
C
            IF (NMHFLG.EQ.1.OR.NMHFLG.EQ.4) LEF(IX)=LEF(IX)+LEVMTH
          ENDIF
C
C  Apply I/M effect LEV emission factor
C
          LEF(IX) = LEF(IX) * PCTREM
C
C  Adjust for offset model (CO, LDGV & LDGT only).
C  Adjust offset for fuel cf as well.
C
          IF (IP.EQ.2) LEF(IX)=LEF(IX)+(OFFCO(IDX,IV)*FUELCF)
C
C  Store the LEV(IX) EF in STRLEV for calculating the AEF
C
          STRLEV(ILEV,IX) = LEF(IX)
C
   20   CONTINUE
C
   30 CONTINUE
C
C  LEV are merged into the MOBILE5 fleet at a rate shown in
C  array LEVIMP in BD05 starting in 1994 and ending in 2005.
C
      LMY=(MY-1994)+1
      IF (LMY.LE.1)  LMY=1
      IF (LMY.GT.12) LMY=12
C
C  Tier 1 uses the Original BERs ( basic emission rates ) from the
C  Sub. BEF, but uses the Model Year 1996 BERs stored in "BEEF"
C
      STRLEV(1,1) = BEEF
      STRLEV(1,2) = STRLEV(1,1)
C
C  Initialize and add up to return final AEF value
C
      AEF    = 0.0
      WGT(1) = 1.0
      WGT(2) = 0.0
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C
      IF (IV.EQ.2.OR.IV.EQ.3) THEN
        WGT(1) = LDTAB(IV-1,1)
        WGT(2) = LDTAB(IV-1,2)
      END IF
C
      DO 50 ILEV=1,7
C
        DO 40 IX=1,2
          ICAR=IV+(IX-1)
          IF (IV.EQ.3) ICAR=IV+IX
          AEF=AEF+STRLEV(ILEV,IX)*LEVIMP(ILEV,LMY,ICAR)*WGT(IX)
   40   CONTINUE
C
   50 CONTINUE
C
      JDX=MAXYRS-IDX+1
C
C     Save AEF as by model year base weighted average
C     of BEF and LEV's.
C
      BYBEF4(IP,JDX,IV)=AEF
      RETURN
      END
