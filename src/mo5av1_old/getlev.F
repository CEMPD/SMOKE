C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETLEV(ICY,MY,IP,IV,VMTAGE,AEF,IDX,BEEF)
C
C  Determine AEF, emissions using 'California' LEV standards
C
C  Called by BEF
C
C  Calls LEVEF, TOGFID, VOCFID and FUEL.
C
C  Input on call: BEEF
C
C    parameter list: ICY,MY,IP,IV,VMTAGE,AEF,IDX,BEEF
C    common blocks:
C
C    /FLAGS1/ OXYFLG
C    /FLAGS4/ NMHFLG
C    /OFFSET/ OFFCO,OFFMTH
C    /OMTCOM/ OMTCF,OMTTAM,OMTCFF
C    /MAXIMA/ MAXYRS
C    /LEVFLG/ LEVFLG,LEVYRS,IMFLG
C    /BASEQ9/ LEVIMP
C    /BASE12/ CAPIV1, CAPIV2
C    /BYMYC2/ BYBEF4
C
C  Output on return:
C
C    parameter list: AEF
C    common blocks: NONE
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
C  ICAP    I     Flag to indicate capped methane offset
C                1=no 2=yes
C  ICAR    I     index into LEV Classes
C  ILEV    I     index into LEVIMP process of LEVs into MOBILE
C  IX      I     index "A" or "B" for LDGT1
C  LEF     R     LEV emission rates
C  LEFADJ  R     THC EF for LEVs
C  LEVMTH  R     LEV Methane
C  LMY     I     index into LEVIMP model year process of LEVs
C  PCTNLV  R     I/M credit for EF that are not LEVs
C  PCTREM  R     LEV I/M credit
C  STRLEV  R     Stores (local) LEV EF's
C  TEMLEV  R     Tempory storage(local) for LEV EF's
C  WGT     R     Tempory storage(local) for Weighting of
C                "A" & "B"
C
C
C  Notes:
C
C  GETLEV added for MOBILE5, now it processes Low Emitting Vehicles.
C  GETLEV was modified in MOBILE5v13 to include a tier 1 with
C  special I/M cap on methane emissions.
C  GETLEV was modified in MOBILE5v15.1 PLUS to change the order
C  of the corrections factors to become more consistent with
C  the order in BEF.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      REAL LEVMTH, STRLEV, LEFADJ
      REAL LEVSTD, LEVIMP, LDTAB
C
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
      COMMON /BASE12/ CAPIV1,CAPIV2
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
C
      REAL LEF
      DIMENSION LEF(2),WGT(2),STRLEV(8,2)
C
C  ICAP is initialized to no cap then checked, when LEV with
C  special I/M for model year greater than LEVYRS.
C
C
      ICAP=1
      IF(MY.GE.LEVYRS.AND.IMFLG.EQ.2) THEN
         IF(IV.EQ.1) THEN
            IF(VMTAGE.GE.CAPIV1) ICAP=2
         ELSE
            IF(VMTAGE.GE.CAPIV2) ICAP=2
         ENDIF
      ENDIF
C
C  loop through the seven LEV types:
C  ( TIER 1, TLEV-Inter,TLEV, LEV-Inter,LEV, ULEV-Inter, and ULEV )
C  note: ILEV = 1 are "Current" vehicles and "8th" are ZEV's
C
      DO 30 ILEV = 2, 7
        JLEV = ILEV - 1
C
C  loop to fill compute the LDGT1 A's and LDGT1 B's
C
          DO 20 IX = 1, 2
C
          LEF(IX) = 0.0
          STRLEV(ILEV,IX) = 0.0
C
C  No "B" LDGV's
C
          IF(IV.EQ.1.AND.IX.EQ.2) GOTO 20
          ICAR=IV+(IX-1)
C
C  Get Exhaust LEV EF - Low Emission Vehicles Emission Factor
C
          CALL LEVEF(JLEV,IP,ICAR,BER,VMTAGE,MY)
C
          LEF(IX) = BER
C
C  Determine non LEV EF
C
          CALL BEFUN(MY,IP,IV,VMTAGE,S1,S2,BEFADJ)
C
C  Calculate LEV methane
C
          IF(IP.EQ.1) THEN
            LEVMTH=1./(1.-(OFFMTH(IDX,IV,ICAP)/BEFADJ))
            LEVMTH=LEVMTH-1.
            LEVMTH=LEVMTH*LEF(IX)
          ENDIF
C
C  Initialize I/M credit to 1 for special I/M Case.
C
          PCTREM = 1.0
          PCTNLV = 1.0
C
C  Calculate I/M Credit for not Special I/M Case
C
          IF(IMFLG.EQ.1) THEN
            PCTNLV = PCLEFT(MY,ICY,IP,IV)
            LEFADJ = BER
            IF(IP.EQ.1) LEFADJ = BER + LEVMTH
            PCTREM = 1.0 - ((1.0-PCTNLV)*(BEFADJ/LEFADJ))
          ENDIF
C
C  Calculate fuel adjustment for IP,IV,MY and not LEV EF.
C  Assume LDGV fuel adjustment (100% three way catalyst)
C
          BEFADJ = BEFADJ * PCTNLV
          FUELCF = FUEL(MY,IDX,IP,1,BEFADJ,ICAP)
C
C  Add methane to LEV emission factor and
C  Also adjust for operating mode, temperature and fuel.
C
          IF(IP.EQ.1) LEF(IX)=LEF(IX)+LEVMTH
          LEF(IX)=LEF(IX)*OMTCF(IDX,IP,IV)*FUELCF
C
          IF(IP.EQ.1) THEN
            LEF(IX)=LEF(IX)-LEVMTH
C
C  Adjust NMHC part first
C
            IF(NMHFLG.EQ.3) LEF(IX)=LEF(IX)*VOCFID(MY,IV)
            IF(NMHFLG.GE.4) LEF(IX)=LEF(IX)*TOGFID(MY,IV,OXYFLG)
C
C  Add back in the methane for THC and TOG
C
            IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.4)
     *        LEF(IX)=LEF(IX)+LEVMTH
          ENDIF
C
C  Apply I/M effect LEV emission factor
C
          LEF(IX) = LEF(IX) * PCTREM
C
C  Adjust for offset model (CO, LDGV & LDGT only).
C  Adjust offset for fuel cf as well.
C
          IF(IP.EQ.2)
     *    LEF(IX)=LEF(IX)+(OFFCO(IDX,IV)*FUELCF)
C
C  Store the LEV(IX) EF in STRLEV for calculating the AEF
C
          STRLEV(ILEV,IX) = LEF(IX)
C
   20   CONTINUE
   30 CONTINUE
C
C  LEV are merged into the MOBILE5 fleet at a rate shown in
C  array LEVIMP in BD05 starting in 1994 and ending in 2005.
C
      LMY=(MY-1994)+1
      IF (LMY.LE.1)  LMY=1
      IF (LMY.GT.12) LMY=12
C
C  Tier 1 uses the Orginal BERs ( basic emission rates ) from the
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
      IF(IV.EQ.2) THEN
         WGT(1) = LDTAB(1,1)
         WGT(2) = LDTAB(1,2)
      ENDIF
C
      DO 50 ILEV=1,7
        DO 40 IX=1,2
          ICAR=IV+(IX-1)
          AEF=AEF+STRLEV(ILEV,IX)*LEVIMP(ILEV,LMY,ICAR)*WGT(IX)
   40   CONTINUE
   50 CONTINUE
C
      JDX=MAXYRS-IDX+1
C
C     Save AEF as as by model year base weighted average
C     of BEF and LEV's.
C
      BYBEF4(IP,JDX,IV)=AEF
C
      RETURN
      END
