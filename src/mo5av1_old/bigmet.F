C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGMET(IV,MY,IDX,JDX)
C
C  BIGMET computes the additive methane offset to the composite
C  HC emission factor with and without a cap for tier 1 and lev.
C
C  Called by BIGCFX.
C
C  Calls EXMPTR, IVPTRT, IEVPTR.
C
C  Input on call:
C
C    parameter list: IV,MY,IDX
C    common blocks:
C    /BASE12/ CAPIV1,CAPIV2
C    /CUMCOM/ CUMMIL
C    /NMETH1/ ZDMTH1,ZDMTH2,ZDMTH3
C    /NMETH2/ MTHMYG
C    /REGION/ IREJN
C    /TEMPC3/ CHFRAC,TCF
C
C  Output on return:
C
C    common blocks:
C    /OFFSET/ OFFMTH
C
C  Local array subscripts
C
C  CRBFNJ(2) - CRBFNJ(IFDS)
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  CRBFNJ   R    fraction of carb and fuel injected vehicles
C  IFDS     I    fuel deleivery system index,
C                1=carb, 2=injection
C  IVTEMP   I    temporary IV index into ZDMTHi arrays
C  OFF      R    by-bag methane offset
C  OFFMAX   R    maximum methane offset calculated in stream
C  OMT      R    unweighted product of operating mode vmt
C                fraction times temperature correction factor
C  VMT      R    local mileage version of VMTAGE
C  VMTTMP   R    temporary storage of capped VMTAGE
C
C  Notes:
C
C  BIGMET was modified for MOBILE4.1 to include methane offset
C  ZML and DR by fuel delivery system and region.
C  BIGMET was modified for MOBILE5v13 to include a cap for
C  methane offset for lev and tier 1 with special I/M.
C
C
      COMMON /BASE12/ CAPIV1,CAPIV2
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /INJECT/ TBI(13,4),PFI(13,4),MAXFIY,MINFIY
      COMMON /NMETH1/ ZDMTH1(2,3,2,7,3,2),ZDMTH2(2,7,2,2),ZDMTH3(7,3,2)
      COMMON /NMETH2/ MYGMTH(7,8,2),MTHMYG
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
C
      DIMENSION CRBFNJ(2)
C
C  Initialize.
C
      OFFMTH(IDX,IV,1)=0.0
      OFFMTH(IDX,IV,2)=0.0
      VMT=CUMMIL(JDX,IV)/10000.
C
C  Set myg pointer.
C
      CALL EXMPTR(MY,IV)
C
C  Determine carb/fuel injections fraction.
C  (LDGV/T only)
C
      IF(IV.LE.3) THEN
      IG5=IEVPTR(MY,IV,2)
      CRBFNJ(2)=TBI(IG5,IV)+PFI(IG5,IV)
      CRBFNJ(1)=1.0-CRBFNJ(2)
      ENDIF
C
C  Set and loop for each bag.
C
      NBAG=3
      IF(IV.EQ.4.OR.IV.EQ.7) NBAG=1
C
      DO  20 IB=1,NBAG
C
C  Set the mode/temp factor.
C  (HDGV and HDDV are assummed 100% hot stabilized.)
C
      OMT=1.0
      IF(IV.NE.4.AND.IV.NE.7) OMT=CHFRAC(IB)*TCF(IB)
C
C  Determine LDGV, LDGT1, and LDGT2 offsets.
C
      IF(IV.LE.3) THEN
      OFF=0.
      OFFMAX=0.
      DO 10 IFDS=1,2
      OFF=OFF+(ZDMTH1(1,IB,IFDS,MTHMYG,IV,IREJN)
     *       + ZDMTH1(2,IB,IFDS,MTHMYG,IV,IREJN)*VMT)
     *       * CRBFNJ(IFDS)
      VMTTMP=CAPIV2
      IF(IV.EQ.1) VMTTMP=CAPIV1
      OFFMAX=OFFMAX+(ZDMTH1(1,IB,IFDS,MTHMYG,IV,IREJN)
     *       + ZDMTH1(2,IB,IFDS,MTHMYG,IV,IREJN)*VMTTMP)
     *       * CRBFNJ(IFDS)
  10  CONTINUE
      ENDIF
C
C  Determine HDGV and MC offsets.
C
      IF(IV.EQ.4.OR.IV.EQ.8) THEN
      IVTEMP=1
      IF(IV.EQ.8) IVTEMP=2
      OFF=ZDMTH2(1,MTHMYG,IVTEMP,IREJN)
     *   +ZDMTH2(2,MTHMYG,IVTEMP,IREJN)*VMT
      OFFMAX=OFF
      ENDIF
C
C  Determine LDDV, LDDT and HDDV offsets.
C
      IF(IV.GE.5.AND.IV.LE.7) THEN
      IVTEMP=IV-4
      OFF=ZDMTH3(MTHMYG,IVTEMP,IREJN)
      OFFMAX=OFF
      ENDIF
C
C  Summate the by-bag methane offsets.
C
      OFFMTH(IDX,IV,1)=OFFMTH(IDX,IV,1)+OFF*OMT
      OFFMTH(IDX,IV,2)=OFFMTH(IDX,IV,2)+OFFMAX*OMT
C
   20 CONTINUE
C
      RETURN
      END
