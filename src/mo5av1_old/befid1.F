C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID1(MY,IDX,IP,IV)
C
C  BEFID1 returns the basic idle emission factor adjusted for I/M,
C  tampering, oxy fuel, and LEVs, depending on the case.
C
C  Called by IDLCAL.
C
C  Calls FUEL, PCLEFT and PCTLEV.
C
C  Input on call:
C
C    parameter list: MY,IDX,IP,IV,VMTAGE
C    common blocks:
C    /FLAGS2/ NEWFLG
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
C    function: BEFID1
C
C  Local array subscripts:
C
C  CHGHRS(4)  -  CHGHRS ( IVREV )
C  CRT(4,11)  -  CRT ( IV, LMY )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  CREDIM   R    I/M adjustment to idle emission factor
C  CRT      R    fraction of current standards vehicles
C  FAC      R    is the adjustment to the current standard idle
C                emission factor to reflect Tier 1 standards
C  FACL     R    is the adjustment to the current standard idle
C                emission factor to reflect LEV standards V
C  ILEV     I    Flag to indicate that LEVs exist for MY, IV combo
C  IGUTIO   I    model year index for untampered idle offsets (UTIDOF),
C                =1 for 1979-, =2 for 1980, ..., =14 for 1992+
C  IVVT     I    index for LDGV (=1) and LDGT1/2 (=2)
C  PTIER1   R    Percent tier 1 vehicles for that MY and IV
C
C
C  Notes:
C
C     BEFID1 is called by IDLCAL when IP is NOx and MY is greater
C     than or equal to 1977 for cars and greater than or equal to
C     1979 for trucks.  BEFID1 calculates idle emissions for IV = 1-3.
C     BEFID1 was created for MOBILE5.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG,ALHFLG
C
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IDLEQ3/ UTIDOF(14,3,2),LDVIDL(25)
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OPMOD1/ BFRCOF(8,8,3,6),BFSTEP(8,12,3,3),BFR50K(4,12,3,3)
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      DIMENSION FAC(12,3,3),FACL(12,3,3),CRT(12,3)
      DATA CRT/      0.6, 0.2,10*0.0,
     2               0.6, 0.2,10*0.0,
     3               1.0, 1.0, 0.50, 9*0.0/
      DATA  FAC/     0.34, 0.679, 0.775, 0.7, 8*0.663,
     2               0.4, 0.8, 10*1.0,
     3               0.16, 0.32, 10*0.4,
     1               0.381, 0.761, 0.868, 0.785, 8*0.743,
     2               0.688, 1.375, 1.482, 1.245, 8*1.127,
     3               0.212, 0.423, 10*0.529,
     1               0.458, 0.916, 1.052, 0.96, 8*0.914,
     2               0.676, 1.352, 1.556, 1.423, 8*1.356,
     3               0.422, 0.844, 0.97, 0.886, 8*0.844/
      DATA FACL/     0.0332, 0.0497, 0.0663, 0.0519, 0.0976,
     *          0.147, 0.193, 0.184, 0.18, 3*0.165,
     2               0.1, 0.15, 0.2, 0.26,
     *          0.49, 0.74, 0.97, 0.925, 0.9, 3*0.825,
     3               0.04, 0.06, 0.08, 0.054,
     *          0.1, 0.15, 0.196, 0.19, 0.19, 3*0.18,
     1               0.0371, 0.0557, 0.0743, 0.059,
     *          0.112, 0.168, 0.221, 0.211, 0.205, 3*0.188,
     2               0.113, 0.169, 0.225, 0.293,
     *          0.552, 0.834, 1.093, 1.042, 1.014, 3*0.929,
     3               0.0529, 0.0794, 0.106, 0.0772,
     *          0.143, 0.215, 0.28, 0.272, 0.272, 3*0.257,
     1          3*0.0, 0.12, 0.225, 0.339, 0.444, 0.425, 0.417, 3*0.385,
     2          3*0.0, 0.353, 0.664, 1.003, 1.315, 1.254, 1.22, 3*1.119,
     3          3*0.0, 0.219, 0.413, 0.623, 0.816, 0.78, 0.761, 3*0.701/
C
      BEFID1=0.0
      PTIER1=1.0
      ILEV=0
      LMY=MY-1993
      IF(LMY.LT.1) LMY=1
      IF(LMY.GT.12) LMY=12
C
      IF(LEVFLG.EQ.2.AND.IV.LE.2.AND.
     *    MY.GE.LEVYRS.AND.NEWFLG.LT.5) THEN
          ILEV=1
          PTIER1=PCTLEV(MY,IV,1)
      ENDIF
C
C  Assign ICY and determine the I/M credit.
C  Check the IMTIER flag for possible Special I/M credit assignment
C  PCLEFT looks at the IMFLAG and if IMFLAG=1 CREDIM is set to 1.0
C
      ICY=MY-IDX+MAXYRS
      CREDIM=PCLEFT(MY,ICY,IP,IV)
      IF(IMTIER.EQ.1.AND.MY.GE.1994) CREDIM=1.0
C
C  In some MY cases use the MOBILE4 revised idle EF calculation
C  algorithm.
C
      IGUTIO=MY-1978
      IF(IGUTIO.LT.1) IGUTIO=1
      IF(IGUTIO.GT.14) IGUTIO=14
      IVVT=IV
      IF(IV.EQ.3) IVVT=2
      BEFID1=UTIDOF(IGUTIO,IP,IVVT)
C
C  Estimate idle emissions of new standard vehicles based on
C  a ratio of new car standards.
C
C  BEFLEV is a function of BEFID1, LEVIMP and standards.
C  Save LDGV idle emissions to use for estimating new
C  standard truck emissions.
C
      IF(IV.EQ.1) LDVIDL(IDX) = BEFID1
      IF(MY.GE.1994.AND.NEWFLG.LT.5)
     *  BEFID1=BEFID1*CRT(LMY,IV)+LDVIDL(IDX)*FAC(LMY,IP,IV)
C
      IF(ILEV.EQ.1) THEN
C
          IF(IV.LT.3) BEFID1=BEFID1*PTIER1
     *                      +LDVIDL(IDX)*FACL(LMY,IP,IV)
C
C         for the levimp 2=TLEV now 3=TLEV
CCC       IF(IV.EQ.3)
CCC  *                BEFID1=BEFID1*PTIER1
CCC  *                      +BEFID1*(LEVIMP(3,LMY,4)*LDTAB(2,1)
CCC  *                              +LEVIMP(3,LMY,5)*LDTAB(2,2))
CCCC
CCCC *   +LDVIDL(IDX)*FACL(LMY,IP,IV)
CCCCC
      ENDIF
C
C  Apply the I/M Credit the same regardless of the value of ILEV
C
      BEFID1 = BEFID1*CREDIM
C
C  Add in the Idle offset for Tampering using PTIER1 if LEV.
C
      IF(ILEV.EQ.1) BEFID1=BEFID1+OFFIDL(IP,IDX,IV)*PTIER1
C
C  Add in the Idle offset for Tampering.
C
      IF(ILEV.EQ.0) BEFID1=BEFID1+OFFIDL(IP,IDX,IV)
C
C  Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID1=BEFID1*(OMTCF(IDX,IP,IV)/OMTCFF(IDX,IP,IV))
C
C  Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID1=BEFID1*SCIADJ(IDX,IP,IV)
C
C  Add RVP Correction Factor to gasoline vehicles and trucks
C
      BEFID1=BEFID1*RVPCF(IDX,IP,IV)
C
C  Determine and apply the oxy fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
C
      IF(ILEV.EQ.1) THEN
        CALL BEFUN(1992,IP,IV,VMTAGE,S1,S2,BEFADJ)
      ELSE
        CALL BEFUN(MY,IP,IV,VMTAGE,S1,S2,BEFADJ)
      ENDIF
C
      BEFADJ=BEFADJ*CREDIM
      FUELCF=FUEL(ICY,IP,IV,MY,BEFADJ)
      BEFID1=BEFID1*FUELCF
C
      RETURN
      END
