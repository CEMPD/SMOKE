C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFIDL(MY,IDX,IP,IV,VMTAGE)
C
C  BEFIDL returns the basic idle emission factor adjusted for I/M,
C  methane, tampering and oxy fuel, depending on the case.
C
C  Called by IDLCAL.
C
C  Calls BEFUN, FUEL, IBFPTR, IDLPTR, PCLEFT, PCTLEV, TOGFID,
C  and VOCFID.
C
C  Input on call:
C
C    parameter list: MY,IDX,IP,IV,VMTAGE
C    common blocks:
C    /FLAGS1/ OXYFLG
C    /FLAGS2/ NEWFLG,IMTIER
C    /FLAGS4/ NMHFLG
C    /IDLEQ3/ UTIDOF
C    /LEVBLK/ LEVYRS
C    /MAXIMA/ MAXYRS
C    /OFFSET/ OFFCO
C    /OMTCOM/ OMTCF,OMTCFF
C    /OPMOD1/ BFRCOF,BFSTEP,BFR50K
C    /OPMOD2/ MYSTEP
C    /REGION/ IREJN
C    /RVPEX2/ RVPCF
C    /SPEED6/ SCFIDL,SCIADJ
C    /TAMID3/ OFFIDL
C
C  Output on return:
C
C    function: BEFIDL
C
C  Local array subscripts:
C
C  CHGHRS(4)  -  CHGHRS ( IVREV )
C  COLIM(14)  -  COLIM  ( IGUTIO )
C  D50   (3)  -  D50    ( I23F )
C  DBF   (3)  -  DBF    ( I23F )
C  MTHIDL(2)  -  MTHIDL ( IREJN )
C  ZBF   (3)  -  ZBF    ( I23F )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  BER      R    untampered basic emission rate
C  BFGT50   R    mileage greater than 50,000 ( BFMILE = 6.0 ->
C                BFGT50 = 1.0 )
C  BFMILE   R    cumulative mileage in ICY for IV from MY
C                scaled by 10**-4
C  CHGHRS   R    factor for converting BER units from gpm to gph = FTP
C                standard speed for IVREV ( BER is returned from BEF
C                uncorrected = assumes FTP conditions)
C  COLIM    R    CO Idle EF LBs for each IGUTIO MYG; keeps low vmt
C                cases > 0.0
C  CREDIM   R    I/M adjustment to idle emission factor
C  CRT      R    fraction of current standards vehicles
C  D50      R    BF deterioration rate for beyond 50K miles
C  DBF      R    BF deterioration rate for the first 50K miles
C  FAC      R    is the adjustment to the current standard idle
C                emission factor to reflect Tier 1 standards
C  FACL     R    is the adjustment to the current standard idle
C                emission factor to reflect LEV standards V
C  I23F     I    index for Bag 2 ( = 1 ), Bag 3 ( = 2 ), and FTP
C                ( = 3 )
C  IB       I    index for Bag 1 ( = 1 ), Bag 2 ( = 2 ), Bag 3
C                ( = 3 ), and FTP ( = 4 )
C  ILEV     I    Flag to indicate that LEVs exist for MY.
C  IGUTIO   I    model year index for untampered idle offsets
C                ( UTIDOF ), = 1 for 1979-, = 2 for 1980, ..., = 14
C                for 1992+
C  IVVT     I    index for LDGV ( = 1 ) and LDGT1/2 ( = 2 )
C  MTHIDL   R    methane offset to idle EF in gph (assumes methane is
C                approx. 6% of total HC over the affected MY range for
C                Bags 2 & 3)
C  OMCFF    R    operating mode correction factor at FTP
C                conditions
C  OMCFH    R    operating mode correction factor at Hot FTP
C                conditions
C  PTIER1   R    Percent tier 1 vehicles for that MY.
C  RATIO    R    ratio of OMCF at Hot FTP to OMCF at FTP
C  ZBF      R    BF zero mile level
C
C
C  Notes:
C
C   OMTCFF was added for MOBILE4.1, version D02  [ OMTCFF is
C   OMTCF @ FTP 100% stable Mode /0/0/0/, 75 deg F and 9.0 RVP ]
C   BEFIDL was modified for MOBILE4.1.
C   BEFIDL was modified in MOBILE5 to include kinks for LDGT1/2 and
C   update mapping assumptions.
C   modified to include a methane offset for LEVs
C   BEFIDL was modified for MOBILE5v13 to include calls to PCTLEV
C
C
      DIMENSION FAC(12,3,3),FACL(12,3,3),CRT(12,3)
      REAL LEVMTH
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER PRTFLG,HCFLAG,COLDFG,ALHFLG
C
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /IDLEQ3/ UTIDOF(14,3,2),LDVIDL(25)
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OPMOD1/ BFRCOF(8,8,3,6),BFSTEP(8,12,3,3),BFR50K(4,12,3,3)
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      DIMENSION DBF(3),ZBF(3),D50(3),CHGHRS(4),COLIM(14)
      REAL MTHIDL(2)
C
      DATA MTHIDL/0.9,1.5/,CHGHRS/19.6,19.6,19.6,20.0/
      DATA COLIM/303.37,56.71,22.56,25.63,35.17,30.36,31.77,
     *            26.23,25.08,24.07,23.71,23.24,22.38,21.78/
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
      BER=0.0
      BEFIDL=0.0
      BEFLEV=0.0
      LMY=MY-1993
      IF(LMY.LT.1) LMY=1
      IF(LMY.GT.12) LMY=12
C
C  Obtain an untampered BER by using the new subroutine BEFUN.
C
      CALL BEFUN(MY,IP,IV,VMTAGE,S1,S2,BER)
C
C  Assign the ILEV flag if the MY is greater than or equal 1994
C  and the LEVFLG is set.
C
      ILEV=0
      IF(LEVFLG.EQ.2.AND.MY.GE.LEVYRS) ILEV=1
      IF(IV.GE.3) ILEV=0
C
C  Assign ICY and determine the I/M credit.
C
      ICY=MY-IDX+MAXYRS
C
C  Check the IMTIER flag for possible Special I/M credit assignment
C  PCLEFT looks at the IMFLAG and if IMFLAG=1 CREDIM is set to 1.0
C
      CREDIM=PCLEFT(MY,ICY,IP,IV)
      IF(IMTIER.EQ.1.AND.MY.GE.1994) CREDIM=1.0
C
C  Revised HC & CO calculation algorithm case:
C
C  LEV emissions starting 1994
C
      PTIER1=1.0
C
      IF(ILEV.EQ.1) PTIER1=PCTLEV(MY,IV,1)
C
C  Set up various bag fraction pointers and mileage values.
C
      MYSTEP=2
      IF(MY.LT.1981.OR.(MY.LT.1984.AND.IV.GT.1)) MYSTEP=1
C
      IGBF=IBFPTR(MY,IP,IV)
      JDX=MAXYRS-IDX+1
C
      DO 25 I23F=1,3
      D50(I23F)=0.0
   25 CONTINUE
C
      BFMILE=CUMMIL(JDX,IV)*.0001
      BFGT50=0.0
C
C  Determine Bag Fraction coefficients for Bag 2, Bag 3, and FTP.
C  Use appropriate coefficients for 1981+ LDGV and 1984+ LDGT1/2.
C
      IF((IV.EQ.1.AND.MY.GE.1981).OR.(IV.GE.2.AND.MY.GE.1984))
     *   GOTO 35
C
      DO 30 I23F=1,3
      IB=I23F*2+2
      ZBF(I23F)=BFRCOF(IB-1,IGBF,IP,IV)
      DBF(I23F)=BFRCOF(IB,IGBF,IP,IV)
   30 CONTINUE
      GOTO 45
C
   35 DO 40 I23F=1,3
      IB=I23F*2+2
      ZBF(I23F)=BFSTEP(IB-1,IGBF,IP,IV)
      DBF(I23F)=BFSTEP(IB,IGBF,IP,IV)
   40 CONTINUE
C
C  Use kinked Bag Fractions for 1981+ LDGV or 1984+ LDGT1/2
C  if mileage is above 50,000.
C
   45 IF((IV.EQ.1.AND.MY.LT.1981).OR.
     *  (IV.GE.2.AND.MY.LT.1984).OR.
     *  (BFMILE.LE.5.0)) GOTO 55
C
      DO 50 I23F=1,3
      IB=I23F+1
      D50(I23F)=BFR50K(IB,IGBF,IP,IV)
   50 CONTINUE
C
      BFGT50=BFMILE-5.0
      BFMILE=5.0
C
C  Compute OMCF at FTP and OMCF at Hot FTP.  Hot FTP is split
C  between hot stabilized ( 0.521 ) and hot start ( 0.479 ) operating
C  conditions.
C
   55 OMCFF=ZBF(3)+DBF(3)*BFMILE+D50(3)*BFGT50
C
      OMCFH=0.521*(ZBF(1)+DBF(1)*BFMILE+D50(1)*BFGT50)
     *      +0.479*(ZBF(2)+DBF(2)*BFMILE+D50(2)*BFGT50)
C
C  Compute the ratio of OMCF at Hot FTP to OMCF at FTP.
C
      RATIO=OMCFH/OMCFF
C
C  Convert the uncorrected (FTP conditions) exhaust ER from gpm to
C  gph and then adjust for no cold start component (idle occurs
C  under Hot FTP conditions).
C
      BEFIDL=BER*CHGHRS(IV)*RATIO
C
C  Identify the pointers and lookup & add in the UnTampered IDle
C  OFfset.
C
      IGUTIO=MY-1978
      IF(IGUTIO.LT.1) IGUTIO=1
      IF(IGUTIO.GT.14) IGUTIO=14
      IVVT=IV
      IF(IV.EQ.3) IVVT=2
      BEFIDL=BEFIDL+UTIDOF(IGUTIO,IP,IVVT)
C
C  If LDGV, raise low VMT CO cases (which can go negative by using
C  the curve) to an average (by the low end data) value.
C
      IF(IV.EQ.1.AND.IP.EQ.2.AND.IREJN.EQ.1)
     *    BEFIDL=AMAX1(BEFIDL,COLIM(IGUTIO))
      IF(IV.EQ.1.AND.IP.EQ.2.AND.IREJN.EQ.2)
     *    BEFIDL=AMAX1(BEFIDL,COLIM(IGUTIO)*1.17)
C
C  BEFLEV is a function of BEFIDL, LEVIMP and standards.
C  Save LDGV idle emissions to use for estimating new
C  standard truck emissions.
C
      IF(IV.EQ.1) LDVIDL(IDX) = BEFIDL
C
C  Estimate idle emissions of new standard vehicles based on
C  a ratio of new car standards.
C
      IF(MY.GE.1994.AND.NEWFLG.LT.5)
     *   BEFIDL=BEFIDL*CRT(LMY,IV)+LDVIDL(IDX)*FAC(LMY,IP,IV)
C
      IF(ILEV.EQ.1) THEN
C
         IF(IV.LT.3) BEFLEV=BEFIDL*PTIER1
     *                     +LDVIDL(IDX)*FACL(LMY,IP,IV)
C
CCC      IF(IV.EQ.3)
CCC  *       BEFLEV=BEFIDL*(LEVIMP(3,LMY,4)*LDTAB(2,1)
CCC  *                     +LEVIMP(3,LMY,5)*LDTAB(2,2)
CCC  *                     +LDVIDL(IDX)*FACL(LMY,IP,IV))
CCCCC
      ENDIF
C
C  Apply the I/M Credit then same way to LEV and tier 1
C
      BEFIDL=BEFIDL*CREDIM
      BEFLEV=BEFLEV*CREDIM
C
C   Add in the Idle offset for Tampering
C
      BEFIDL=BEFIDL+OFFIDL(IP,IDX,IV)
C
C   Reduce OFFIDL by TIER1 percent
C
      BEFLEV=BEFLEV+OFFIDL(IP,IDX,IV)*PTIER1
C
C   HC Only Subtract off the Methane( LD gas and MC )
C   Compute idle VOC, TOG, or NMOG, when necessary.
C
      IF(IP.EQ.1) THEN
          LEVMTH = BEFLEV*(MTHIDL(IREJN)/BEFIDL)
          BEFIDL = BEFIDL - MTHIDL(IREJN)
          BEFLEV = BEFLEV - LEVMTH
          IF(NMHFLG.EQ.3) BEFIDL = BEFIDL * VOCFID(MY,IV)
          IF(NMHFLG.EQ.3) BEFLEV = BEFLEV * VOCFID(MY,IV)
          IF(NMHFLG.GE.4) BEFIDL = BEFIDL * TOGFID(MY,IV,OXYFLG)
          IF(NMHFLG.GE.4) BEFLEV = BEFLEV * TOGFID(MY,IV,OXYFLG)
C
C   Add back in the methane for THC and NMOG.
C
          IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.5) THEN
             BEFIDL = BEFIDL + MTHIDL(IREJN)
             BEFLEV = BEFLEV + LEVMTH
          ENDIF
C
      ENDIF
      IF(ILEV.EQ.1) BEFIDL = BEFLEV
C     IF(BEFIDL.LT.0.0) BEFIDL=0.0
C
C   Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFIDL=BEFIDL*(OMTCF(IDX,IP,IV)/OMTCFF(IDX,IP,IV))
C
C   Apply the Idle Speed Correction Factor Adjustment.
C
      BEFIDL=BEFIDL*SCIADJ(IDX,IP,IV)
C
C   Add RVP Correction Factor to gasoline vehicles and trucks
C
      BEFIDL=BEFIDL*RVPCF(IDX,IP,IV)
C
C   Add CO Offset correction to LD gas
C   SCFIDL=> Speed Correction Factor @ 2.5 mph [ BIGSAL ]
C   2.5   => 2.5 mph converts g/mile to g/hour
C   OMTCFF=> OMTCF @ 100% STAB /0/0/0/ 9RVP 75 F
C
      IF(IP.EQ.2) THEN
         IF ((MY.GE.1981.AND.IV.EQ.1) .OR.
     *       (MY.GE.1984.AND.IV.EQ.2) .OR.
     *       (MY.GE.1984.AND.IV.EQ.3))
     *   BEFIDL=BEFIDL+( OFFCO(IDX,IV)*SCFIDL(IDX,IV)*2.5
     *   /OMTCFF(IDX,IP,IV) )
      ENDIF
C
C  Determine and apply the fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
      BEF92=BER
C
      IF(ILEV.EQ.1)
     *CALL BEFUN(1992,IP,IV,VMTAGE,S1,S2,BEF92)
C
      BEFADJ=BEF92
C
      BEFADJ=BEFADJ*CREDIM
      FUELCF=FUEL(ICY,IP,IV,MY,BEFADJ)
      BEFIDL=BEFIDL*FUELCF
C
      RETURN
      END
