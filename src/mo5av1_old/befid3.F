C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION BEFID3(MY,IDX,VMTAGE)
C
C  BEFID3 returns the basic idle emission factor adjusted for I/M,
C  methane, tampering and oxy fuel.
C
C  Called by IDLCAL.
C
C  Calls BEFUN, IBFPTR, IDLPTR, PCLEFT, TOGFID, VOCFID and FUEL.
C
C  Input on call:
C
C    parameter list: MY,IDX,VMTAGE
C    common blocks:
C    /FLAGS1/ OXYFLG
C    /FLAGS4/ NMHFLG
C    /IDLEQ1/ ZMLIDL,DRIDL
C    /MAXIMA/ MAXYRS
C    /OMTCOM/ OMTCF,OMTCFF
C    /OPMOD1/ BFSTEP,BFR50K
C    /OPMOD2/ MYSTEP
C    /REGION/ IREJN
C    /RVPEX2/ RVPCF
C    /SPEED6/ SCIADJ
C    /TAMID3/ OFFIDL
C
C  Output on return:
C
C    function: BEFID3
C
C  Local array subscripts:
C
C  D50   (3)  -  D50    ( I23F )
C  DBF   (3)  -  DBF    ( I23F )
C  MTHIDL(2)  -  MTHIDL ( IREJN )
C  ZBF   (3)  -  ZBF    ( I23F )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BER      R    untampered basic emission rate
C  BFGT50   R    mileage greater than 50,000 (BFMILE=6.0 -> BFGT50=1.0)
C  BFMILE   R    cumulative mileage in ICY for IV from MY
C                scaled by 10**-4
C  CREDIM   R    I/M adjustment to idle emission factor
C  D50      R    BF deterioration rate for beyond 50K miles
C  DBF      R    BF deterioration rate for the first 50K miles
C  I23F     I    index for Bag 2 ( = 1 ), Bag 3 ( = 2 ), and
C                FTP ( = 3 )
C  IB       I    index for Bag 1 ( = 1 ), Bag 2 ( = 2 ), Bag 3 ( = 3 ),
C                Bag 3 ( = 3 ), and FTP ( = 4 )
C  MTHIDL   R    methane offset to idle EF in gph (assumes methane is
C                approx. 6% of total HC over the affected MY range for
C                Bags 2 & 3)
C  OMCFF    R    operating mode correction factor at FTP
C                conditions
C  OMCFH    R    operating mode correction factor at Hot FTP
C                conditions
C  RATIO    R    ratio of OMCF at Hot FTP to OMCF at FTP
C  ZBF      R    BF zero mile level
C
C
C  Notes:
C
C     BEFID3 is called by IDLCAL and calculates idle emissions
C     for HC HDGV.  BEFID3 was created for MOBILE5.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OPMOD1/ BFRCOF(8,8,3,6),BFSTEP(8,12,3,3),BFR50K(4,12,3,3)
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /TAMID3/ OFFIDL(3,25,4)
C
      DIMENSION DBF(3),ZBF(3),D50(3)
      REAL MTHIDL(2)
      DATA MTHIDL/0.9,1.5/
      BER=0.0
      BEFID3=0.0
      RATIO=1.0
C
C  Assign ICY and determine the I/M credit.
C
      ICY=MY-IDX+MAXYRS
      CREDIM=PCLEFT(MY,ICY,1,4)
C
      IF(MY.LT.1984) GOTO 50
C
C  Obtain an untampered Basic Emission Rate
C  HDGVs use LDGT2 BER for all HC 1984+
C
      CALL BEFUN(MY,1,3,VMTAGE,S1,S2,BER)
C
C  Set up various bag fractions pointers and mileage values.
C
      MYSTEP=2
      IGBF=IBFPTR(MY,1,3)
      JDX=MAXYRS-IDX+1
C
      DO 10 I23F=1,3
        D50(I23F)=0.0
   10 CONTINUE
      BFMILE=CUMMIL(JDX,4)*.0001
      BFGT50=0.0
C
C  Determine Bag Fraction coefficients for Bag 2, Bag 3, and FTP.
C  Use appropriate coefficients for 1984+ HDGV.
C
      DO 20 I23F=1,3
        IB=I23F*2+2
        ZBF(I23F)=BFSTEP(IB-1,IGBF,1,3)
        DBF(I23F)=BFSTEP(IB,IGBF,1,3)
   20 CONTINUE
C
C  Use kinked Bag Fractions for 1984+ HDGV if mileage is above
C  50,000.
C
      IF(BFMILE.LE.5.0) GOTO 40
      DO 30 I23F=1,3
        IB=I23F+1
        D50(I23F)=BFR50K(IB,IGBF,1,3)
   30 CONTINUE
      BFGT50=BFMILE-5.0
      BFMILE=5.0
C
C  Compute OMCF at FTP and OMCF at Hot FTP.  Hot FTP is split
C  between hot stabilized (.521) and hot start (.479) operating
C  conditions.
C
   40 OMCFH=0.521*(ZBF(1)+DBF(1)*BFMILE+D50(1)*BFGT50)
     *      +0.479*(ZBF(2)+DBF(2)*BFMILE+D50(2)*BFGT50)
      OMCFF=ZBF(3)+DBF(3)*BFMILE+D50(3)*BFGT50
C
C  Compute the ratio of OMCF at Hot FTP to OMCF at FTP.
C
      RATIO=OMCFH/OMCFF
C
C  Convert the uncorrected (FTP conditions) exhaust from gpm to
C  gph and then adjust for no cold start component (idle occurs
C  under Hot FTP conditions).  CHGHRS=20.0
C
      BEFID3=BER*20.0*RATIO
      IF(MY.GE.1986) BEFID3=BEFID3+1.28
      GOTO 60
C
C  In the pre-1984 model year case use the MOBILE4 revised idle EF
C  calculation algorithm.
C
   50 IGIDL=IDLPTR(MY,1,4)
      BEFID3=ZMLIDL(IGIDL,1,4,IREJN)+DRIDL(IGIDL,1,4,IREJN)*VMTAGE
C
C  Apply the I/M credit.
C
   60 BEFID3=BEFID3*CREDIM
C
C  Add in the Idle offset for Tampering
C
      BEFID3=BEFID3+OFFIDL(1,IDX,4)
C
C  Subtract off the Methane
C  Compute idle VOC, TOG, or NMOG, when necessary.
C
      BEFID3=BEFID3-MTHIDL(IREJN)
C
      IF(NMHFLG.EQ.3) BEFID3=BEFID3*VOCFID(MY,4)
      IF(NMHFLG.GE.4) BEFID3=BEFID3*TOGFID(MY,4,OXYFLG)
C
C  Add back in the methane for THC and NMOG.
C
      IF(NMHFLG.EQ.1.OR.NMHFLG.EQ.5) BEFID3=BEFID3+MTHIDL(IREJN)
C     IF(BEFID3.LT.0.0) BEFID3=0.0
C
C  Mult the OMTCFF ( Operating Mode & Temperature CF )
C
      BEFID3=BEFID3*(OMTCF(IDX,1,4)/OMTCFF(IDX,1,4))
C
C  Apply the Idle Speed Correction Factor Adjustment.
C
      BEFID3=BEFID3*SCIADJ(IDX,1,4)
C
C  Add RVP Correction Factor
C
      BEFID3=BEFID3*RVPCF(IDX,1,4)
C
C  Determine and apply the oxy fuel correction factor.
C  Use BEFUN to get full (non-idle) ef-based OXYCF.
C
      BER=BER*CREDIM
      FUELCF=FUEL(ICY,1,4,MY,BER)
      BEFID3=BEFID3*FUELCF
      RETURN
      END
