C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE REFUEL(ICY)
C
C  REFUEL computes the refueling losses HC in grams / mile by relative
C  (IDX ordered) model year for the 4 vehicle types affected.
C
C  Called by EFCALX.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /CITRV1/ RVPUWX
C    /FLAGS3/ RLFLAG
C    /MAXIMA/ MAXYRS
C    /RLCOM1/ ROADFE,SPILL,OBED,OBES,OBDF
C    /RLCOM2/ IOBMY,IVOB,IS2SY,NPHASE,S2EFF(4)
C    /TAMOU2/ TOB
C    /TEMPS/  AMBT,TEMMIN,TEMMAX,TEMEVP
C
C  Output on return:
C
C    common blocks:
C    /RLCOM3/ RLRATE
C
C  Local array subscripts:
C
C  S2LEFT(4)       -  S2LEFT ( IVTAM )
C
C  Local variable dictionary:
C
C  Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IALLYR   I    all stations year: year by which phase in of Stage II vrs
C                controls is completed => 100% of volume of gas passes through
C                Stage II vrs.
C  IOBCHK   I    check for ob vrs installation, given ob vrs program exists:
C                    1 = no ob vrs    2 = ob vrs installed
C  DFTEMP   R    dispensed fuel temperature (F).
C  MDX      I    model year index into refueling loss rate arrays
C  RVP      R    RVP of the fuel in use that calendar year
C  S2LEFT   R    percentage left of uncontrolled refueling losses after Stage
C                II controls are applied.
C  S2VOL    R    percentage of gasoline pumped that has S2 controls applied.
C  TDFDIF   R    difference in temperature of fuel in vehicle
C                tank and dispensed fuel
C
C  Notes:
C
C  RLFLAG impact:
C
C    1 or (3 & my/iv  -  use the uncontrolled rl rate (assume no recovery
C    not in OB program)  systems installed)
C
C    2 or (4 & my/iv  -  compute and assign Stage II (S2) only vrs controlled
C    not in OB program)
C
C    3 or 4 and my/iv -  weight uncontrolled & controlled rate by best estimate
C    in OB program       of ob vrs tamper rate & its complement, respectively,
C                        to get ob only vrs controlled refueling losses.
C                        Assume 100% installation.
C
C    4                -  compute and assign S2 & ob vrs controlled losses.
C                        See above.
C
C    5                -  assign 0.0 (no losses case: results = release MOBILE3).
C
C  The factors returned are weighted by their corresponding TF's and summed
C  across IDX to get the refueling loss for ICY for each IVTAM.
C
C  REFUEL was changed in MOBILE4.1 to use
C  an equation to compute Displacement.
C
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
C
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /RLCOM1/ ROADFE(32,4),DISPL,SPILL,OBED,OBES,OBDF(4)
      COMMON /RLCOM2/ IOBMY,IVOB(4),IS2SY,NPHASE,S2EFF(4)
      COMMON /RLCOM3/ RLRATE(25,4)
      COMMON /TAMOU2/ TOB(25,4),GCONLY(25,4)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION S2LEFT(4)
C
C  Initialize
C
      DO 10 IVTAM=1,4
      DO 10 IDX=1,MAXYRS
      RLRATE(IDX,IVTAM)=0.0
   10 CONTINUE
C
C  Refueling losses are zero if the no losses case is selected.
C
      IF(RLFLAG.EQ.5) RETURN
C
C  If Stage II vrs requested, compute the corresponding cf by IVTAM.
C
      IF(RLFLAG.EQ.1.OR.RLFLAG.EQ.3) GOTO 20
C
      IALLYR=IS2SY+NPHASE
      IF(ICY.LT.IS2SY) S2VOL=0.0
      IF(ICY.GE.IS2SY.AND.ICY.LT.IALLYR) S2VOL=(ICY-IS2SY)/(NPHASE*1.0)
      IF(ICY.GE.IALLYR) S2VOL=1.0
C
      S2LEFT(1)=1.0-S2EFF(1)*S2VOL
      S2LEFT(2)=S2LEFT(1)
      S2LEFT(3)=S2LEFT(2)
      S2LEFT(4)=1.0-S2EFF(4)*S2VOL
C
C  Set and limit the Dispensed Fuel Temperature (F).
C
   20 DFTEMP=AMBT
      IF(TEMFLG.EQ.1) DFTEMP=(TEMMAX+TEMMIN)/2.
C
      IF(DFTEMP.LT.20.0) DFTEMP=20.0
      IF(DFTEMP.GT.95.0) DFTEMP=95.0
C
C  Set and limit the difference in temperature of the
C  fuel in vehicle tank and dispensed fuel.
C
      TDFDIF=0.418*DFTEMP-16.6
C
      IF(TDFDIF.GT.20.0) TDFDIF=20.0
C
C  Set the displacement RVP.
C
      RVP=RVPUWX
C
C  Calculate and limit the Displacement (g/gal).
C
      DISPL=-5.909-0.0949*TDFDIF+0.0884*DFTEMP+0.485*RVP
C
      IF(DISPL.LT.1.80) DISPL=1.80
C
C  REFUEL is only called once per scenario, computing and saving for later
C  use all IVTAM * IDX onboard rl offsets.  If RLFLAG = 3, the results are
C  program applies.
C
      DO 40 IVTAM=1,4
      DO 30 IDX=1,MAXYRS
C
C  Set the road fuel economy MY pointer.
C
      MY=ICY+IDX-MAXYRS
C
      IF(MY.LE.1969)MDX=1
      IF(MY.GT.1969.AND.MY.LT.2000)MDX=MY-1969+1
      IF(MY.GE.2000)MDX=32
C
C  RLFLAG = 1: Use uncontrolled refueling loss.
C
      IOBCHK=2
      IF(MY.LT.IOBMY.OR.IVOB(IVTAM).EQ.1.OR.IDX.EQ.1) IOBCHK=1
      IF(RLFLAG.EQ.1.OR.(RLFLAG.EQ.3.AND.IOBCHK.EQ.1))
     *    RLRATE(IDX,IVTAM)=(DISPL+SPILL)/ROADFE(MDX,IVTAM)
C
C  RLFLAG = 2: Calculate Stage II only vrs controlled loss.
C
      IF(RLFLAG.EQ.2.OR.(RLFLAG.EQ.4.AND.IOBCHK.EQ.1))
     *  RLRATE(IDX,IVTAM)=(S2LEFT(IVTAM)*DISPL+SPILL)/
     *                     ROADFE(MDX,IVTAM)
C
C  RLFLAG = 3: Calculate an onboard only vrs controlled refueling loss.
C
C  Use the evap cannister only (ID = 6) probability of tampering +
C  100% of the given IV x MY fleet is equipped.  Obtain the onboard loss by
C  weighting the uncontrolled and controlled refueling loss rates by the
C  tampering rate and its complement, respectively.
C
      IF(RLFLAG.GE.3.AND.IOBCHK.EQ.2) RLRATE(IDX,IVTAM)=
     * ((1.0-(1.0-TOB(IDX,IVTAM))*OBED)*DISPL+OBES*SPILL)/
     * ROADFE(MDX,IVTAM)+OBDF(IVTAM)
C
C  RLFLAG = 4: Calculate S2 & ob vrs controlled refueling loss:
C
C  Handled above.  Reduces to either RLFLAG=2 or RLFLAG=3 case, depending on
C  IOBCHK.  Onboard takes precedent over Stage II and you can't get credit for
C  both.
C
      CONTINUE
C
   30 CONTINUE
   40 CONTINUE
C
      RETURN
      END
