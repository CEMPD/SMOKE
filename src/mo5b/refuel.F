C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE REFUEL(ICY,IY)
C
C  REFUEL computes the refueling losses HC in grams / mile by relative
C  (IDX ordered) model year for the 4 vehicle types affected.
C
C  Called by EFCALX.
C
C  Input on call:
C
C    parameter list: ICY, IY
C    common blocks:
C    /CITRV1/ RVPUWX
C    /FLAGS3/ RLFLAG
C    /MAXIMA/ MAXYRS
C    /RLCOM1/ ROADFE,SPILL,OBED,OBES,OBDF,OBPIP,OBPIY
C    /RLCOM2/ IVOB,IS2SY,NPHASE,S2EFF(4)
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
C  IDXOB    I    Index to the ob vrs phase-in tables
C
C  Notes:
C
C  RLFLAG impact:
C
C    1 or (3 & my/iv  -  maps to RLFLAG = 3
C    not in OB program)
C
C    2 or (4 & my/iv  -  maps to RLFLAG = 4
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
C  5-May-1994 @ CSC-tmm Request 2-411
C    Code added to implement 1994 Onboard Refueling Vapor Recovery Regulations
C
C  21-Nov-1994 @ CSC-ked Request 0-411
C    Allow printing of Warning message 154 only once for July runs
C  09-June-1995 @ CSC-wcs Request # 0-411 ORVR Fix (Allow printing of warning
C    message 156 only once for July runs.)
C
C  08-Feb-1996 @ DynTel Request 2-615 ORVR effects are disabled when NEWFLG 
C                is set to 5 or 6.
C
      INCLUDE 'CITRV1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'RLCOM1.I'
      INCLUDE 'RLCOM2.I'
      INCLUDE 'RLCOM3.I'
      INCLUDE 'TAMOU2.I'
      INCLUDE 'TEMPS.I'
C
      DIMENSION S2LEFT(4)
C
C  Initialize
C
      S2LEFT(1) = 1.0
      S2LEFT(2) = 1.0
      S2LEFT(3) = 1.0
      S2LEFT(4) = 1.0
C
      DO 10 IVTAM=1,4
      DO 10 IDX=1,MAXYRS
        RLRATE(IDX,IVTAM)=0.0
   10 CONTINUE
C
C  03/25/94 TMM@CSC REQUEST #2-411
C    Turn on onboard vapor recovery unconditionally
C      and implement the phase-in scheme in OBPIY and OBPIP
C
C    RLFLAG value 1 now maps to 3; 2 maps to 4;
C
C  Refueling losses are zero if the no losses case is selected.
C
      IF(RLFLAG.EQ.5) RETURN
C
C   If Clean Air Act is disabled (NEWFLG=5 or 6),
C     disable Onboard Vapor Recovery.
C
      IF (NEWFLG.GE.5) THEN
C
C  09-June-1995 @ CSC-wcs Request # 0-411 ORVR Fix (Allow printing of warning
C    message 156 only once for July runs.)
C
        IF (IY.EQ.1) THEN
          CALL QUITER(0.,0,156,INERR)
        END IF
C
        IF (RLFLAG.EQ.3) THEN
          RLFLAG = 1
        ELSEIF (RLFLAG.EQ.4) THEN 
          RLFLAG = 2
        ENDIF
C
      ELSE
C
        IF (RLFLAG.EQ.1) THEN
          RLFLAG = 3
        ELSEIF (RLFLAG.EQ.2) THEN
          RLFLAG = 4
        ENDIF
C
      ENDIF
C
C  Write warning for ob vrs regulations (1994)
C
      IF (OUTFMT.GE.3.AND.OUTFMT.LE.5.AND.NEWFLG.LT.5.AND.IY.EQ.1)
     *    CALL QUITER(0.,0,154,INERR)
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
C
      IDXOB = 0
      IOBCHK = 1
      IF (IVTAM.LE.3) THEN
        DO 50 I=1,3
          IF (MY.GE.OBPIY(I,IVTAM)) THEN
            IDXOB=I
            IOBCHK=2
          ENDIF
 50     CONTINUE
      ENDIF
C
      IF(RLFLAG.EQ.1.OR.(RLFLAG.EQ.3.AND.IOBCHK.EQ.1)) THEN
          RLRATE(IDX,IVTAM)=(DISPL+SPILL)/ROADFE(MDX,IVTAM)
C
C  RLFLAG = 2: Calculate Stage II only vrs controlled loss.
C
      ELSE IF((RLFLAG.EQ.2.OR.(RLFLAG.EQ.4.AND.IOBCHK.EQ.1))) THEN
        RLRATE(IDX,IVTAM)=(S2LEFT(IVTAM)*DISPL+SPILL)/
     *                     ROADFE(MDX,IVTAM)
C
C  RLFLAG = 3: Calculate an onboard only vrs controlled refueling loss.
C
C  Use the evap cannister only (ID = 6) probability of tampering +
C  100% of the given IV x MY fleet is equipped.  Obtain the onboard loss by
C  weighting the uncontrolled and controlled refueling loss rates by the
C  tampering rate and its complement, respectively.
C
C  Weight the calculation by the phase-in of onboard refueling vapor
C  recovery regulations of 1994.
C
      ELSE IF(RLFLAG.GE.3.AND.IOBCHK.EQ.2) THEN
       RLRATE(IDX,IVTAM)=
     *  (DISPL/ROADFE(MDX,IVTAM))*((1-OBED)*(1-TOB(IDX,IVTAM))
     *                                     *OBPIP(IDXOB,IVTAM)
     *                 +TOB(IDX,IVTAM)*OBPIP(IDXOB,IVTAM)*S2LEFT(IVTAM)
     *                 +(1-OBPIP(IDXOB,IVTAM))*S2LEFT(IVTAM) )
     *      +(SPILL/ROADFE(MDX,IVTAM))*((1-OBES)*OBPIP(IDXOB,IVTAM)
     *                             +(1-OBPIP(IDXOB,IVTAM)))
      ENDIF
C
C  RLFLAG = 4: Calculate S2 & ob vrs controlled refueling loss:
C
C  Handled above.  Reduces to either RLFLAG=2 or RLFLAG=3 case, depending on
C  IOBCHK.  Onboard takes precedent over Stage II and you can't get credit for
C  both.
C
C
   30 CONTINUE
   40 CONTINUE
C
      RETURN
      END
