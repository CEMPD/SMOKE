C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION CALUDI(RVPW,TMIN,TMAX,FILLED)
C
C  CALUDI uses passed in fuel RVP levels, minimum and maximum fuel tank
C  temperatures and fleet average percent of fuel tank filled to estimate
C  an Uncontrolled Diurnal emission rate.
C
C  Called by LOCAL
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list:
C    RVPW,TMIN,TMAX,FILLED
C
C    common blocks:
C    /REGION/ IREJN
C
C  Output on return:
C
C    function: CALUDI
C
C  Local array subscripts :
C
C  AIRPRE(2)  -  AIRPRE ( IREJN )
C  FUELT(2)   -  FUELT  ( JV )
C  VAPOR(2)   -  VAPOR  ( JV )
C
C  Local variable / array dictionary:
C
C  Name   Type              Description
C  ------ ----  ----------------------------------------------------
C  A       R    coefficient used to calculate VAPOR(I)
C  AIRPRE  R    air pressure
C  A100    R    coefficient used to calculate A
C  C       R    constant, used to calculate A100
C  DENSTY  R    fuel density, a function of RVPW
C  FILLED  R    percent of fuel tank filled
C  FUELT   R    incremental fuel tank temperature in degrees Fahrenheit
C  GTP     R    grams of HC loss for temperature pair FUELT
C  PI      R    Greek pi = ratio of circumference of a circle to its diameter.
C               Used in computing X.
C  RVPW    R    fuel RVPs used to compute the components of the Wade Index
C  TMAX    R    maximum temperature
C  TMIN    R    minimum temperature
C  TP1     R    first coefficient for Wade equation
C  TP2     R    second coefficient for Wade equation
C  TP3     R    third coefficient for Wade equation
C  UDISUM  R    incremental and, eventually, total Uncontrolled Diurnal rate
C  VAPOR   R    vapor pressure at FUELT (2 temperatures, 1 degree apart)
C  VP100   R    vapor pressure at 100F, a function of RVPW
C  VSPACE  R    vapor space in cubic feet, a function of FILLED
C  WMOLEC  R    molecular weight, a function of RVPW and fuel temp
C  X       R    coefficient used to calculate A100, a function of VP100
C
C  Notes:
C
C  CALUDI is the sum of HC loss at each fuel temperature increment under the
C  input conditions.
C  CALUDI temperature range "cuts" code changed for MOBILE4.1 so as to
C  not lose the last difference pair due to round off.  Also added an
C  immediate return if TMIN=TMAX.
C
      INCLUDE 'REGION.I'
C
      DIMENSION AIRPRE(2),FUELT(2),VAPOR(2)
C
      DATA AIRPRE/14.696,12.5/,PI/3.14159/
C
      UDISUM=0.0
C
      IF(TMIN.EQ.TMAX) GOTO 99
C
C  To calculate the CALUDI value, first compute several parameters and then
C  move stepwise through the fuel tank temperature range, adding each
C  degree difference pair's contribution to the sum UDI of the grams HC
C  loss over the entire range.
C
C  Calculate fuel density for given RVPW.
C
      DENSTY=6.4-0.01977*RVPW
C
C   Calculate vapor space under given percent of fuel tank filled.
C
      VSPACE=2.4062-0.02139*FILLED
C
C  Calculate vapor pressure at 100 (VP100) for given RVPW.
C
      VP100=1.0223*RVPW
     *    +(0.0357*RVPW)/(1.0-0.0368*RVPW)
C
C  Calculate A100 according to VP100.
C
      IF(VP100.LT.14.18) GOTO 10
      C=80.861
      X=0.11*COS((4.0*VP100-9.0)*PI/14.0)
     *  +5.4*ALOG(VP100)
      GOTO 20
C
   10 C=66.561
      X=0.12*COS((VP100-6.0)*PI/4.0)
     * -0.21*SIN(2.0*PI/7.5*(VP100-4.0))
C
   20 A100=C-12.822*VP100
     *      +1.3291*VP100**2
     *      -0.07991*VP100**3
     *      +1.9017E-03*VP100**4-X
C
C  Initialize fuel tank temperature pair.
C
      FUELT(1)=TMIN
      FUELT(2)=FUELT(1)+1.0
      IF(FUELT(2).GT.TMAX) FUELT(2)=TMAX
C
C  Iteration starts here.
C
   30 CONTINUE
C
C  Calculate molecular weight.
C
      WMOLEC=69.69-1.274*RVPW
     *         +0.059*(FUELT(1)+FUELT(2))/2.0
C
C  Calculate vapor pressures.
C
      DO 40 JV=1,2
      A=A100+(100.0-FUELT(JV))
     *      *((262.0/(A100/6.0+560.0))-0.01328)
C
C  pass JV's A  < 0.0  =>  CALUDI < 0.0  =>  diurnal evap < 0.0  =>  fatal error
C  Technically, high RVPW (in A100) and high temperature (in FUELT(JV)) does not
C  make sense anyway: the gas tank would blow up.
C
      IF(A.LT.0.0) CALL QUITER(A,JV,97,INERR)
C
      VAPOR(JV)=14.696
     *        -0.53059*A
     *        +7.6961E-03*A**2
     *        -5.4907E-05*A**3
     *        +1.7044E-07*A**4
   40 CONTINUE
C
C  Apply Wade equation.
C
      TP1=VSPACE*118040.0*DENSTY/(690.0-4.0*WMOLEC)
      TP2=VAPOR(1)/(AIRPRE(IREJN)-VAPOR(1))
     *   +VAPOR(2)/(AIRPRE(IREJN)-VAPOR(2))
      TP3=(AIRPRE(IREJN)-VAPOR(1))/(FUELT(1)+460.0)
     *   -(AIRPRE(IREJN)-VAPOR(2))/(FUELT(2)+460.0)
      GTP=TP1*TP2*TP3
C
      UDISUM=UDISUM+GTP
C
      FUELT(1)=FUELT(1)+1.0
      FUELT(2)=FUELT(2)+1.0
      IF(FUELT(2).GT.TMAX) FUELT(2)=TMAX
      IF(FUELT(1).LT.TMAX) GOTO 30
C
   99 CALUDI=UDISUM
C
      RETURN
      END
