C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CALTEM
C
C  CALTEM calculates exhaust and evaporative temperatures.
C
C  Called by LOCAL.
C
C  Calls OUTPOL and QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS3/ TEMFLG
C    /FLAGS4/ PRTFLG
C    /MAXIMA/ MAXPOL
C    /TEMPS/  AMBT,TEMMIN,TEMMAX
C
C  Output on return:
C
C    common block:
C    /TEMPS/ TEMEXH,TEMEVP
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ITEM     I    QUITER temperature name number
C  TRISE    R    temperature rise term common to many equations
C
C  Notes:
C
C  None.
C
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
C  Determine pollutants.
C
      CALL OUTPOL(IP1,IP2)
C
C  Compute the temperature rise term.
C  Warn user that diurnal evap efs may be inaccurate if TRISE>40F.
C  Skip warning if CO-only run.
C
      TRISE=TEMMAX-TEMMIN
C
      IF(TRISE.GT.40.0.AND.PRTFLG.NE.2) CALL QUITER(TRISE,0,98,INFORM)
C
C  Compute exhaust temperature for each pollutant.
C
      DO 30 IP=IP1,IP2
C
C  Select proper equation.
C
      GOTO (10,15,20),IP
C
C  Calculate HC exhaust temperature if TEMFLG=1, otherwise
C  use the ambient temperature.
C
   10 TEMEXH(IP)=AMBT
C
      IF(TEMFLG.EQ.1)
     *   TEMEXH(IP)=2.2857+0.97674*TEMMIN+0.56881*TRISE
     *             +0.0024642*TEMMIN*TRISE
C
      GOTO 25
C
C  Calculate CO exhaust temperature if TEMFLG=1, otherwise
C  use the ambient temperature.
C
   15 TEMEXH(IP)=AMBT
C
      IF(TEMFLG.EQ.1)
     *   TEMEXH(IP)=2.2857+0.97674*TEMMIN+0.56881*TRISE
     *             +0.0024642*TEMMIN*TRISE
C
      GOTO 25
C
C  Calculate NOx exhaust temperature if TEMFLG=1, otherwise
C  use the ambient temperature.
C
   20 TEMEXH(IP)=AMBT
C
      IF(TEMFLG.EQ.1)
     *   TEMEXH(IP)=2.2857+0.97674*TEMMIN+0.56881*TRISE
     *             +0.0024642*TEMMIN*TRISE
C
C  Round the exhaust temperature to the nearest 0.1F.
C
   25 TEMEXH(IP)=FLOAT(IFIX(TEMEXH(IP)*10.+.5))/10.
C
C  Check to see if TEMMIN<=TEMEXH<=TEMMAX.
C
      IF(TEMEXH(IP).LT.TEMMIN.OR.TEMEXH(IP).GT.TEMMAX)
     *   CALL QUITER(TEMEXH(IP),2,111,INERR)
C
   30 CONTINUE
C
C  Compute evaporative temperature for each evap case.
C
      DO 90 IHRD=1,6
C
C  Select proper equation.
C
      GOTO (45,50,55,60,65,70),IHRD
C
C  Calculate the hot soak temperature if TEMFLG=1 and TEMMAX is
C  at or above 40F, otherwise use the ambient temperature.
C
   45 TEMEVP(IHRD)=AMBT
C
      IF(TEMFLG.EQ.1.AND.TEMMAX.GE.40.)
     *   TEMEVP(IHRD)=-1.7474+1.0290*TEMMIN+0.99202*TRISE
     *               -0.0025173*TEMMIN*TRISE
C
      GOTO 80
C
C  Calculate running loss temperature if TEMFLG=1 and TEMMAX is
C  at or above 40F, otherwise use the ambient temperature.
C
   50 TEMEVP(IHRD)=AMBT
C
      IF(TEMFLG.EQ.1.AND.TEMMAX.GE.40.)
     *   TEMEVP(IHRD)=-1.1977+1.0205*TEMMIN+1.0181*TRISE
     *               -0.0023797*TEMMIN*TRISE
C
      GOTO 80
C
C  Calculate 8 AM (10AM) diurnal tank temperature.
C  Tank temperature lags ambient temperature by about 2 hours.
C
   55 TEMEVP(IHRD)=TEMMIN+0.53*TRISE
C
      GOTO 80
C
C  Calculate 10-11 AM (12 noon) diurnal tank temperature.
C  Tank temperature lags ambient temperature by about 2 hours.
C
   60 TEMEVP(IHRD)=TEMMIN+0.71*TRISE
C
      GOTO 80
C
C  Calculate 2 PM diurnal tank temperature.
C  Tank temperature lags ambient temperature by about 2 hours.
C
   65 TEMEVP(IHRD)=TEMMIN+0.94*TRISE
C
      GOTO 80
C
C  Calculate resting loss temperature if TEMFLG=1 and TEMMAX is
C  at or above 40F, otherwise use the ambient temperature.
C
   70 TEMEVP(IHRD)=AMBT
      IF(TEMFLG.EQ.1.AND.TEMMAX.GE.40.)
     *    TEMEVP(IHRD)=-2.4187+0.99031*TEMMIN+0.67839*TRISE
C
C  Round the evaporative temperature to the nearest 0.1F.
C
   80 TEMEVP(IHRD)=FLOAT(IFIX(TEMEVP(IHRD)*10.+.5))/10.
C
C  Check to see if TEMMAX<=TEMEVP<=TEMMAX.
C
      ITEM=2+IHRD
C
      IF(TEMEVP(IHRD).LT.TEMMIN.OR.TEMEVP(IHRD).GT.TEMMAX)
     *   CALL QUITER(TEMEVP(IHRD),ITEM,111,INERR)
C
   90 CONTINUE
C
C  If any evap temperature is 40F or less, or if minimum temp is 25F
C  or less, warn the user that no evap calculations will occur.
C  Skip warning if CO-only run.
C
      IF((TEMMAX   .LE.40.0 .OR.
     *    AMBT     .LE.40.0.AND.TEMFLG.EQ.2 .OR.
     *    TEMEVP(1).LE.40.0 .OR.
     *    TEMEVP(2).LE.40.0 .OR.
     *    TEMMIN   .LE.25.0).AND.
     *    PRTFLG   .NE.2) CALL QUITER(0.,0,83,INFORM)
C
      RETURN
      END
