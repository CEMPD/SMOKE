C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE RVPEXH(ICY)
C
C  MOBILE4 calculates a correction factor array to adjust the base emission
C  factor for fuel volatility (= Reid Vapor Pressure = RVP).  This
C  multiplicative cf is neutral (= 1.0) only when (1) the actual RVP <= 9.0
C  (= the RVP of the fuel used in the tests upon which the ber slopes and
C  intercepts in Block Data are based) and / or (2) the ambient temperature
C  (TEMEXH) is <= 45.0 (below 45 degrees there is no volatility effect) & / or
C  the vehicle class is not equipped in the given MY with the technology known
C  by testing to be subject to an RVP exhaust effect on the given pollutant.
C  Thus the effect depends on IP and MY as well as RVP and TEMEXH.
C
C  There are 3 nonneutral cases.  RVPEXH computes the MYG 1 HC & CO case,
C  as well as initializing the cf array to 1.0 to cover the neutral cases.
C  BIGTCF computes the 2 operating mode ("bag") dependent cases: MYGs 2 & 3
C  >= 75 degrees and MYGs 2 & 3 < 75 degrees.
C
C  Called by EFCALX.
C
C  Calls AMIN1, IGRVPT and OUTPOL.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /CITRV2/ RVPX,RVP090
C    /MAXIMA/ MAXPOL,MAXYRS
C    /RVPEX1/ R9G1,TLL
C    /TEMPS/  TEMEXH
C
C  Output on return:
C
C    common blocks:
C    /RVPEX2/ RVPCF
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CFB      R    RVP cf at begin point of temperature interpolation interval
C  CFL      R    RVP cf at end   point of temperature interpolation interval
C  RADJCF   R    RVP cf assessed at FTP (75.0) temperature - accounts only for
C                RVP deviation from norm, not temperature.
C  REAT     R    RVP exhaust IP temperature
C
C  Notes:
C
C  Support for oxygenated and cert. fuels from MOBILE4 was removed in 4.1.
C
C
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /RVPEX1/ R9G1(3,2),RG23(3,3,3,2),TLL(2)
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
C  Determine pollutants.
C
      CALL OUTPOL(IP1,IP2)
C
C  Initialize the RVP exhaust correction factor to the neutral 1.0 and depart,
C  if no IVTAM x IP x IDX cases will be in the "not-by-bag" computed category.
C
      DO 10 IVTAM=1,4
      DO 10 IP=IP1,IP2
      DO 10 IDX=1,MAXYRS
      RVPCF(IDX,IP,IVTAM)=1.0
   10 CONTINUE
C
C  The RVPCF algorithm applies only to HC & CO (not NOx) and MYG 1 (1970-79
C  for LDGV, 1970-80 for LDGT1, 1978-80 for LDGT2, no MYG 1 for HDGV).
C
      IF((IP1.EQ.3.AND.IP2.EQ.3).OR.ICY-MAXYRS.GT.1980) RETURN
C
C  Reduce supplied RVP to 11.7 for the RVP cf calculation only.
C
      RVP1=AMIN1(RVPX,11.7)
C
C  Now can identify and compute the "not-by-bag" cases for LDGV/T.
C
      DO 30 IVTAM=1,3
      DO 30 IP=IP1,IP2
      IF(IP.EQ.3) GOTO 30
C
C  If IP's temperature is <= 45 degrees, then skip this IP iteration.
C  If it is above 95 degrees, then reduce it to 95 for this computation.
C
      IF(TEMEXH(IP).LE.45.0) GOTO 30
      REAT=AMIN1(TEMEXH(IP),95.0)
C
      DO 20 IDX=1,MAXYRS
      MY=ICY+IDX-MAXYRS
      IG=IGRVPT(MY,IVTAM)
C
C  If MYG 0, then no computation is needed for this MY.  Do next IDX iteration.
C  If MYG 2 has been reached, further calculations for this IVTAM x IP case
C  will be done in BIGTCF, so exit IDX loop.
C
      IF(IG.EQ.0) GOTO 20
      IF(IG.GT.1) GOTO 30
C
C  Use base or in-use RVP.  If the selected RVP is <= 9.0, skip the
C  computation.  Once RIVAL <= 9.0, exit IDX loop, not just 1 iteration.
C
      RIVAL=RVP1
      IF(RIVAL.LE.RVP090) GOTO 30
C
C  Compute the cf.  There are 2 cases, regarding temperature.
C
C  (1) Temperature is neutral.
C
      RADJCF=(R9G1(1,IP)+R9G1(2,IP)*RIVAL)/R9G1(3,IP)
      IF(REAT.EQ.75.0) RVPCF(IDX,IP,IVTAM)=RADJCF
      IF(REAT.EQ.75.0) GOTO 20
C
C  (2) Temperature is a variable in the RVP cf calculation.  Interpolate
C      between 2 cf's computed at temperatures for which data exists.
C
C      (a) First assign the endpoints at 45 and 75 degrees.
C
      CFB=1.0
      CFL=RADJCF
C
C      (b) Now interpolate between the 2 endpoints (CFB,CFL) or extrapolate
C          beyond CFL (up to 95), using actual (REAT) and (TLL) temperatures.
C
      RVPCF(IDX,IP,IVTAM)=
     *    CFB+(CFL-CFB)*(REAT-TLL(1))/TLL(2)
C
   20 CONTINUE
   30 CONTINUE
C
      RETURN
      END
