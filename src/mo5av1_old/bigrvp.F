C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGRVP(ICY,IP,IG)
C
C  BIGRVP computes the combined temperature and RVP correction factor algorithm
C  used when the pollutant's temperature >= 75 degrees (FTP standard) and
C  the RVP cf alone when T < 75 degrees.
C
C  See RVPEXH for an overall description of the RVP cf cases.
C
C  Called by BIGTCF.
C
C  Calls AMIN1 AND AMAX1.
C
C  Input on call:
C
C    parameter list: IP,IG
C    common blocks:
C    /CITRV2/ RVPX,RVP090
C    /RVPEX1/ RG23,TTL
C    /TEMPC3/ MDLOHI
C    /TEMPS/  TEMEXH
C
C  Output on return:
C
C    common blocks:
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CFB      R    RVP cf at begin point of temperature interpolation interval
C  CFL      R    RVP cf at end   point of temperature interpolation interval
C  RDIFA    R    the larger of 0.0 and (applicable RVP - 9.0)
C  REAT     R    adjusted RVP exhaust IP temperature
C  REUT     R    RVP exhaust IP unadjusted temperature
C  TDIFA    R    the larger of 0.0 & (adjusted T capped at 95.0) - 75.0)
C  TDIFU    R    (unadjusted exhaust T capped at 110.0) - 75.0
C
C  Notes:
C
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER RFGSYR
      INTEGER RFGFLG,WINFLG
      CHARACTER*1 ASTMCL
      LOGICAL RFGON
C
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
      COMMON /RVPEX1/ R9G1(3,2),RG23(3,3,3,2),TLL(2)
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
C  Reduce supplied RVP to 11.7 for the RVP cf calculation only.
C
      RVP1=AMIN1(RVPX,11.7)
C
C  If RVP x T temp is .GT. 95 degrees, then reduce it to 95 for here only.
C  If T term temp is .GT. 110 degrees, then reduce it to 110 for here only.
C
      REAT=AMIN1(TEMEXH(IP),95.0)
      TDIFA=REAT-75.0
      TDIFA=AMAX1(TDIFA,0.0)
      REUT=AMIN1(TEMEXH(IP),110.0)
      TDIFU=REUT-75.0
C
C  Identify RVP to be used. RVP can drop to 7.0 if >75 degrees, else 9.0
C
      RIVAL=RVP1
C
C  Maintain zero NOx effect for oxygenated and reformulated gasoline
C
      IF(IP.EQ.3) THEN
         INEW=1
         IF(OXYFLG.EQ.2.OR.RFGON.AND.WINFLG.EQ.2) THEN
           IF(ICY.LT.IUSESY) RIVAL=RVPBAS
           IF(ICY.GE.IUSESY) RIVAL=RVPIUS
           INEW=2
         ENDIF
         IF(RFGON.AND.WINFLG.EQ.1) THEN
           RIVAL=8.7
           INEW=2
         ENDIF
C
C  Apply weathering to the NOx RVP. Equation uses the calculated Hot
C  Soak temperature. Check that this weathered RVP is within limits.
C
         RVPLSS=-2.4908
     &          +0.026196*TEMEVP(1)
     &          +0.00076898*RIVAL*TEMEVP(1)
            IF(RVPLSS.LT.0.0) RVPLSS=0.0
         IF(INEW.EQ.2) RIVAL=RIVAL-RVPLSS
      ENDIF
C
C  Determine RVP difference from 9.0, set limits.
C
      RDIFA=RIVAL-RVP090
      IF(REAT.GT.75.0) RDIFA=AMAX1(RDIFA,-2.0)
      IF(REAT.LE.75.0) RDIFA=AMAX1(RDIFA,0.0)
C
C  Compute the cf for each bag.
C
      IGR=IG-1
C
C  Case 1: one equation produces a combined RVP & high/neutral T cf.
C
      IF(MDLOHI(IP).EQ.1) GOTO 20
C
C  Evaluate the combined RVP & high T cf formula.
C  Except for NOx, combined RVP & T cf is never less than 1.0.
C
      DO 10 IB=1,3
      TCF(IB)=RG23(1,IB,IP,IGR)*RDIFA+
     *               RG23(2,IB,IP,IGR)*TDIFU+
     *               RG23(3,IB,IP,IGR)*RDIFA*TDIFA
      IF(IP.NE.3) TCF(IB)=AMAX1(TCF(IB),0.0)
      TCF(IB)=EXP(TCF(IB))
   10 CONTINUE
      RETURN
C
C  Case 2: interpolation produces a RVP cf for temperatures from 45-75 degrees.
C
C      (a) First assign the endpoints at 45 and 75 degrees.
C
   20 CFB=1.0
      DO 30 IB=1,3
      CFL=RG23(1,IB,IP,IGR)*RDIFA
      CFL=EXP(CFL)
C
C      (b) Now interpolate between the 2 endpoints (CFB,CFL).
C
      TCF(IB)=CFB+(CFL-CFB)*(REAT-TLL(1))/TLL(2)
C
   30 CONTINUE
C
      RETURN
      END
