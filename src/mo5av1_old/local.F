C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE LOCAL(ICY,INERR,*)
C
C  LOCAL initializes the local (city specific) data common blocks.
C
C  Called by PARSEC.
C
C  Calls CALUDI, CALTEM, QUITER and RDCITY.
C
C  Input after call to RDCITY:
C
C    parameter list: ICY
C    common blocks:
C    /CITRV1/ RVPBAS,RVPIUS,IUSESY
C    /CITRV2/ RVPHS,RVP090
C    /FLAGS1/ OXYFLG
C    /FLAGS3/ LOCFLG
C    /OXY1/ SHRMKT,IGASHW,OXYCNT
C    /OXY2/ NFUEL
C    /RFORM1/ RFGFLG,ASTMCL
C    /RFORM2/ RFGSYR,RFGRVP
C    /RVPNAT/ RVPLIM,RVPDIU,PFUL
C    /TEMPS/  AMBT,TEMMIN,TEMMAX,TEMEVP
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /CITCIN/ UDI
C    /CITRV1/ RVPICY
C    /CITRV2/ RVPHS
C    /RFORM1/ IASTM,IPHASE,WINFLG,RFGON
C
C  Local array subscripts :
C
C  RVPOXY(3)  -  RVPOXY ( NFUEL )
C
C  Local variable / array dictionary:
C
C  Name    Type              Description
C  ------  ----  ----------------------------------------------------
C  LETTER   C    Used to set ASTM class
C  RVPLOS   R    RVP loss due to weathering
C  RVPOXY   R    Fuel RVP boosted by an adjustment by fuel type.
C
C  Notes:
C
C  LOCAL was changed in MOBILE4.1 to remove the ASTM-based
C  RVPUNC calculations.  Also, oxygenated and certification fuels
C  support from MOBILE4 was removed.
C  Reformulated fuel effects were added for MOBILE5.
C
      CHARACTER*1 ASTMCL,LETTER(5)
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER RFGFLG,WINFLG,RFGSYR
      LOGICAL RFGON
C
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /EVAL/ MEVAL,MONTH(12)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /OXY2/ NFUEL,MXOXMY,MAXOXY,OXYVAL(25,18,4),OXYBF(25,18,4)
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
      COMMON /RVPNAT/ RVPLIM(2),RVPDIU(2),PFUL
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION RVPOXY(3)
C
      DATA LETTER / 'A','B','C','D','E' /
C
C
C  Read in and validate the second required scenario record, which contains
C  the user supplied local data, if the user has not already provided this
C  information in one-time data.
C
      IF(LOCFLG.EQ.1) CALL RDCITY(INERR,*99)
C
C  At this point, regardless of LOCFLG, AMBT, TEMMIN & TEMMAX have been read in.
C  Check if AMBT is within the (TEMMIN,TEMMAX) range and then calculate
C  the exhaust and evaporative temperatures.
C
      IF(AMBT.LT.TEMMIN.OR.AMBT.GT.TEMMAX) CALL QUITER(AMBT,0,95,INERR)
C
      CALL CALTEM
C
C  Set RVP related parameters.
C
C  The evaluation year affects the RVP(s) to be used for this scenario:
C
      IF(ICY.LT.IUSESY) RVPICY=RVPBAS
      IF(ICY.GE.IUSESY) RVPICY=RVPIUS
C
C  Select the appropriate indices for ASTM class and phase.
C  Adjust RVP to reflect the effects of reformulated gasoline.
C  Default value for phase is 1.
C
      IPHASE=1
      IASTM =3
      RFGON=.FALSE.
      IF(RFGFLG.EQ.2) THEN
        IASTM=0
        IF(ASTMCL.EQ.LETTER(1)) THEN
          IASTM=1
        ELSEIF(ASTMCL.EQ.LETTER(2)) THEN
          IASTM=2
        ELSEIF(ASTMCL.EQ.LETTER(3)) THEN
          IASTM=3
        ELSEIF(ASTMCL.EQ.LETTER(4)) THEN
          IASTM=4
        ELSEIF(ASTMCL.EQ.LETTER(5)) THEN
          IASTM=5
        ENDIF
C
C  If user selects a reform area; yet does not specify an ASTM class
C  => error.
C
        IF(IASTM.EQ.0) THEN
          CALL QUITER(0.,0,149,INERR)
          IASTM=3
        ELSE
          ASTMCL=LETTER(IASTM)
        ENDIF
C
        DO 30 IPH=2,3
          IF(ICY.GE.RFGSYR(IPH)) IPHASE=IPH
   30   CONTINUE
C
C  Check if summer reformulated gas area or winter reformulated gas area
C   w or w/o oxygenated fuels.
C
C   WINFLG = 1 : Summer
C          = 2 : Winter
C
        IF(WINFLG.NE.1.AND.WINFLG.NE.2.OR.LOCFLG.EQ.2) THEN
          IF(MEVAL.NE.1.AND.MEVAL.NE.7) CALL QUITER(0.,MEVAL,120,INERR)
          IF(MEVAL.EQ.1) THEN
            WINFLG=2
          ELSE
            WINFLG=1
          ENDIF
        ENDIF
C
        IF(ICY.GE.RFGSYR(1)) THEN
          RFGON=.TRUE.
          IF(WINFLG.EQ.1) THEN
            RVPICY=RFGRVP(IASTM,IPHASE)
            SHRMKT(1)=0.00
            SHRMKT(2)=1.00
            SHRMKT(3)=0.00
            OXYCNT(1)=.021
            OXYCNT(2)=.000
          ELSEIF(WINFLG.EQ.2) THEN
            IF(OXYFLG.EQ.2) THEN
              IF(SHRMKT(1).GT.0.0) CALL QUITER(0.,0,147,INERR)
              IF(OXYCNT(1).EQ.0.021.AND.OXYCNT(2).EQ.0.021) GOTO 80
C
C  Oxy is not both 2.1%
C
              CHKOXY=SHRMKT(2)*OXYCNT(1)+SHRMKT(3)*OXYCNT(2)
              IF(CHKOXY.LT.0.021) CALL QUITER(0.,0,148,INERR)
C
   80 CONTINUE
            ELSEIF(OXYFLG.EQ.1) THEN
              IGASHW=2
              SHRMKT(1)=.00
              SHRMKT(2)=.70
              SHRMKT(3)=.30
              OXYCNT(1)=.015
              OXYCNT(2)=.035
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C
C  Determine the RVP boosts from oxy fuels, if necessary.
C
      IF(OXYFLG.NE.2.AND..NOT.RFGON) GOTO 20
C
      BSTRVP=0.
      BSTRVP=1.5532-0.07508*RVPICY
C
C  Determine oxy fuel commingling effect.
C
      BSTMAX=(BSTRVP/0.76)*0.0754
C
      BSTGCM=BSTMAX*(1.0-((SHRMKT(3)**2-SHRMKT(3)+0.25)/0.25))
C
      BSTCOM=0.0
      FACMCM=0.5
      IF(SHRMKT(3).GT.0.0.AND.SHRMKT(3).LT.1.0)
     *       BSTCOM=BSTGCM*(FACMCM*SHRMKT(2)+SHRMKT(1))/(1.0-SHRMKT(3))
C
C  Compute the gasoline, ether blend, and alcohol blend RVPs.
C  Include the alcohol blend RVP waiver, as necessary
C
      RVPOXY(1)=RVPICY+BSTCOM
      RVPOXY(2)=RVPOXY(1)
      RVPOXY(3)=RVPICY+BSTCOM+BSTRVP*FLOAT(IGASHW-1)
C
C  Compute a composite in-use RVP.
C
      RVPICY=0.0
      DO 10 IFUEL=1,NFUEL
       RVPICY=RVPICY+RVPOXY(IFUEL)*SHRMKT(IFUEL)
   10 CONTINUE
C
C  Save the basic, comingled, but un-weathered, RVP.
C
   20 RVPUWX=RVPICY
C
C  Apply weathering to the in-use RVP. Equation uses the calculated Hot
C  Soak temperature. Check that this weathered RVP is within limits,
C  reset and warn if necessary.
C
      RVPLOS=-2.4908+0.026196*TEMEVP(1)+0.00076898*RVPICY*TEMEVP(1)
C
      IF(RVPLOS.LT.0.0) RVPLOS=0.0
C
      RVPICY=RVPICY-RVPLOS
C
      IF(RVPICY.LT.RVPDIU(1)) CALL QUITER(RVPDIU(1),4,90,INERR)
      IF(RVPICY.LT.RVPDIU(1)) RVPICY=RVPDIU(1)
C
      IF(RVPICY.GT.RVPDIU(2)) CALL QUITER(RVPDIU(2),4,90,INERR)
      IF(RVPICY.GT.RVPDIU(2)) RVPICY=RVPDIU(2)
C
C  Calculate Uncontrolled Diurnals (UDI) for these conditions:
C  1 = FTP, 2 = In-Use, 3 = 8AM-11AM, 4 = 10AM-3PM, 5 = 8AM-2PM
C
      UDI(1)=CALUDI(RVP090,      60.,      84., 40.)
      UDI(2)=CALUDI(RVPICY,   TEMMIN,   TEMMAX,PFUL)
      UDI(3)=CALUDI(RVPICY,TEMEVP(3),TEMEVP(4),PFUL)
      UDI(4)=CALUDI(RVPICY,TEMEVP(4),   TEMMAX,PFUL)
      UDI(5)=CALUDI(RVPICY,TEMEVP(3),TEMEVP(5),PFUL)
C
C  Assign the hot soak and exhaust effect RVPs
C
      RVPHS=RVPICY
      RVPX=RVPICY
C
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
      END
