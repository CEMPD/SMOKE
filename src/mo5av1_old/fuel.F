C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION FUEL(MYR,IDX2,IP,IV,BEFADJ,ICAP)
C
C  FUEL computes the oxy fuel, effect of industry average fuel and
C  reformulated gas adjustments.  Due to the introduction of oxidation and
C  3-way catalyst on LDTs, special handling of those technologies was needed.
C
C  Called by BEF.
C
C  Calls ITAMPT
C
C  Input on call:
C
C    parameter list: MYR,IDX2,IP,IV,BEFADJ,ICAP
C    common blocks:
C     /EGSCAL/ TGS
C     /FLAGS1/ OXYFLG
C     /LOOKUP/ IVTAM,IQG
C     /MYCODE/ MY
C     /OFFSET/ OFFMTH(25,8,2)
C     /OXY1/   SHRMKT,OXYCNT
C     /OXY2/   MXOXMY,MAXOXY,OXYVAL,OXYBF
C     /RFORM1/ IASTM,IPHASE,IAFCNT,WINFLG,RFGON
C     /RFORM2/ RFGSYR
C
C  Output on return:
C
C    function: FUEL
C    common blocks:
C     /LOOKUP/ IVTAM,IQG
C     /MYCODE/ MY
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  AVGFAC   R    percent decrease in emissions from industry average fuel to
C                Indolene fuel
C  BENFAC   R    benefit factor for determining oxy adjustment
C  CAT      R    fraction of vehicles that are catalyst equipped
C  FUEL1    R    1st oxy fuel adjusment factor: pre-81 LDV, HDGV, MC and
C                pre-91 LDTs
C  FUEL2    R    2nd oxy fuel adjusment factor: 1981+ LDV and 86+ LDTs
C  IAST     I    ASTM mapping
C  IAFCON   I    local variable for base fuel indicator
C  IGOXY    I    my group for oxy fuels
C  IPOXY    I    pointer to the oxy value and benefit arrays
C  IPH      I    reformulated gas phase
C  IVI      I    reformulated gas phase vehicle mapping
C  MYMAP    I    model year for mapping into the oxy value arrays
C  NOCAT    R    fraction of vehicles that are not catalyst equipped
C  OFMTH    R    value holder for the methane offset to basic emission rate
C  OXYADJ   R    value holder for the oxygenated fuel benefit adjustment
C  OXYCT    R    local variable for the oxygen content
C  PHSFAC   R    HC reduction for reformulated gas
C  P3W      R    percent of vehicles with 3-way catalyst
C  POX      R    percent of vehicles with ox catalyst
C  TOG1     R    value holder for the TOG FID correction factor
C  TOG2     R    value holder for the TOG FID oxygen correction factor
C  WT1      R    1st fuel adjustment factor weighting factor
C  WT2      R    2nd fuel adjustment factor weighting factor
C  X1       R    lower x-value (IPOXY) for interpolation
C  X2       R    upper x-value (IPOXY) for interpolation
C  Y1       R    lower y-value (benefit factor) for interpolation
C  Y2       R    upper y-value (benefit factor) for interpolation
C
C  Notes:
C
C  FUEL was added for MOBILE4.1 as OXYADJ and renamed for MOBILE5.
C  FUEL was MODIFIED for MOBILE5v03 to include the effect of industry average
C  fuel and reformulated gasoline fuel adjustments.
C
      CHARACTER*1 ASTMCL
      INTEGER OXYFLG
      INTEGER WINFLG,RFGSYR
      LOGICAL RFGON
C
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /OXY2/ NFUEL,MXOXMY,MAXOXY,OXYVAL(25,18,4),OXYBF(25,18,4)
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
C
      REAL NOCAT
C
      DIMENSION PHSFAC(2,2,3),AVGFAC(3),OXYCT(2)
      DATA AVGFAC /.136, .080, .138/
C
C  Continuation code is IV
C                      ASTM         ASTM
C                        B            C
      DATA PHSFAC / 1.122,1.801, 1.117,1.779,
     2              1.076,1.771, 1.067,1.738,
     3              1.051,1.720, 1.044,1.693/
C
C  Initialize
C
      FUEL=1.0
C
C  If vehicle is diesel fueled, no adjustment is necessary.
C
      IF(IV.GE.5.AND.IV.LE.7) RETURN
C
C  Reset.
C
      IAFCON=IAFCNT
      MY=MYR
C
C  Local variable for oxygen content cannot be greater than 2.7% for HC.
C
      DO 10 IFUEL=2,3
        OXYCT(IFUEL-1)=OXYCNT(IFUEL-1)
        IF(IP.EQ.1.AND.OXYCNT(IFUEL-1).GT..027) OXYCT(IFUEL-1)=.027
   10 CONTINUE
C
C  Compute the number of catalyst and non-catalyst vehicles,
C
      P3W=0.0
      POX=0.0
C
      IF(IV.LE.4) THEN
        IVTAM=IV
C
        DO 15 IQG=5,6
        IG1=ITAMPT(1)
        P3W=P3W+TGS(IG1,IQG,IVTAM)
   15   CONTINUE
C
        DO 20 IQG=2,3
        IG1=ITAMPT(1)
        POX=POX+TGS(IG1,IQG,IVTAM)
   20   CONTINUE
        POX=POX-P3W
      ENDIF
C
      CAT=P3W+POX
      NOCAT=1.0-CAT
C
C  Compute the industry average fuel adjustment factor.
C  For NOx, only use the percentage of 3-way catalysts.
C
      IF(RFGON.AND.IPHASE.GE.2) IAFCON=339
C
      IF(IP.EQ.3) THEN
        FUEL=FUEL+(FUEL*P3W*(IAFCON-50)/(339.0-50.0)*
     *    (AVGFAC(IP)/(1.0-AVGFAC(IP))))
      ELSE
        FUEL=FUEL+(FUEL*CAT*(IAFCON-50)/(339.0-50.0)*
     *    (AVGFAC(IP)/(1.0-AVGFAC(IP))))
      ENDIF
C
C  Return if no reformulated or oxygenated fuel effects are to be calculated.
C
      IF(.NOT.RFGON.AND.OXYFLG.EQ.1) GOTO 90
C
C  Skip RFG and Oxy effects if NOx
C
      IF(IP.EQ.3) GOTO 90
C
C  Initialize oxy/reform variables.
C
      FUEL1=0.0
      FUEL2=0.0
      WT1=1.0
      WT2=1.0
C
C  Compute the benefit factor and the 1st oxy fuel adjustment factor,
C  as necessary.
C
      IF((IV.EQ.1.AND.MY.LT.1981).OR.
     *  ((IV.EQ.2.OR.IV.EQ.3).AND.MY.LE.1990).OR.IV.GE.4) THEN
        BENFAC=NOCAT*OXYFAC(IP,1)+CAT*OXYFAC(IP,2)
        FUEL1=SHRMKT(1)
     *    +SHRMKT(2)*(1.0-OXYCT(1)*BENFAC)
     *    +SHRMKT(3)*(1.0-OXYCT(2)*BENFAC)
      ENDIF
C
C  Compute the 2nd oxy factor, as necessary.
C
      IF((IV.EQ.1.AND.MY.LT.1981).OR.
     * ((IV.EQ.2.OR.IV.EQ.3).AND.MY.LT.1986).OR.
     *   IV.GE.4) GOTO 80
C
C  Adjust the basic emisson factor value for base fuel.
C
      BEFADJ=BEFADJ*FUEL
C
C  For 1986-87 LDTs, map the model year.
C
      MYMAP=MY
      IF((IV.EQ.2.OR.IV.EQ.3).AND.
     *   (MY.EQ.1986.OR.MY.EQ.1987)) MYMAP=MYMAP+1
C
C  Look up the bracketing values and set up for interpolation.
C
      IGOXY=MYMAP-1981+1
      IF(IGOXY.GT.MXOXMY) IGOXY=MXOXMY
C
C  Loop if pollutant is HC to take care of separate ether and ethanol curves.
C
      IF(IP.GT.1) THEN
        IPOXY=IP+1
      ELSE
        IPOXY=1
      ENDIF
C
   40 DO 50 IG=1,MAXOXY
      IF(BEFADJ.LE.OXYVAL(IG,IGOXY,IPOXY)) GOTO 60
   50 CONTINUE
C
      BENFAC=OXYBF(MAXOXY,IGOXY,IPOXY)
      GOTO 70
C
   60 IF(IG.EQ.1) THEN
        BENFAC=OXYBF(IG,IGOXY,IPOXY)
        GOTO 70
      ENDIF
C
      IF(IG.GT.1) THEN
      X1=OXYVAL(IG-1,IGOXY,IPOXY)
      Y1=OXYBF(IG-1,IGOXY,IPOXY)
      X2=OXYVAL(IG,IGOXY,IPOXY)
      Y2=OXYBF(IG,IGOXY,IPOXY)
      ENDIF
C
C  Interpolate the benefit factor.
C
      BENFAC=(BEFADJ*(Y2-Y1)+X2*Y1-X1*Y2)/(X2-X1)
C
C  Compute the 2nd oxy fuel adjustment factor.
C
   70 IF(IPOXY.EQ.1) THEN
        FUEL2=
     *     SHRMKT(2)*(1.0-OXYCT(1)*BENFAC)
        IPOXY=2
C
C  HC Loop.
C
        GOTO 40
      ENDIF
      IF(IPOXY.EQ.2) THEN
        FUEL2=SHRMKT(1)
     *    +FUEL2
     *    +SHRMKT(3)*(1.0-OXYCT(2)*BENFAC)
      ELSE
        FUEL2=SHRMKT(1)
     *    +SHRMKT(2)*(1.0-OXYCT(1)*BENFAC)
     *    +SHRMKT(3)*(1.0-OXYCT(2)*BENFAC)
      ENDIF
C
C  Determine the 1986-90 LDT oxy weighting factors.
C
      IF((IV.EQ.2.OR.IV.EQ.3).AND.(MY.GE.1986.AND.MY.LE.1990)) THEN
      WT1=1.0-P3W
      WT2=P3W
      ENDIF
C
C  Apply the summer reformulated fuel phase factor to the fuel factors, if
C  they were set since initialization.  Map to the appropriate phase factor.
C  Compute the weighted oxy fuel adjustment factor.
C
   80 IF(IP.EQ.1.AND.RFGON.AND.WINFLG.EQ.1) THEN
        IPH=1
        IF(IPHASE.EQ.3) IPH=2
        IAST=1
        IF(IASTM.GE.3) IAST=2
        IVI=IV
        IF(IV.EQ.4) IVI=3
        IF(IV.EQ.8) IVI=1
        IF(FUEL1.NE.0.0) FUEL1=1.0-((1.0-FUEL1)*PHSFAC(IPH,IAST,IVI))
        IF(FUEL2.NE.0.0) FUEL2=1.0-((1.0-FUEL2)*PHSFAC(IPH,IAST,IVI))
      ENDIF
C
C  Compute the weighted oxy fuel adjustment factor.
C  Convert HC from TOG to THC
C
      OXYADJ=WT1*FUEL1+WT2*FUEL2
      IF(IP.EQ.1) THEN
        OFMTH=OFFMTH(IDX2,IV,ICAP)
        TOG1=TOGFID(MY,IV,1)
        TOG2=TOGFID(MY,IV,2)
        FUEL=FUEL*
     *    (((TOG1*(BEFADJ-OFMTH)+OFMTH)*OXYADJ-OFMTH)/TOG2+OFMTH)/BEFADJ
      ELSE
        FUEL=FUEL*OXYADJ
      ENDIF
C
   90 RETURN
      END
