C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION FUEL(MYR,IDX2,IP,IV,ICY,BEFADJ,IY)
C
C  FUEL computes the oxy fuel, effect of industry average fuel and
C  reformulated gas adjustments.  Due to the introduction of oxidation and
C  3-way catalyst on LDTs, special handling of those technologies was needed.
C
C  Called by BEF.
C
C  Calls IEVPTR, ITAMPT, P3WPOX, and TOGFID
C
C  Input on call:
C
C    parameter list: MYR,IDX2,IP,IV,ICY,BEFADJ,IY
C    common blocks:
C     /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDFRC,FIDEFF,IVDFRC,IVDEFF
C     /EGSCAL/ TGS
C     /FLAGS1/ DTGFLG, OXYFLG
C     /FLAGS2/ NEWFLG
C     /INJECT/ PFI, TBI
C     /LOOKUP/ IVTAM,IQG
C     /MYCODE/ MY
C     /OFFSET/ OFFMTH
C     /OXY1/   SHRMKT,OXYCNT
C     /OXY2/   MXOXMY,MAXOXY,OXYVAL,OXYBF
C     /RFORM1/ IASTM,IPHASE,IAFCNT,SEAFLG,RFGON
C     /RFORM2/ RFGSYR
C
C  Output on return:
C
C    function: FUEL
C    common blocks:
C     /EVAPAR/ IFDS, ISTDT, ISTDP, FINJ
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
C  DGADJF   R    detergent gasoline adjustment to the basic emission factor
C  DPHASE   R    local detergent gasoline phase-in factor
C  FUEL1    R    1st oxy fuel adjusment factor: pre-81 LDV, HDGV, MC and
C                pre-91 LDTs
C  FUEL2    R    2nd oxy fuel adjustment factor: 1981+ LDV and 86+ LDTs
C  IAFCON   I    local variable for base fuel indicator
C  IAST     I    ASTM mapping
C  IG5      I    pointer to fuel injected sales percentage for given my
C  IGOXY    I    my group for oxy fuels
C  IPH      I    reformulated gas phase
C  IPOXY    I    pointer to the oxy value and benefit arrays
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
C
C  Notes:
C
C  FUEL was added for MOBILE4.1 as OXYADJ and renamed for MOBILE5.
C  FUEL was MODIFIED for MOBILE5v03 to include the effect of industry average
C  fuel and reformulated gasoline fuel adjustments.
C
C  22-April-1994 @ CSC-tmm Request #2-421 Added code for detergent gasoline adjustment
C     Added ICY to parameter list
C  7-May-1994 @ CSC-tmm Request #2-441 Modified code to permit the detergent gasoline
C     parameters to be read in and stored in COMMON.
C
C  29-July-1994 @ CSC-tmm Request #2-457 Replaced WINFLG with SEAFLG
C  18-August-1994 @ CSC-rwk for tmm Request #0-441, Added a modification for
C     July runs for detergent gas.  Also added a call to QUITER if the User
C     inputs calculated an unreasonable detergent gas factor.
C  24-August-1994 @CSC-pme request 2-446. The third dimension of OFFMTH array and
C                 ICAP argument were removed.
C  22-Sep-1994 @ CSC-bsg Request 426  Modularized the calculation of the
C  3-Way oxidation values, by calling a function to calculate them.  Updated
C  PHSFAC.
C  12-Dec-1994 @ CSC-ked Request 0-441 Detergent Gas Fix for July and smooth
C     error control.
C  08-Feb-1996 @ DynTel Request 2-615 Detergent gasoline effects are disabled
C                when NEWFLG is set to 5 or 6.
C  22-Jul-1996 @ DynTel-wcs Request 2-627 Alter PHSFAC values in order
C                meet VOC HC Exhaust target values for LDGV & LDGT
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in effect
C               was updated by refining phase-in period
C
      INCLUDE 'DTGGAS.I'
      INCLUDE 'EGSCAL.I'
      INCLUDE 'EVAPAR.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'INJECT.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'MYCODE.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'OXY2.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'RFORM2.I'
C
      REAL NOCAT
C
      DIMENSION PHSFAC(2,2,3),AVGFAC(3),OXYCT(2)
      DATA AVGFAC /.136, .080, .138/
C
C  Sep-22-1994 @ CSC-bsg Request 426 Updated PHSFAC.
C  Continuation code is IV
C  22-Jul-1996 @ DynTel-wcs Request 2-627 Alter PHSFAC values in order
C                meet VOC HC Exhaust target values for LDGV & LDGT
C
C                      ASTM         ASTM
C                        B            C
      DATA PHSFAC / 1.000,1.506, 1.000,1.506,
     2              1.000,1.477, 1.000,1.477,
     3              1.000,1.440, 1.000,1.440/
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
C  Begin detergent gasoline effects calculations
C    22-April-1994 @ CSC-tmm Request #2-421
C  Modified to only calculate when the new DTGFLG is equal to 1 or 2
C    7-May-1994 @ CSC-tmm Request #2-441
C
C  Look up HS/DU model year group pointers:
C
C  (1) 1 emission standard code for each fuel delivery system (IFDS) -
C          1 = carbureted (ISTDC)
C          2 = TBI (ISTDT)
C          3 = PFI (ISTDP)
C
C  (2) 1 for the percentage of vehicle type's fleet that is TBI (FINJ(1))
C      and PBI (FINJ(2)). (IG5 dimension of TBI and PFI).
C
C  8-Feb-1996 @ DynTel Request 2-615 Disables detergent gasoline effects when
C               NEWFLG is 5 or 6.
C
      IF (DTGFLG.GT.0 .AND. NEWFLG.LT.5) THEN
C
        IF(IV.GE.4.OR.MY.LT.1981) THEN
          FINJ(1) = 0.0
          FINJ(2) = 0.0
        ELSE
C
C  TBI
          IFDS=2
          ISTDT=IEVPTR(MY,IV,1)
C
          IG5=IEVPTR(MY,IV,2)
          FINJ(1)=TBI(IG5,IV)
C
C  PFI
          IFDS=3
          ISTDP=IEVPTR(MY,IV,1)
C
          IG5=IEVPTR(MY,IV,2)
          FINJ(2)=PFI(IG5,IV)
        ENDIF
C
C  Calculate the phase-in factor
C
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in year was
C              divided into phase-in start year(DGPISY) and end year(DGPIEY).
C
        IF (IY.EQ.1) THEN
          IF (ICY.LT.DGPISY) THEN
            DPHASE = 0.0
          ELSE IF (ICY.GE.DGPISY .AND. ICY.LT.DGPIEY) THEN
            DPHASE = DGPHSF
          ELSE
            DPHASE = 1.0
          ENDIF
        ELSE
C
C         July case (IY .EQ. 2  &  ICY = ICY + 1)
C
          IF (ICY.LE.DGPISY) THEN
            DPHASE = 0.0
          ELSE IF (ICY.GT.DGPISY .AND. ICY.LE.DGPIEY) THEN
            DPHASE = DGPHSF
          ELSE
            DPHASE = 1.0
          ENDIF
        ENDIF
C
C  Compute the detergent gasoline adjustment to the basic emission factor
C
        DGADJF = 1.0 - ( (FINJ(1) + FINJ(2)) * FIDFRC * (FIDEFF(IP)
     *                     * DPHASE) )
     *               - ( IVDFRC * (IVDEFF(IP) * DPHASE) )
C
C       If DGADJF is not reasonable, call QUITER and abort the run.
C
        IF (DGADJF .GT. 1.0 .OR. DGADJF .LT. 0.001) THEN
            CALL QUITER(DGADJF,IP,85,INERR)
        ENDIF
C
        FUEL = FUEL * DGADJF
C
      ENDIF
C
C  End of Detergent gasoline effects calculations
C  Sep-22-1994 @ CSC-bsg Request 426  Modularized the calculation of the
C  3-Way oxidation values, by calling a function to calculate them.
C  Compute the number of catalyst and non-catalyst vehicles,
C
      CALL P3WPOX(IV,P3W,POX)
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
   80 IF(IP.EQ.1.AND.RFGON.AND.SEAFLG.EQ.1) THEN
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
        OFMTH=OFFMTH(IDX2,IV)
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
