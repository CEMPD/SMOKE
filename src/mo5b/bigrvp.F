C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGRVP(ICY,IP,IG,IVTEMP)
C
C  BIGRVP computes the combined temperature and RVP correction factor algorithm
C  used when the pollutant's temperature >= 75 degrees (FTP standard) and
C  the RVP cf alone when T < 75 degrees.
C
C  See RVPEXH for an overall description of the RVP cf cases.
C
C  Called by BIGTCF.
C
C  Calls AMIN1, AMAX1, and P3WPOX.
C
C  Input on call:
C
C    parameter list: ICY,IP,IG,IVTAM
C    common blocks:
C    /CITRV2/ RVPX,RVP090
C    /RFORM1/ IPHASE,SEAFLG,RFGON
C    /RVPEX1/ RG23,TTL
C    /TEMPC3/ MDLOHI
C    /TEMPS/  TEMEXH
C
C  Output on return:
C
C    common blocks:
C    /TEMPC3/ TCF
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CFB      R    RVP cf at begin point of temperature interpolation interval
C  CFL      R    RVP cf at end   point of temperature interpolation interval
C  P3W      R    Percent 3-way catalyst
C  RDIFA    R    the larger of 0.0 and (applicable RVP - 9.0)
C  REAT     R    adjusted RVP exhaust IP temperature
C  REUT     R    RVP exhaust IP unadjusted temperature
C  RFGNOX   R    6.8% reduction to RFG NOx percent three way emissions, for
C                Summer, RFG applies, and in Phase II.
C  RIVAL    R    RVP used in TVF calculations
C  RVP1     R    New RVP
C  TDIFA    R    the larger of 0.0 & (adjusted T capped at 95.0) - 75.0)
C  TDIFU    R    (unadjusted exhaust T capped at 110.0) - 75.0
C
C
C  Notes:
C
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1.
C
C  29 July 1994 @ CSC-tmm Request #2-457 Replaced WINFLG with SEAFLG
C  of the RFG NOx effects for Phase II.
C  Sep-26-1994 @ CSC-bsg Request 426  Added IVTAM parameter for calling 
C  P3WPOX.
C  Oct-11-1994 @ CSC-bsg Request 426  Clarified documentation.  Also replaced
C  greater than or equal to 2000 with IPHASE equal 3 to avoid July
C  interpolation problems.
C  Oct-13-1994 @ CSC-bsg Request 426  Clarified documentation.  Also relocated
C  and reworked the RVP reassignment based on INOX loop instead of RFGNOX
C  and only occurring in Summer.
C  Oct-14-1994 @ CSC-bsg Request 426 Further restriction of the definition of
C  RFGNOX=2 was added.  Now RFGNOX=2 when it is Summer (SEAFLG=1).
C  Oct-18-1994 @ CSC-bsg Request 426  Changed the definitions of what RVP TCF
C  will depend on.  The 6.8% reductions will apply only to 8.7 RVP TCF
C  Oct-19-1994 @ CSC-bsg Request 426  Confirmed logic for weathering and
C  applied weathering in appropriate cases.
C  Oct-28-1994 @ CSC-bsg Request 426  Changed definitions of when to use
C  InUse RVP for NOX.  Also whenever using InUse RVP the RVP will need to be
C  weathered.
C  Dec-09-1994 @ CSC-ked Request 475  (RFG NOX Clean-up)  Modified RFG NOX
C  further to include Phase I and Winter.
C
      INCLUDE 'CITRV1.I'
      INCLUDE 'CITRV2.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'MYCODE.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'RFORM2.I'
      INCLUDE 'RVPEX1.I'
      INCLUDE 'TEMPC3.I'
      INCLUDE 'TEMPS.I'
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
C  Summer Reform will use 8.7 RVP for HC and NOx.
C
      IF(IP.NE.2) THEN
        IF(RFGON.AND.SEAFLG.EQ.1) THEN
          RIVAL=8.7
C
C  Apply weathering to the NOx and HC RVP.  Equation uses the calculated Hot
C  Soak temperature.  Check that this weathered RVP is within limits.
C
          RVPLSS=-2.4908
     &           +0.026196*TEMEVP(1)
     &           +0.00076898*RIVAL*TEMEVP(1)
          IF(RVPLSS.GE.0.0) RIVAL=RIVAL-RVPLSS
        ENDIF
C
C  Apply weathering to the NOx in an Oxy Program (and not RFGON).  First the
C  RVP needs to be reset to the user input and then weathered.
C
        IF(IP.EQ.3) THEN
          IF((OXYFLG.EQ.2.AND..NOT.RFGON).OR.
     &       (RFGON.AND.SEAFLG.EQ.2)) THEN
            IF(ICY.LT.IUSESY) RIVAL=RVPBAS
            IF(ICY.GE.IUSESY) RIVAL=RVPIUS
            RVPLSS=-2.4908
     &            +0.026196*TEMEVP(1)
     &            +0.00076898*RIVAL*TEMEVP(1)
            IF(RVPLSS.GE.0.0) RIVAL=RIVAL-RVPLSS
          ENDIF
        ENDIF
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
      GOTO 40
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
C  Apply an RFG NOX effect of a 6.8% reduction to the percent three way
C  only.  This is accomplished by weighting together the percent three way and
C  the complement (1-P3W) with no reduction.  TCF is adjusted only -- when
C  RFG applies and calculating NOx emissions.
C
   40 IF(RFGON .AND. (IP .EQ. 3)) THEN
        IF ((SEAFLG .EQ. 1) .AND. (IPHASE .EQ. 3)) THEN
            RFGNOX = 1.0 - 0.068
        ELSE
            RFGNOX = 1.0 - 0.015
        ENDIF
        CALL P3WPOX(IVTEMP,P3W,POX)
        DO 50 IB = 1, 3
           TCF(IB) = ( P3W * RFGNOX + ( 1.0 - P3W ) ) * TCF(IB)
   50   CONTINUE
      ENDIF
C
      RETURN
      END
