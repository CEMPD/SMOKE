      SUBROUTINE FUELPAR(ICY,INERR) 
C 
C  FUELPAR initializes fuel specific variables. 
C 
C  Called by PRCSCN. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  17 Jul 01 @EPA-bag Bug , Changed the subroutine name to avoid conflict with 
C            array name in BDSPEED9.I 
C  07 Mar 01 @EPA-djb Bug 318, Added conditional to setting of S_PPM for RFG. 
C  23 Aug 00 @EPA-djb Bug 220, Initialized RFGON to .FALSE. and 
C            allowed RFG oxy to affect RVP. 
C  21 Jul 00 @EPA-bag  running loss; removed CITCIN.I 
C  26 Jun 00 @EPA-elg  Set JPHASE conditional statement. 
C  12 Jun 00 @EPA-djb Bug99 Added RFG RVP calculation. 
C  19 May 00 AIR Task 03: Added WXRVP_FLAG (Weathered RVP). 
C  17 May 00 AIR Task 03: Removed the TEMPHRLY dependency for weathered RVP. 
C  31 May 00 @EPA-ddj Bug102, Removed call to QUITER for warning 
C            message 170 (detergent additive regulations). 
C  19 May 00 @EPA        revised S_FUEL model year index, IYR 
C  28 Mar 00 @DynTel-rwk 1-028 reinserted RFG flag. 
C  25 Feb 00 @DynTel-ddj 1-028, Removed obsolete ASTM code. 
C  04 Apr 00: AIR Task 02: Removed ASTM class code, replaced with a 
C             hardcoded RFG Region 2 (ASTM of C). This code must be fixed 
C             when the RFG Task is installed. Also removed warning for  
C             gasoline detergent additive regulations (1994) 
C  03 Dec 98 @DynTel-BAG 2-664 Replaced AMBT with TEMHRLY(IH) and deleted 
C            calls to CALUDI. Added hourly array for hot soak hourly RVP's. 
C  24 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input after call to RDCITY: 
C 
C    argument list: ICY 
C    common blocks: 
C    /CITRV1/ RVPBAS 
C    /CITRV2/ RVPHS 
C    /EVAL/   MEVAL 
C    /FLAGS1/ OXYFLG 
C    /FLAGS5/ WXRVP_FLAG 
C    /OXY1/   SHRMKT,IGASHW,OXYCNT 
C    /OXY2/   NFUEL 
C    /RFORM1/ JPHASE,SEAFLG,RFGFLG 
C    /RFORM2/ RFGSYR,RFGRVP 
C    /RVPNAT/ PFUL 
C    /TEMPS/  TEMHRLY,TEMMIN,TEMMAX 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /CITRV1/ RVPICY,RVPUWX 
C    /CITRV2/ RVPHS,RVPNOX,OXYRVP 
C    /RFORM1/ IPHASE,RFGON,SEAFLG 
C 
C  Local array subscripts : 
C 
C  RVPOXY(3)  -  RVPOXY ( NFUEL ) 
C 
C  Local variable / array dictionary: 
C 
C  Name        Type              Description 
C  ------      ----  ---------------------------------------------------- 
C  IH           I    Hourly loop index 
C  NOXLOS       R    RVP loss due to weathering for NOX calculations 
C                    for RFG and OXY cases. 
C  RFG_DEFAULT  I    Flag indicting whether default RFG fuel parameters are 
C                    to be retained for RFG calculation. 
C  RVPLOS       R    RVP loss due to weathering 
C  RVPTMP       R    temporary RVP value used to determine hourly weathered 
C                    hot soak RVP 
C  RVPOXY       R    Fuel RVP boosted by an adjustment by fuel type. 
C 
C  Notes: 
C 
C  LOCAL was changed in MOBILE4.1 to remove the ASTM-based 
C  RVPUNC calculations.  Also, oxygenated and certification fuels 
C  support from MOBILE4 was removed. 
C  Reformulated fuel effects were added for MOBILE5. 
C 
C  TMM@ARC, Request 2-366, 1/4/94 
C    Display warning about Winter Reformulated Gasoline rules 
C 
C  29 July 1994 @ CSC-tmm Request #2-457 
C     Added code to preserve the value of WINFLG as it was supplied 
C     by the user. A new flag, SEAFLG, was added to contain the 
C     calculated value of the season. If WINFLG is zero, then  
C     SEAFLG is determined from MEVAL. If WINFLG is 1 or 2, then SEAFLG 
C     is set to the value of WINFLG. This will permit the use of WINFLG 
C     when LOCFLG is equal to 2. Previously, WINFLG was ignored if 
C     LOCFLG was equal to 2. 
C   Oct-21-97 @DynTel-gjr 2-645 Removed RVPIUS and IUSESY. 
C   4 Jun 98 @DynTel-MLA 2-659  All input functions, including a call to 
C            the obsolete (in Mobile6) fuction RDCITY have been removed 
C            from this subroutine. This subroutine is called only after 
C            all of the scenario-specific data has been read. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'OXY2.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'RVPNAT.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'TEMPS.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: INERR 
C 
      INTEGER IFUEL 
      INTEGER IYR 
      INTEGER RFG_DEFAULT 
      INTEGER, PARAMETER :: INO = 0 
C 
      REAL CHKOXY 
      REAL BSTMAX 
      REAL BSTRVP 
      REAL BSTGCM 
      REAL BSTCOM 
      REAL FACMCM 
      REAL RVPOXY 
      REAL RVPTMP 
      REAL RVPLOS 
      REAL TRISE 
C 
      DIMENSION RVPOXY(3) 
C 
C  Set RVP related parameters. 
C 
      RVPICY=RVPBAS 
C 
C     Begin RFG block 
C 
      RFGON=.FALSE. 
      IF(RFGFLG.EQ.2) THEN  ! RFG has been selected by the user 
C 
C  Select the appropriate indices for phase. 
C  Adjust RVP to reflect the effects of reformulated gasoline. 
C  Default value for phase is 1. 
C 
      IF(ICY.LT.RFGSYR(2)) THEN 
          JPHASE = 1 
      ELSE 
         JPHASE = ICY - 1998 
         IF(JPHASE.GE.18) JPHASE = 18 
      END IF 
C 
      IF(ICY.GE.2015) THEN 
         IYR = 16 
      ELSEIF(ICY.GT.2000) THEN 
         IYR = ICY - 1999 
      ELSE 
         IYR = 1 
      END IF 
C  Check if summer reformulated gas area or winter reformulated gas area 
C   w or w/o oxygenated fuels. 
C 
C   WINFLG = 1 : Summer 
C          = 2 : Winter 
C 
      IF (WINFLG.EQ.0) THEN 
        IF(MEVAL.NE.1.AND.MEVAL.NE.7) CALL QUITER(0.,MEVAL,120,INERR) 
        IF(MEVAL.EQ.1) THEN 
           SEAFLG=2 
        ELSE 
           SEAFLG=1 
        ENDIF 
      ELSE 
         SEAFLG = WINFLG 
      ENDIF 
C 
      IF(ICY.GE.RFGSYR(1)) THEN 
         RFGON=.TRUE. 
         IF(SEAFLG.EQ.1) THEN        ! If summer 
            RVPICY = RFGRVP(RFG_REGION,JPHASE) 
            IF(ICY.GE.2000) S_PPM = S_FUEL(1,2,IYR,2) 
            SHRMKT(1)=0.00 
            SHRMKT(2)=1.00 
            SHRMKT(3)=0.00 
            OXYCNT(1)=.021 
            OXYCNT(2)=.000 
         ELSEIF(SEAFLG.EQ.2) THEN    ! Else if winter 
C 
C   TMM@ARC, Request 2-366, 1/4/94 
C     Display warning about winter reformulated gasoline rules 
C 
            CALL QUITER(0.,0,112,INERR) 
C 
            RFG_DEFAULT = 0 
C 
            IF(OXYFLG.EQ.2) THEN      ! OXY values adjusted 
              IF(SHRMKT(1).GT.0.0) THEN 
                 RFG_DEFAULT = 1 
                 CALL QUITER(0.,0,147,INERR) 
              END IF 
              CHKOXY=SHRMKT(2)*OXYCNT(1)+SHRMKT(3)*OXYCNT(2) 
              IF(CHKOXY.LT.0.021) THEN 
                 RFG_DEFAULT = 1 
                 CALL QUITER(0.,0,148,INERR) 
              END IF 
            ELSEIF(OXYFLG.EQ.1) THEN  ! Retain default OXY values 
                  RFG_DEFAULT = 1 
            END IF 
C 
            IF(RFG_DEFAULT .EQ. 1) THEN 
              IGASHW=2 
              SHRMKT(1)=.00 
              SHRMKT(2)=.70 
              SHRMKT(3)=.30 
              OXYCNT(1)=.015 
              OXYCNT(2)=.035 
            ENDIF 
         ENDIF   !Season block 
      ELSE 
        RFGON=.FALSE. 
      ENDIF      !calendar year check 
      ENDIF      !RFG block 
C 
C  Determine the RVP boosts from oxy fuels, if necessary. 
C 
      IF(OXYFLG.EQ.2.OR.RFGON) THEN 
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
        DO IFUEL=1,NFUEL 
          RVPICY=RVPICY+RVPOXY(IFUEL)*SHRMKT(IFUEL) 
        END DO 
      ENDIF 
C 
C  Save the basic, comingled, but un-weathered, RVP.(Used in REFUEL.) 
C 
      RVPUWX=RVPICY 
C 
C  Compute the hot soak temperature from the daily maximum 
C  and minimum temperatures. 
C 
      TEMEVP=(TEMMAX+TEMMIN)/2. 
      TRISE=TEMMAX-TEMMIN 
C 
      IF(TEMMAX.GE.40..AND.TRISE.GT.0.) 
     *   TEMEVP=-1.7474+1.0290*TEMMIN+0.99202*TRISE 
     *          -0.0025173*TEMMIN*TRISE 
C 
C  Apply weathering to the in-use RVP. Equation uses the calculated 
C  daily temperature based on the Hot Soak temperature method. Check that 
C  this weathered RVP is within limits, reset and warn if necessary. 
C 
      IF(WXRVP_FLAG==INO) THEN 
C 
C  No fuel weathering case. 
C 
        RVPNOX=8.7 
        OXYRVP=RVPBAS 
      ELSE 
C 
C  Base weathered RVP calculation. 
C 
        RVPLOS = -2.4908+0.026196*TEMEVP+ 
     *            0.00076898*RVPICY*TEMEVP 
        IF(RVPLOS.LT.0.0) RVPLOS=0.0 
        RVPICY=RVPICY-RVPLOS 
C 
C  A weathered 8.7 psi RVP value is needed for reformulated fuel effects 
C  for HC and NOx emissions. 
C 
        RVPLOS = -2.4908+0.026196*TEMEVP+ 
     *            0.00076898*8.7*TEMEVP 
        RVPNOX=8.7-RVPLOS 
C 
C  A weathered user supplied RVP value without commingling effects is 
C  needed for oxygenated fuel effects for NOx emissions. 
C 
        RVPLOS = -2.4908+0.026196*TEMEVP+ 
     *            0.00076898*RVPBAS*TEMEVP 
        OXYRVP=RVPBAS-RVPLOS 
      ENDIF 
C 
C  Determine Hot Soak RVP value.  Start with base weathered RVP value. 
C  RVPHS is limited by RVPDIU. 
C 
      RVPTMP=RVPICY 
C 
      IF(RVPTMP.LT.RVPDIU(1)) CALL QUITER(RVPDIU(1),4,90,INERR) 
      IF(RVPTMP.LT.RVPDIU(1)) RVPTMP=RVPDIU(1) 
C 
      IF(RVPTMP.GT.RVPDIU(2)) CALL QUITER(RVPDIU(2),4,90,INERR) 
      IF(RVPTMP.GT.RVPDIU(2)) RVPTMP=RVPDIU(2) 
C 
C     Assign the hot soak effect RVPs 
C 
      RVPHS=RVPTMP 
C 
      RETURN 
      END 
