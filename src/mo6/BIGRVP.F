      SUBROUTINE BIGRVP(IP,IG) 
C 
C  BIGRVP computes the combined temperature and RVP correction factor 
C  algorithm used when the pollutant's temperature >= 75 degrees (FTP 
C  standard) and the RVP cf alone when T < 75 degrees. 
C 
C  See RVPEXH for an overall description of the RVP cf cases. 
C 
C  Called by BIGTCF. 
C 
C  Calls AMIN1, AMAX1, and P3WPOX. 
C 
C  Changes : 
C 
C  03 Jul 01 @EPA-djb Bug329, removed RFGNOX adjustment (this is now 
C            done in FUELCF).  Added text to explain RFG/Oxy code. 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; 
C            removed MYCODE.I and added model year, MYBGRVP, to parameter 
C            list 
C  26 Jun 00 @EPA-elg  Set JPHASE conditional statement. 
C  12 Jun 00 @EPA-djb Bug99 Move weathered RVP calculations to subroutine LOCAL. 
C  01 Jun 00 AIR Task 03: Restored code to MOBILE5b then updated to M6. 
C  16 May 00 @EPA-ddj BUG27, Changed RVP weathering calculation to be 
C            consistent with MOBILE5b, except for hourly temperature 
C            rather than weighted temperature.  Added to bug list 
C            as BUG99. 
C  14 Apr 99 @DynTel-ddj 2-000 Changed all temperatures to hourly 
C            temperature.  Added IH to parameters being passed to 
C            BIGRVP. 
C  16 Mar 99 @DynTel-bg 2-664  Removed evaporative temperature array, TEMEVP, 
C            and temporarily replaced it with a daily average temperature. 
C  14 Aug 98 @DynTel-zk 2-000  Explicit typing 
C  18 November 1997 Dyntel@MMS 2-645 removed IP index of MDLOHI 
C  12 November 1997 Dyntel@MMS 2-645 removed index of TEMEXH(IP) 
C  Oct-21-97 DynTel@gjr 2-645 Removed RVPIUS and IUSESY. ICY no longer 
C            needed in this module. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified 
C            the code to use the symbolic vehicle type names and vehicle 
C            type vectors that were introduced in Mobile6. 
C  Dec-09-1994 @ CSC-ked Request 475  (RFG NOX Clean-up)  Modified RFG 
C            NOX further to include Phase I and Winter. 
C  Oct-28-1994 @ CSC-bsg Request 426  Changed definitions of when to use 
C            InUse RVP for NOX.  Also whenever using InUse RVP the RVP 
C            will need to be weathered. 
C  Oct-19-1994 @ CSC-bsg Request 426  Confirmed logic for weathering and 
C            applied weathering in appropriate cases. 
C  Oct-18-1994 @ CSC-bsg Request 426  Changed the definitions of what 
C            RVP TCF will depend on.  The 6.8% reductions will 
C            apply only to 8.7 RVP TCF 
C  Oct-14-1994 @ CSC-bsg Request 426 Further restriction of the 
C            definition of RFGNOX=2 was added.  Now RFGNOX=2 when it is 
C            Summer (SEAFLG=1). 
C  Oct-13-1994 @ CSC-bsg Request 426  Clarified documentation.  Also 
C            relocated and reworked the RVP reassignment based on INOX 
C            loop instead of RFGNOX and only occurring in Summer. 
C  Oct-11-1994 @ CSC-bsg Request 426  Clarified documentation.  Also 
C            replaced greater than or equal to 2000 with IPHASE equal 3 
C            to avoid July interpolation problems. 
C  Sep-26-1994 @ CSC-bsg Request 426  Added IVTAM parameter for calling 
C            P3WPOX. 
C  29 July 1994 @ CSC-tmm Request #2-457 Replaced WINFLG with SEAFLG 
C            of the RFG NOx effects for Phase II. 
C       
C  Input on call: 
C 
C    argument list: IP,IG,IVTAM 
C    common blocks: 
C    /CITRV2/ RVPICY,RVP090 
C    /RFORM1/ JPHASE,SEAFLG,RFGON 
C    /RVPEX1/ RG23,TTL 
C    /TEMPC3/ MDLOHI 
C    /TEMPS/  TEMHRLY 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
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
C  CFB      R    RVP cf at begin point of temperature interpolation 
C                interval 
C  CFL      R    RVP cf at end   point of temperature interpolation 
C                interval 
C  P3W      R    Percent 3-way catalyst 
C  RDIFA    R    the larger of 0.0 and (applicable RVP - 9.0) 
C  REAT     R    adjusted RVP exhaust IP temperature 
C  REUT     R    RVP exhaust IP unadjusted temperature 
C  RIVAL    R    RVP used in TVF calculations 
C  RVP1     R    New RVP 
C  TDIFA    R    the larger of 0.0 & (adjusted T capped at 95.0) - 75.0) 
C  TDIFU    R    (unadjusted exhaust T capped at 110.0) - 75.0 
C 
C 
C  Notes: 
C 
C  Oxygenated and cert. fuels support from MOBILE4 was removed in 
C  MOBILE4.1. 
C 
C 
      IMPLICIT NONE  
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'RVPEX1.I' 
      INCLUDE 'TEMPC3.I' 
      INCLUDE 'TEMPS.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IG 
C 
C     Declare local variables 
C 
      INTEGER IB 
      INTEGER IGR 
      REAL CFB 
      REAL CFL 
      REAL RDIFA 
      REAL REAT 
      REAL RIVAL 
      REAL RVP1 
      REAL REUT 
      REAL TDIFA 
      REAL TDIFU 
C 
C  Reduce supplied RVP to 11.7 for the RVP cf calculation only. 
C 
      RVP1 = AMIN1(RVPICY,11.7) 
C 
C  If RVP x T temp is .GT. 95 degrees, then reduce it to 95 for here 
C  only.  If T term temp is .GT. 110 degrees, then reduce it to 110 for 
C  here only. 
C 
      REAT = AMIN1(TEMEXH,95.0) 
      TDIFA = REAT - 75.0 
      TDIFA = AMAX1(TDIFA,0.0) 
      REUT = AMIN1(TEMEXH,110.0) 
      TDIFU = REUT - 75.0 
C 
C  Identify RVP to be used. RVP can drop to 7.0 if >75 degrees, else 9.0 
C 
      RIVAL = RVP1 
C 
C  Special cases for HC and NOx with Oxy or RFG. 
C 
      IF(IP.NE.2) THEN 
C 
C  Use the weathered 8.7 RVP for NOx and HC for cases with RFG 
C  in the summer season. 
C  This assures that summer RFG has no effect on NOx emissions compared 
C  to a base fuel with an 8.7 RVP. 
C 
        IF(RFGON.AND.SEAFLG.EQ.1) RIVAL = RVPNOX 
C 
C  Use the weathered user supplied RVP for NOx with an Oxy Program 
C  and without RFG.  This also applies to cases with RFG in winter. 
C  This assures that oxy fuels have no effect on NOx emissions. 
C 
        IF(IP.EQ.3) THEN 
          IF((OXYFLG.EQ.2.AND..NOT.RFGON).OR. 
     &       (RFGON.AND.SEAFLG.EQ.2)) THEN 
            RIVAL=OXYRVP 
          ENDIF 
        ENDIF 
      ENDIF 
C 
C  Determine RVP difference from 9.0, set limits. 
C 
      RDIFA=RIVAL-RVP090 
      IF(REAT.GT.75.0) THEN 
         RDIFA = AMAX1(RDIFA,-2.0) 
      ELSE 
         RDIFA = AMAX1(RDIFA,0.0) 
      END IF 
C 
C  Compute the cf for each bag. 
C 
      IGR=IG-1 
C 
C  Case 1: one equation produces a combined RVP & high/neutral T cf. 
C 
      IF(MDLOHI.EQ.1) GOTO 20 
C 
C  Evaluate the combined RVP & high T cf formula. 
C  Except for NOx, combined RVP & T cf is never less than 1.0. 
C 
      DO 10 IB = 1, 3 
      TCF(IB) = RG23(1,IB,IP,IGR)*RDIFA + 
     &          RG23(2,IB,IP,IGR)*TDIFU + 
     &          RG23(3,IB,IP,IGR)*RDIFA*TDIFA 
      IF(IP.NE.3) TCF(IB) = AMAX1(TCF(IB),0.0) 
      TCF(IB) = EXP(TCF(IB)) 
   10 CONTINUE 
      GOTO 40 
C 
C  Case 2: interpolation produces a RVP cf for temperatures from 45-75 
C          degrees. 
C 
C      (a) First assign the endpoints at 45 and 75 degrees. 
C 
   20 CFB = 1.0 
      DO 30 IB = 1, 3 
      CFL = RG23(1,IB,IP,IGR)*RDIFA 
      CFL = EXP(CFL) 
C 
C      (b) Now interpolate between the 2 endpoints (CFB,CFL). 
C 
      TCF(IB) = CFB + (CFL - CFB)*(REAT - TLL(1)) / TLL(2) 
C 
   30 CONTINUE 
C 
   40 RETURN 
      END 
