      REAL FUNCTION EF_POST93(MY,IP,IV,VMTAGE,ISR,INH,ISTD) 
C 
C  EF_POST93 returns the normal and high emission factors for 1994+ light 
C  duty gasoline vehicles and trucks. 
C 
C  Called by ARBEF, EF_DSPTCH, EF_HI_FRCN, IMPROG, IMBENEFIT 
C 
C  Calls EF81_93 
C 
C  Changes: (Last change first) 
C 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C  26 Sep 00 Bug 64: Modified HC Start/Running Splits (JK) 
C  07 Jul 00 AIR Task 05: Removed unused JDX from parameter list. 
C  15 Jun 00 @EPA-bag, Bug73 Removed tier0 requirement for application 
C                 of CRT50 ratio to emissions which was missed in 
C                 M6000614 integration 
C   8 Jun 00 @EPA-ddj, Bug73 Applied CRT50 ratio to emissions. 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C   9 MAR 00 @EPA-elg           HRLOOP updates 
C   7 Jan 00 @DynTel-ddj 1-024, Modified Tier1 and LEV EF_POST93 
C            equations for HC and NOx. 
C  25 Oct 99 @DynTel-bag 2-694 Fixed Tier0 emissions call to ef81_93 for HC. Rather than 
C            having HC in the parameter list CO was in the parameter list. 
C  12 Oct 99 @DynTel-bag 2-694 Added parameters for ISTD and vehicle type; added comments. 
C  8  Sep 99 @DynTel-bag 2-694 Renamed to EF_POST93 
C  26 AUG 99 @DynTel-HXQ Deleted the non-used flag IGM 
C  28 Jul 99 @DynTel-bg  2-694 Comments, rearranged variables declaration list. 
C  27 Jul 99 @DynTel-MLA 2-694 Added code to test for user input. 
C  29 May 99 @DynTel-ddj 2-700,  
C  18 May 99 @Dyntel-HXQ 2-650 delete the code of changing unit of basic 
C                        emission from gram/start to gram/mile 
C  29 Apr 99 @DynTel-ddj 2-694, Updated regression coefficients 
C            to reflect revised EPA documentation.  Also streamlined 
C            coding. 
C  22 Feb 99 @DynTel-RJD 2-694 New function for Mobile 6. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,IP,IV,VMTAGE,ISR,INH,ISTD 
C    common blocks: 
C 
C  /CEVBMY/ BMYMPD,BMYTPD 
C  /LDGBER/ CRT050K,CRT100K 
C 
C  Output on return: 
C 
C    function: EF_POST93 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C  SLOPE1   R    slope of ef curve 
C  ZPOINT   R    zero mileage point = intercept of ef curve 
C  KMILES   R    vehicle miles traveled in kilomiles. 
C 
      IMPLICIT NONE 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LDGBER.I' 
C 
C     Declare external functions 
C 
      REAL EF81_93 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: INH 
      INTEGER, INTENT(IN) :: ISTD 
      REAL   , INTENT(IN) :: VMTAGE 
C 
C     Declare local variables 
C 
      REAL     CERTRATIO 
      INTEGER, PARAMETER :: RUN=1,START=2 
      INTEGER, PARAMETER :: NORMAL=1,HIGH=2 
      INTEGER, PARAMETER :: HC=1,CO=2,NOX=3 
      INTEGER, PARAMETER :: TIER0=1,TIER1=3,TIER2=4 
      INTEGER, PARAMETER :: MY1993=1993 
      INTEGER, PARAMETER :: PORTED_FI=3 
C 
      REAL    EF_FTP 
      REAL    KMILES 
      REAL    USEREF  ! user supplied emissions passed through CHK94UEF 
      LOGICAL USERINP ! determines whether or not user emissions have been supplied 
C 
      CERTRATIO = 1.0 
C 
C Test for BD40 block data input or user input 
C 
      KMILES = VMTAGE * 10.0 
      CALL CHK94UEF(IV,IP,MY,ISR,INH,ISTD,KMILES,USERINP,USEREF) 
      IF (USERINP) THEN 
        EF_POST93 = USEREF 
        GOTO 95 
      ENDIF 
C 
C Calculate tier 0 basic emission factors for NMHC emissions 
C 
      POLLUTANT : SELECT CASE (IP) 
C 
      CASE(HC) POLLUTANT 
C 
        IF ( ISTD.EQ.TIER0 ) THEN 
C 
C       Tier 0 HC basic emission factors are calculated from EF81_93 
C       with the settings IFDS=3(ported fuel injection) and MY = 1993. 
C 
          EF_POST93 = EF81_93(MY1993,HC,IV,VMTAGE,PORTED_FI,ISR,INH) 
C 
        ELSE 
C 
C Find the certification ratio for HC emissions 
C 
         IF(TIER2_FLAG.NE.1 .OR. MY.LE.2003 .OR. VVLDGAS(IV).EQ.0 
     &      .OR. ISTD.NE.TIER2) THEN 
C                                  see eqs.11a&b, p.10 M6.EXH.007 
            CERTRATIO = CRT050K(ISTD,HC,IV)/CRT050K(TIER0,HC,VTLDGV) 
          END IF 
C 
C         The following normal and high tier0 vehicle basic emissions 
C         lines in grams/mile are used as a basis for tier1 and later 
C         basic emissions lines. Ratios of appropriate certification 
C         standards to the tier0 standard of 0.41grams/mile (from BD40.FOR) 
C         are used to determine the emissions for all standards 
C         classifications (see BD40.FOR and p.12 of M6.EXH.007). 
C 
          IF (INH.EQ.NORMAL) THEN 
             EF_FTP = 0.16 + 0.0186*VMTAGE    ! eq.13, p.19 M6.EXH.007 
          ELSE 
             EF_FTP = 2.076                   ! eq.14, p.19 M6.EXH.007 
          END IF 
C 
C         Running emissions for Tier 1 and LEV vehicles (gm/mile) 
          IF (ISR.EQ.RUN) THEN 
             EF_POST93 = (0.25+0.066*VMTAGE-0.0032*VMTAGE**2 
     *                    +0.00006*VMTAGE**3)*EF_FTP 
          ELSE ! START Emission 
             EF_POST93 = (10.75-0.95*VMTAGE+0.0474*VMTAGE**2 
     *                    -0.0008*VMTAGE**3)*EF_FTP 
          END IF 
C 
        END IF 
C 
C ******** END IP = HC 
C 
      CASE(CO) POLLUTANT 
C 
        IF (ISTD.EQ.TIER0) THEN 
C 
C       Tier 0 CO basic emission factors are calculated from EF81_93 
C       with the settings IFDS=3(ported fuel injection) and MY=1993. 
C 
          EF_POST93 = EF81_93(MY1993,CO,IV,VMTAGE,PORTED_FI,ISR,INH) 
C 
        ELSE 
C 
         IF(TIER2_FLAG.NE.1 .OR. MY.LE.2003 .OR. VVLDGAS(IV).EQ.0 
     &      .OR. ISTD.NE.TIER2) THEN 
C 
C         Find the certification ratio for CO emissions 
C 
             CERTRATIO = CRT050K(ISTD,CO,IV)/CRT050K(TIER1,CO,IV) 
          END IF 
C 
C         Calculate Tier 1 and later FTP emission factors for CO pollutant 
C 
C         First do Tier1 through LEV standards. 
C 
C         HIGH distinguished from NORMAL emissions  
          IF (INH.EQ.NORMAL) THEN 
            IF (ISTD.LE.8) THEN 
              IF (ISR.EQ.RUN.AND.IV.LE.2) THEN         !running normal emissions cars and trucks1 
                EF_POST93 = 0.2821 + 0.2293*VMTAGE     ! eq.5a p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.3.OR.IV.EQ.4) ) THEN  !running normal emissions trucks2&3 
                EF_POST93 = 0.3219 + 0.2678*VMTAGE     ! eq.5b p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.5) ) THEN !running normal emissions trucks4 
                EF_POST93 = 0.366 + 0.2678*VMTAGE     ! eq.5b p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.IV.LE.2) THEN  !start normal emissions cars and trucks1 
                EF_POST93 = 15.176 + 0.0703*VMTAGE     ! eq.5c p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.(IV.EQ.3.OR.IV.EQ.4) ) THEN   !start normal emissions trucks2,3,&4 
                EF_POST93 = 21.884 + 0.168*VMTAGE      ! eq.5d p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.(IV.EQ.5) ) THEN  !start normal emissions trucks2,3,&4 
                EF_POST93 = 24.868 + 0.168*VMTAGE      ! eq.5d p.11, M6.EXH.009 
C 
              END IF 
C 
            ELSE !ultra-low and zero emissions vehicles (ISTD.GE.8) 
              IF (ISR.EQ.RUN.AND.IV.LE.2) THEN         !running normal emissions cars and trucks1 
                EF_POST93 = 0.1411 + 0.2293*VMTAGE     ! eq.6a p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.3.OR.IV.EQ.4)) THEN  !running normal emissions trucks2,3,&4 
                EF_POST93 = 0.1609 + 0.2678*VMTAGE     ! eq.6b p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.5)) THEN  !running normal emissions trucks2,3,&4 
                EF_POST93 = 0.183 + 0.2678*VMTAGE     ! eq.6b p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.IV.LE.2) THEN  !start normal emissions cars and trucks1 
                EF_POST93 = 7.588 + 0.0703*VMTAGE      ! eq.6c p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.(IV.EQ.3.OR.IV.EQ.4)) THEN   !start normal emissions trucks2,3,&4 
                EF_POST93 = 10.942 + 0.168*VMTAGE      ! eq.6d p.11, M6.EXH.009 
C 
              ELSE IF (ISR.EQ.START.AND.(IV.EQ.5)) THEN  !running normal emissions trucks2,3,&4 
                EF_POST93 = 12.434 + 0.168*VMTAGE        ! eq.6b p.11, M6.EXH.009 
C 
              END IF 
            END IF 
C 
C 
C         HIGH distinguished from NORMAL emissions  
          ELSE 
C 
            IF (ISR.EQ.RUN.AND.IV.LE.2) THEN           !running high emissions cars and trucks1 
              EF_POST93 = 36.106                       ! eq.7a p.12, M6.EXH.009 
C 
            ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.3.OR.IV.EQ.4)) THEN   !running high emissions trucks2 & 3 
              EF_POST93 = 33.283                       ! eq.7b p.12, M6.EXH.009 
C 
            ELSE IF (ISR.EQ.RUN.AND.(IV.EQ.5)) THEN    !running normal emissions trucks2,3,&4 
              EF_POST93 = 33.283                       ! eq.6b p.11, M6.EXH.009 
C 
            ELSE IF (ISR.EQ.START.AND.IV.LE.2) THEN    !start high emissions cars and trucks1 
              EF_POST93 = 38.060                       ! eq.7c p.12, M6.EXH.009 
C 
            ELSE IF (ISR.EQ.START.AND.(IV.EQ.3.OR.IV.EQ.4)) THEN  !start start emissions trucks2,3,&4 
              EF_POST93 = 83.862                       ! eq.7d p.12, M6.EXH.009 
C 
            ELSE IF (ISR.EQ.START.AND.(IV.EQ.5)) THEN  !start start emissions trucks2,3,&4 
              EF_POST93 = 83.862                       ! eq.7d p.12, M6.EXH.009 
C 
            END IF 
C 
          END IF 
C         HIGH distinguished from NORMAL emissions  
        END IF 
C       Tier 0 distinguished from Tier 1 and later FTP emission factors 
C 
C ******** END IP = CO 
C 
      CASE(NOX) POLLUTANT 
C 
C 
        IF (ISTD.EQ.TIER0) THEN 
C 
C       Tier 0 CO and NOX basic emission factors are calculated from EF81_93 
C       with the settings IFDS=3(ported fuel injection) and MY=1993. 
C 
          EF_POST93=EF81_93(MY1993,NOX,IV,VMTAGE,PORTED_FI,ISR,INH) 
C 
        ELSE 
C 
         IF(TIER2_FLAG.NE.1 .OR. MY.LE.2003 .OR. VVLDGAS(IV).EQ.0 
     &      .OR. ISTD.NE.TIER2) THEN 
C 
C       Calculate the certification ratio for NOX emissions 
C 
            CERTRATIO = CRT050K(ISTD,NOX,IV)/CRT050K(TIER1,NOX,VTLDGV)  !see eqs.11a&b, p.10 M6.EXH.007 
         END IF 
C 
C       Calculate Tier 1 FTP emission factors for NOx pollutant 
C 
          IF(INH .EQ. NORMAL) THEN 
             EF_FTP = 0.153 + 0.02941*VMTAGE   ! see eq.5 p.8, M6.EXH.007 
          ELSE  ! HIGH Emitter 
             EF_FTP = 1.294                    ! see eq.6 p.9, M6.EXH.007 
          END IF 
C 
          IF(ISR .EQ. RUN) THEN 
             EF_POST93 = 0.9*EF_FTP            ! see eq.23 p.27, M6.EXH.007 
          ELSE  ! START Emission 
             EF_POST93 = 1.37*EF_FTP           ! see eq.26 p.28, M6.EXH.007 
C                                NOTE:  Equation 26 gives 1.27 as the 
C                                       coefficient, however, Table B-2 
C                                       used 1.37.  EPA stated that 1.37 
C                                       is the correct value. 
          END IF 
C 
        END IF 
C 
C ******** END IP = NOX 
C 
      END SELECT POLLUTANT 
C 
C     Adjust emissions for FSVs using certification standards 
C 
C     FIXED ISTD.GT.2 BUG  03/09/00 ELG 
C 
      IF (IP.NE.CO .AND. INH.EQ.HIGH .AND. ISTD.GT.1) THEN 
        EF_POST93 = 0.5*(1.0 + CERTRATIO)*EF_POST93 
      ELSE 
        EF_POST93 = CERTRATIO * EF_POST93 
      END IF 
C 
   95 CONTINUE 
C 
      RETURN 
      END 
