      SUBROUTINE TESTIM(INERR,RC) 
C 
C  TESTIM  checks the parameters that were entered for each I/M program. 
C 
C  Called by CHKINP 
C 
C  Calls QUITER, CKOBDIM 
C 
C  Changes: (Last change first) 
C 
C 20 Jun 01 @EPA-djb Bug 326, Removed call to QUITER for warning M521. 
C 15 Sep 00 AIR Task 08: Updated to check exhaust and evap I/M separately. 
C 07 Jul 00 AIR Task 05: Removed unused RC from CKOBDIM parameter list. 
C 22 May 00 @EPA-elg   Added call to CKOBDIM. 
C 16 Aug 99 @DynTel-bg  2-663 Renamed SIMFLG and SIMONLY to IMFLG_NP and 
C           IMONLY_NP, respectively. These are for the 'non-periodic' COIM 
C           and OBD programs. 
C 21 Jun 99 @DynTel-bg  2-705  Added a warning for user input of OBD I/M 
C           vehicle models less than 1994. Initial vehicle model year, 
C           MODYR(1,IMPGM), reset to 1994 for this case. 
C 16 Jun 99 @DynTel-MLA 2-663  Replaced INTYP(IMPGM) with IMCHECK(23,IMPGM) 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list:  
C    common blocks: 
C    /CHKLST/ IMCHECK,MAXIMPAR,MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /FLAGS2/ IMFLAG,EVIMFLAG 
C    /IMPAR1/ ICYIM,MODYR 
C    /IMPAR2/ ITEST,ILDT 
C    /IMPAR6/ IFREQ,INTYP 
C    /LDGSIM/ IMFLG_NP,IMONLY_NP 
C    /RSDINF/ RSDCLEAN,RSDDIRTY 
C    /TTCINF/ IASM 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------------- 
C  ALLFOUND   L   Tests for a complete sets of inputs. 
C  COVEROP1   L   Was RSD coverage option 1 selected?(The coverage or fraction of vehicles 
C                 measured by remote sensing is determined with a Poisson distribution, 
C                 the number of tests that will be performed, and the total number of 
C                 vehicles.) 
C  COVEROP2   L   Was RSD coverage option 2 selected?(The coverage or fraction of vehicles 
C                 measured by remote sensing is the user input.) 
C  COVEROP3   L   Was RSD coverage option 3 selected?(The user enters the total number of 
C                 vehicles and the number of confirmed remote sensing test failures.) 
C  CSIMPLIED  L   Was clean screening implied because more than one 
C                 RSD coverage option was selected? 
C  IMINPUT    I   Loop control variable. 
C  IMPGM      I   The I/M program number. 
C  IVIM       I   ILDT pointer 
C 
C  Notes: 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'RSDINF.I' 
      INCLUDE 'TTCINF.I' 
C 
      INTEGER        INERR 
      INTEGER        RC 
C 
      LOGICAL        ALLFOUND 
      LOGICAL        COVEROP1 
      LOGICAL        COVEROP2 
      LOGICAL        COVEROP3 
      LOGICAL        CSIMPLIED 
      INTEGER        IMINPUT 
      INTEGER        IMPGM 
      INTEGER        IVIM 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     IMFLAG, in this version of the model, simply counts the I/M 
C     programs defined. (The number of programs is actually IMFLAG-1.) 
C     Get a count of the I/M programs before we call CKOVLP.  Count 
C     exhaust and evap I/M programs separately. 
C 
      IMFLAG = 1 
      EVIMFLAG = 1 
C 
      DO IMPGM = 1, MAXIMPGM 
C 
C       The first item on the checklist tells whether any data has been 
C       entered for a specific I/M program.          
C 
        IF (IMCHECK(1,IMPGM).EQ.1) CYCLE 
C 
        IF(ITEST(IMPGM).LT.20) THEN 
          IMFLAG = IMFLAG + 1 
        ELSE 
          EVIMFLAG = EVIMFLAG + 1 
        ENDIF 
C 
      END DO 
C 
C     Call CKOBDIM to determine whether exhaust I/M exists for 
C     post-95 MYs or OBD I/M exists for pre-1996 MYs. 
C 
      CALL CKOBDIM(INERR) 
C 
C     Call CKOVLP to find out whether any of the programs overlap. 
C 
      CALL CKOVLP(INERR,RC) 
C 
C     Check the consistency of the rest of the I/M program parameters. 
C     Only do the following checks for exhaust I/M programs. 
C 
      IF(IMFLAG.GE.2) THEN 
C 
      DO IMPGM = 1, MAXIMPGM 
C 
        IF (IMCHECK(1,IMPGM).EQ.1 .OR. ITEST(IMPGM).GE.20) CYCLE 
C 
C       Look at the program data. First, check all the required items. 
C       Items 2 through 6 are entered on the I/M PROGRAM input record. 
C 
        ALLFOUND = .TRUE. 
        DO IMINPUT = 2,6 
          IF (IMCHECK(IMINPUT,IMPGM).EQ.1) THEN 
            ALLFOUND = .FALSE. 
            EXIT 
          ENDIF 
        END DO 
C 
        IF (.NOT.ALLFOUND) THEN 
          CALL QUITER(0.0,IMPGM,187,INERR) 
          RC = 2 
        ENDIF 
C 
C       Make sure cutpoints were entered if the prog type is IM240. 
C 
        IF (ITEST(IMPGM).EQ.4.AND.IMCHECK(7,IMPGM).EQ.1) THEN 
          CALL QUITER(0.0,IMPGM,188,INERR) 
          RC = 2 
        ENDIF 
C 
C       Check the other required inputs - model years covered, vehicles 
C       covered, stringency, compliance rate, and waiver rates. 
C 
        DO IMINPUT = 8,12 
          IF (IMCHECK(IMINPUT,IMPGM).EQ.1) THEN 
            CALL QUITER(1.0*(IMINPUT-7),IMPGM,189,INERR) 
            RC = 2 
          ENDIF 
        END DO 
C 
C       Check all of the RSD inputs. They are all optional, but they 
C       must all be consistent. First, check for coverage option 1. 
C 
        COVEROP1 = (IMCHECK(13,IMPGM).EQ.2.OR. 
     &              IMCHECK(14,IMPGM).EQ.2.OR. 
     &              IMCHECK(15,IMPGM).EQ.2) 
        ALLFOUND = (IMCHECK(13,IMPGM).EQ.2.AND. 
     &              IMCHECK(14,IMPGM).EQ.2.AND. 
     &              IMCHECK(15,IMPGM).EQ.2) 
C 
        IF (COVEROP1.AND..NOT.ALLFOUND) THEN 
          CALL QUITER(0.0,IMPGM,190,INERR) 
          RC = 2 
        ENDIF 
C 
C       Check for coverage option 2. 
C 
        COVEROP2 = (IMCHECK(16,IMPGM).EQ.2) 
C 
C       Make sure that not both options 1 and 2 were selected. 
C 
        IF (COVEROP1.AND.COVEROP2) THEN 
          CALL QUITER(0.0,IMPGM,501,INERR) 
          RC = 2 
        ENDIF 
C 
C       Check for coverage option 3. 
C 
        COVEROP3 = (IMCHECK(17,IMPGM).EQ.2.OR.IMCHECK(18,IMPGM).EQ.2) 
        ALLFOUND = (IMCHECK(17,IMPGM).EQ.2.AND.IMCHECK(18,IMPGM).EQ.2) 
C 
        IF (COVEROP3.AND..NOT.ALLFOUND) THEN 
          CALL QUITER(0.0,IMPGM,502,INERR) 
          RC = 2 
        ENDIF 
C 
C       It is an error if the I/M program type is RSD ONLY and there 
C       is no RSD program modeled! 
C 
        IF (IMCHECK(23,IMPGM).EQ.2.AND..NOT.COVEROP1.AND. 
     &      .NOT.COVEROP2.AND..NOT.COVEROP3) THEN 
          CALL QUITER(0.0,IMPGM,503,INERR) 
          RC = 2 
        ENDIF 
C 
C       Coverage options 1 and 2 must be used in conjunction with 
C       an RSD CO CUTPOINT record (for dirty screening) or an 
C       RSD NOX CUTPOINT record (for clean screening). (If coverage 
C       option 3 is selected for dirty screening or the RSD LOSSES 
C       record is provided for clean screening, coverage options 1 
C       or 2 should not also be selected.) 
C 
        IF ((COVEROP1.OR.COVEROP2).AND.IMCHECK(19,IMPGM).EQ.1.AND. 
     &      IMCHECK(20,IMPGM).EQ.1) THEN 
          CALL QUITER(0.0,IMPGM,504,INERR) 
          RC = 2 
        ENDIF 
C 
C       Coverage option 3 may be selected along with either coverage 
C       option 1 or option 2 if and only if both dirty and clean 
C       screening programs are modeled and the clean screening program 
C       uses coverage option 1 or 2 and the dirty screening program 
C       uses option 3. (The clean screening program is not allowed to 
C       use coverage option 3.) These choices will be assumed if more 
C       than one coverage option is selected. If both clean and dirty 
C       screening programs are modeled and coverage option 3 is not 
C       selected for the dirty screening program, both programs must 
C       use coverage option 1 or both programs must use option 2. 
C 
        CSIMPLIED = (COVEROP3.AND.(COVEROP1.OR.COVEROP2)) 
C 
C       The CO CUTPOINT input record is appropriate only when coverage 
C       options 1 or 2 are selected for the RSD dirty screening program. 
C 
        IF (IMCHECK(19,IMPGM).EQ.2.AND..NOT.COVEROP1.AND. 
     &      .NOT.COVEROP2) THEN 
          CALL QUITER(0.0,IMPGM,505,INERR) 
          RC = 2 
        ENDIF 
C 
        IF (IMCHECK(19,IMPGM).EQ.2.AND.COVEROP3) THEN 
          CALL QUITER(0.0,IMPGM,506,INERR) 
          RC = 2 
        ENDIF 
C 
C       The user has two options for specifying a clean screening 
C       program. They may enter values for the NOX CUTPOINT and for 
C       RSD clean screening coverage, in which case the model will 
C       calculate the clean screening credit losses, or they may 
C       supply their own values for the clean screening credit 
C       losses, in which case neither of the NOX CUTPOINT or the 
C       RSD clean screening coverage values are required. (But the 
C       user may still enter values for RSD coverage, which implies 
C       that a dirty screening program is also modeled and the RSD 
C       coverage values should be applied to that program.) 
C 
        IF (CSIMPLIED.AND.IMCHECK(20,IMPGM).EQ.1) THEN 
          CALL QUITER(0.0,IMPGM,507,INERR) 
          RC = 2 
        ELSE IF (IMCHECK(20,IMPGM).EQ.2) THEN 
          IF (.NOT.COVEROP1.AND..NOT.COVEROP2) THEN 
            CALL QUITER(0.0,IMPGM,508,INERR) 
            RC = 2 
          ENDIF 
        ENDIF 
C 
C       If the RSD LOSSES record is used, the RSD NOX CUTPOINT record 
C       should not be included and if coverage option 3 is selected 
C       for the dirty screening program, coverage options 1 or 2 
C       should not also be selected. 
C 
        IF (IMCHECK(21,IMPGM).EQ.2) THEN 
          IF (IMCHECK(20,IMPGM).EQ.2) THEN 
            CALL QUITER(0.0,IMPGM,509,INERR) 
            RC = 2 
          ENDIF 
          IF (CSIMPLIED) THEN 
            CALL QUITER(0.0,IMPGM,510,INERR) 
            RC = 2 
          ENDIF 
        ENDIF 
C 
C       Set IMFLG_NP, RSDDIRTY, and RSDCLEAN for this I/M program. The settings 
C       of RSDDIRTY and RSDCLEAN describe the coverage options the user has 
C       supplied for the RSD clean and dirty screening programs. 
C       RSD coverage options are described in M6.IM.001 pages 38-40 
C 
        IF (COVEROP3) THEN 
          RSDDIRTY(IMPGM) = 3 
        ELSE IF (COVEROP2.AND.IMCHECK(19,IMPGM).EQ.2) THEN 
          RSDDIRTY(IMPGM) = 2 
        ELSE IF (COVEROP1.AND.IMCHECK(19,IMPGM).EQ.2) THEN 
          RSDDIRTY(IMPGM) = 1 
        ELSE 
          RSDDIRTY(IMPGM) = 0 
        ENDIF 
C 
        IF (IMCHECK(21,IMPGM).EQ.2) THEN 
          RSDCLEAN(IMPGM) = 3 
        ELSE IF (COVEROP2.AND.IMCHECK(20,IMPGM).EQ.2) THEN 
          RSDCLEAN(IMPGM) = 2 
        ELSE IF (COVEROP1.AND.IMCHECK(20,IMPGM).EQ.2) THEN 
          RSDCLEAN(IMPGM) = 1 
        ELSE 
          RSDCLEAN(IMPGM) = 0 
        ENDIF 
C 
C       IMFLG_NP is set to 0 to show no COIM or RSD program, 1 to show RSD 
C       only, 2 to show COIM only, and 3 to show both RSD and COIM. 
C 
        IF (RSDCLEAN(IMPGM).GT.0.OR.RSDDIRTY(IMPGM).GT.0) THEN 
          IF (IMCHECK(22,IMPGM).EQ.2.OR. 
     &        IMCHECK(24,IMPGM).EQ.2) THEN 
            IMFLG_NP(IMPGM) = 3 
          ELSE 
            IMFLG_NP(IMPGM) = 1 
          ENDIF 
        ELSE 
          IF (IMCHECK(22,IMPGM).EQ.2.OR. 
     &        IMCHECK(24,IMPGM).EQ.2) THEN 
            IMFLG_NP(IMPGM) = 2 
          ELSE 
            IMFLG_NP(IMPGM) = 0 
          ENDIF 
        ENDIF 
C 
C       Write a warning message to the output file if the input designates 
C       selective I/M only (COIM ONLY or RSD ONLY). Set the grace period 
C       for the I/M program to 25 years (all model years are exempt from 
C       mandatory periodic inspections). (Not an error - only a warning.) 
C 
        IF (IMCHECK(23,IMPGM).EQ.2.OR.IMCHECK(24,IMPGM).EQ.2) THEN 
          CALL QUITER(0.0,IMPGM,517,INERR) 
          IMGRACE(IMPGM) = 25 
          IMONLY_NP(IMPGM) = 2 
        ENDIF 
C 
C       It's an error if the user indicates selective I/M (RSD only) 
C       and there is no RSD program defined. 
C 
        IF (IMCHECK(23,IMPGM).EQ.2.AND.IMFLG_NP(IMPGM).NE.1.AND. 
     &      IMFLG_NP(IMPGM).NE.3) THEN 
          CALL QUITER(0.0,IMPGM,518,INERR) 
          RC = 2 
        ENDIF 
C 
C       It's an error if the user indicates selective I/M (COIM only) 
C       and there is no COIM program defined. 
C 
        IF (IMCHECK(24,IMPGM).EQ.2.AND.IMFLG_NP(IMPGM).NE.2.AND. 
     &      IMFLG_NP(IMPGM).NE.3) THEN 
          CALL QUITER(0.0,IMPGM,519,INERR) 
          RC = 2 
        ENDIF 
C 
C       Issue a warning if an OBD I/M test is specified and the I/M 
C       vehicles model years include pre-1994 vehicles. Set the first 
C       covered model to 1994. 
C 
        IF ( ITEST(IMPGM).EQ.12.AND.MODYR(1,IMPGM).LT.1994) THEN 
          CALL QUITER(0.0,IMPGM,523,INERR) 
          MODYR(1,IMPGM) = 1994 
        ENDIF 
C 
      END DO 
C 
      END IF 
C 
C     Now check the consistency of the the evap I/M program parameters. 
C 
      IF(EVIMFLAG.GE.2) THEN 
C 
      DO IMPGM = 1, MAXIMPGM 
C 
        IF (IMCHECK(1,IMPGM).EQ.1 .OR. ITEST(IMPGM).LT.20) CYCLE 
C 
C       Look at the program data. First, check all the required items. 
C       Items 2 through 6 are entered on the evap I/M PROGRAM input record. 
C 
        ALLFOUND = .TRUE. 
        DO IMINPUT = 2,6 
          IF (IMCHECK(IMINPUT,IMPGM).EQ.1) THEN 
            ALLFOUND = .FALSE. 
            EXIT 
          ENDIF 
        END DO 
C 
        IF (.NOT.ALLFOUND) THEN 
          CALL QUITER(0.0,IMPGM,559,INERR) 
          RC = 2 
        ENDIF 
C 
C       Make sure that the required model year inputs were entered. 
C 
        IF (IMCHECK(8,IMPGM).EQ.1) THEN 
          CALL QUITER(1.,IMPGM,560,INERR) 
          RC = 2 
        ENDIF 
C 
C       Make sure that at least one evap vehicle class was entered 
C 
        IF (IMCHECK(9,IMPGM).EQ.1) THEN 
          CALL QUITER(2.,IMPGM,560,INERR) 
          RC = 2 
        ELSE 
          DO IVIM=1,14 
          IF(ILDT(IVIM,IMPGM).GT.1) GOTO 20 
          END DO 
        CALL QUITER(2.,IMPGM,560,INERR) 
        RC = 2 
   20   CONTINUE 
        ENDIF 
C 
      END DO 
C 
      END IF 
C 
   99 RETURN 
      END 
