       REAL FUNCTION EVWTPRE(MY,ACTUAL_MY,JKAGE,ITST,HOT_DIURN_RUN ) 
C 
C      EVWTPRE calculates the sawtooth averaged weighting factors  
C      for model years up to and including 1995 
C      calling the pre-95 function SAWLOG  
C      It then applies the gross liquid 
C      leaker fractions calling the function GLL_WGT. 
C      The evap I/M adjustments are included. 
C 
C     The values for test-status are 
C     1=pass pressure, pass purge 
C     2=pass pressure, fail purge 
C     3=fail pressure pass purge 
C     4=fail pressure and purge 
C 
C      The gross leaker fraction will be obtained by setting ITST = 5 
C 
C      Called by RSTLOS  
C 
C      Calls GLL_WGT 
C 
C  Changes: (Last change first) 
C 
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch 
C  15 Sep 00 AIR Task 08: Replaced calculations with array lookup 
C            from EVAPCALC; moved calculations done in SAWLOG to 
C            EVAPCALC called outside model year loop in HCCALX. 
C  18 Jul 00 AIR Task 06: Added code to allow substitution of  
C            user-supplied Fixed P/P Rates. Also, fixed outdated 
C            coding which left unitialized values and unused arrays. 
C 
C      Input on call: 
C 
C      Argument list 
C      IICY  (calendar year) 
C      JKAGE (age) 
C      ITST  (test-status =1-4 for pressure/purge, =5 for GLL) 
C      HOT_DIURN_RUN . Steering index 1=HOTSOAK, 
C      2=DIURNAL, 3=RUNNINGL 
C      /EVAPACT/ FPPRUSER, FPPRFLAG 
C 
C      Output: 
C      Function value 
C 
C      Local array/variable dictionary 
C 
C      Name        Type          Description 
C 
C      ANS          R            holds the weighting factors 
C                                while the liquid leaker fraction is 
C                                applied 
C      GLL_PREF     R            holds 1- gros liq lea fracn 
C 
C      JPHASE       I             age label for phase-in factors 
C      MY           I             Vehicle model year 
C 
C 
       USE LEAKYDAT, ONLY : GLL_COEF,PHASE_IN 
       IMPLICIT NONE 
C 
       INCLUDE 'EVAPACT.I' 
       INCLUDE 'EVIM.I' 
C 
C      Declare parameters 
C 
       INTEGER, INTENT(IN) :: ITST 
       INTEGER, INTENT(IN) :: JKAGE 
       INTEGER, INTENT(IN) :: ACTUAL_MY 
       INTEGER, INTENT(IN) :: MY 
       INTEGER, INTENT(IN) :: HOT_DIURN_RUN 
C 
C      Declare external functions 
C        
       REAL GLL_WGT 
C 
C      Declare local variables 
C 
       REAL ANS 
       REAL GLL_PREF 
C 
C      If we only want the leaker fraction bypass the  
C      pressure/purge calculation 
C 
       IF(ITST.EQ.5) GO TO 200 
C 
C      Otherwise begin the pressure/purge calculation by calling 
C      the relevant sawtooth routine according to the model yr. 
C      However, when requested, apply the user-supplied  
C      Fixed P/P Rates to all model years and ages. 
C 
       IF(FPPRFLAG.EQ.1) THEN 
         ANS=FPPRUSER(ITST) 
       ELSE 
         IF(MY.LE.1998) THEN 
           ANS=EVSAW(ITST,JKAGE,1) 
         ELSE 
           ANS=0.0 
         END IF 
       END IF 
C 
C Reduce by the gross liquid leaker factor, adjusted by user, if required 
C 
       GLL_PREF=1.0-(GLL_WGT(MY,ACTUAL_MY,HOT_DIURN_RUN,JKAGE-1) 
     *         *(1.0-GLLRDC)) 
C 
       ANS=ANS*GLL_PREF 
C 
 200   CONTINUE 
C 
C Return the function value  
C 
C For pressure purge strata 
C 
      IF(ITST.GE.1 .AND. ITST.LE.4) THEN 
C  
        EVWTPRE = ANS 
C       
C For the gross liquid leakers 
C 
      ELSE IF(ITST.EQ.5) THEN 
C 
        EVWTPRE=GLL_WGT(MY,ACTUAL_MY,HOT_DIURN_RUN,JKAGE-1)*(1.0-GLLRDC) 
C 
      END IF 
C 
      RETURN 
C 
      END 
