      BLOCK DATA BD42 
C 
C  BLOCK DATA Subprogram 42: /STDIST/ 
C 
C  BLOCK DATA 42 initializes the weekday and weekend start distributions. 
C 
C  Changes: (Last change first) 
C 
C   10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C   23 May 00 @EPA-djb Added SDIST initialization. 
C    7 Feb 00 @DynTel-ddj 1-032, Modified STARTSPERDAY for vehicle type 
C             28 (LDDT34) to be equal to vehicle type 15(LDDT12). 
C   19 May 99 @DynTel-HXQ 2-650 Change the array STARTSPERDAY 
C   11 May 99 @DynTel-bag 2-650 STARTDIST now a fractions rather than percentages. 
C             This eliminates a multiplication by 0.01 in HRLOOP 
C    6 May 99 @DynTel-MLA 2-000 Expanded STARTDIST to include 24 intervals 
C             for the day instead of only 14 by dividing the last interval 
C             into 11 equal parts. 
C   26 Apr 99 @DynTel-bg 2-000 Explicit Typing 
C 
C  Common block array subscripts: 
C 
C  STARTDIST(14,2)    -  STARTDIST (IH,IWEEK) 
C  STARTSPERDAY(28,2)  -  STARTSPERDAY (IV,IWEEK) 
C 
      IMPLICIT NONE 
      INCLUDE 'STDIST.I' 
C 
      DATA SDIST / 0 / 
C 
C  /STDIST/: 
C 
C  Hardcoded start distribution data 
C 
        DATA STARTDIST/ 
C 
C  weekdays  
C  
     * 0.02040816327,0.05541601256,0.06028257457,0.04725274725, 
     * 0.05164835165,0.0671899529 ,0.08069073783,0.07299843014, 
     * 0.08037676609,0.08979591837,0.084144427  ,0.07723704867, 
     * 0.0601255887 ,11*0.013857571, 
C 
C  weekends  
C 
     * 0.00914169629,0.01929913662,0.03098019299,0.06449974606, 
     * 0.06907059421,0.07973590655,0.10157440325,0.07262569832, 
     * 0.08887760284,0.07364144236,0.08024377857,0.07110208228, 
     * 0.06145251397,11*0.01615956415/ 
C 
C  Hardcoded starts per day data 
C 
C  At present, we only have data for LDGV and LDGT. As a temporary 
C  measure, until other data is provided, we will use the LDGV value 
C  for LDDV and MC and we will use the LDGT value for LDDT. 
C 
        DATA STARTSPERDAY/ 
C 
C  weekdays  
C  
     *25*7.28, 100*8.06, 200*6.88, 25*7.28, 25*8.06, 200*6.65, 
     *25*1.35,  25*6.88,  50*6.65, 25*8.06, 
C 
C  weekends  
C 
     *25*5.41, 100*5.68, 200*6.88, 25*5.41, 25*5.68, 200*6.65, 
     *25*1.35,  25*6.88,  50*6.65, 25*5.68/ 
C 
      END 
