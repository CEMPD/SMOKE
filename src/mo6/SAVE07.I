C 
C  Changes: (Last change first) 
C  30 Apr 01 AIR Task 12: Added VA022 for T2PM_FRAC. 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased 
C         size of VA016A array 
C  14 Jan 00 @Dyntel jwr 1-017 Added VA016A for T2_EVAP 
C  30 Dec 99 @DynTel-ddj 1-021, Added VA016 and VA017 for T2_BIN_FRAC 
C            and TIER2_50K respectively. 
C  22 Jun 99 @DynTel-MLA 2-694  Moved VA015 here (from BASEQ9) and removed 
C            VA015 and VA016 (for obsolete P94REF and P94SEF). 
C 
C  The data is initialized in BDSAVE.FOR 
C 
      REAL VA015, VA016, VA017, VA016A, VA022 
      COMMON /SAVE07/ VA015(11,32,5), VA016(12,12,3,5), VA017(5,12,3) 
     & , VA016A(12,15), VA022(12,12,8) 
