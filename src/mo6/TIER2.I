C 
C  Changes: (Last change first) 
C 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased 
C         size of T2_EVAP array 
C  2 Feb 00 @DynTel-ddj 1-032, Added array CRT50KT2. 
C 14 Jan 00 @Dyntel jwr  1-017 Added T2_EVAP(12,5) and T2_RED(15) 
C 30 Dec 99 @DynTel-ddj  1-021  New include FILE for Tier2 arrays. 
C 
C The initialization of the data could be found in BD_TIER2. 
C 
      REAL TIER2_50K, T2_BIN_FRAC, T2CRT50K, T2_EVAP, T2_RED 
      COMMON /TIER2/ TIER2_50K(5,12,3), T2_BIN_FRAC(12,12,3,5), 
     &               T2CRT50K(5,3,25), T2_EVAP(12,15), T2_RED(15) 
