C 
C These arrays hold the results of the soak time correction factor  
C calculations. 
C 
C  Changes: (Last change first) 
C 
C  21 May 99 @DynTel-ddj, 2-700  Added array RATIO_A(3). 
C   3 May 99 @DynTel-bag  2-000 MINUTES now an array rather than determining 
C            each value in SOAKFCN for each ST_INDX and IP 
C  22 Oct 98 @DynTel-bag  Restructured soak time corrections calculations. 
C            Moved fit parameters from subroutines to SOAK common blocks, 
C            removed SOAKCF array, and reduced SOAKFUN to a two dimensional 
C            array dependent on soak time and catalyst type. 
C 
C 
C They are initialized in block data 41. 
C 
      INTEGER UPR_RANG_LO_LMT 
      REAL RATIO, RATIO_A 
      REAL CSA, CSB, CSC, SOAKFUN, SOK_TIME 
      COMMON /SOAK/ CSA(2,3,2), CSB(2,3,2),    CSC(2,3,2), 
     &              RATIO(3),   SOAKFUN(70,2), SOK_TIME(70), 
     &              RATIO_A(3), UPR_RANG_LO_LMT(3,2) 
