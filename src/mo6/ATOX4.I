C 
C  Changes: (Last change first) 
C 
C 28 Nov 01: AIR Task 28: Moved DBAT to ATWORK Module. 
C 10 Jul 01: AIR Task 15: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      REAL ATCONV 
      REAL ATEMFAC 
      REAL ATHEF 
      REAL ATHEV 
      REAL ATNHWT 
      REAL EFAT 
      REAL VAT 
C 
      COMMON /ATOX4/ ATCONV(2,2,6), ATNHWT(2), ATHEF(24,5,2), 
     *               ATEMFAC(2,10,3,2,2),ATHEV(24,5,5), 
     *               EFAT(6,16,28), VAT(6,16) 
