C 
C  Changes: (Last change first) 
C 
C 09 Mar 02: AIR Task 28: Moved EFAT and VAT to ATWORK due 
C            to expanded ADD HAPS. 
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
C 
      COMMON /ATOX4/ ATCONV(2,2,6), ATNHWT(2), ATHEF(24,5,2), 
     *               ATEMFAC(2,10,3,2,2),ATHEV(24,5,5) 
