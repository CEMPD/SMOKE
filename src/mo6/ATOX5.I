C 
C  Changes: (Last change first) 
C 
C 08 Mar 02: AIR Task 28: Moved ADDITIONAL HAPS related arrays to ATWORK. 
C 28 Nov 01: AIR Task 28: Moved ATHBER to ATWORK. 
C 04 Sep 01: AIR Task 18: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      INTEGER    NATH 
      INTEGER    NMAP 
      LOGICAL    ATHAPS 
C 
      COMMON /ATOX5/ NATH,NMAP,ATHAPS 
