C 
C  Changes: (Last change first) 
C 
C 04 Sep 01: AIR Task 18: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      CHARACTER  ATHEKIND*4 
      CHARACTER  ATHPNAME*32 
      CHARACTER  ATHRATIO*8 
      INTEGER    ATHBEGMY 
      INTEGER    ATHENDMY 
      INTEGER    ATHPMAP 
      INTEGER    ATHPOL 
      INTEGER    NATH 
      LOGICAL    ATHAPS 
      REAL       ATHBER 
      REAL       ATHDET 
      REAL       ATHRATZML 
      REAL       ATHVTYPE 
C 
      COMMON /ATOX5/ ATHBEGMY(50),ATHDET(50),ATHENDMY(50),ATHPOL(50), 
     &               ATHRATZML(50),ATHVTYPE(50),NATH,ATHEKIND(50), 
     &               ATHPMAP(50),ATHBER(2,25,2,4,28,50), 
     &               ATHRATIO(50),ATHPNAME(50),ATHAPS 
