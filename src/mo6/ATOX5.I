C 
C  Changes: (Last change first) 
C 
C 28 Nov 01: AIR Task 28: Moved ATHBER to ATWORK. 
C 04 Sep 01: AIR Task 18: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      CHARACTER  ATHEKIND*4 
      CHARACTER  ATHPNAME*32 
      CHARACTER  ATHRATIO*8 
      CHARACTER  ATHSNAME*12 
      INTEGER    ATHBEGMY 
      INTEGER    ATHENDMY 
      INTEGER    ATHPMAP 
      INTEGER    ATHPOL 
      INTEGER    NATH 
      INTEGER    NMAP 
      LOGICAL    ATHAPS 
      LOGICAL    ATHIV 
      REAL       ATHDET 
      REAL       ATHRATZML 
      REAL       ATHVTYPE 
C 
      COMMON /ATOX5/ ATHBEGMY(50),ATHDET(50),ATHENDMY(50),ATHPOL(50), 
     &               ATHRATZML(50),ATHVTYPE(50),NATH,ATHEKIND(50), 
     &               ATHPMAP(50),NMAP,ATHRATIO(50),ATHPNAME(50), 
     &               ATHSNAME(50),ATHAPS,ATHIV(2,4,10,28) 
