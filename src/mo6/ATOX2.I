C 
C  Changes: (Last change first) 
C 
C 22 Jun 01: AIR Task 15: New INCLUDE. 
C 
C The data is initialized in BD47. 
C 
      REAL ATACET 
      REAL ATACETBASE 
      REAL ATACROLEIN 
      REAL ATBASE 
      REAL ATBENZ 
      REAL ATBENZBASE 
      REAL ATBUTA 
      REAL ATBUTABASE 
      REAL ATCENTER 
      REAL ATCENTERA 
      REAL ATFORM 
      REAL ATFORMBASE 
      REAL ATMTBE 
      REAL ATTRGT 
      REAL ATVOC 
      REAL ATVOCBASE 
      REAL ATWTOXY 
      REAL CMHC 
      REAL VOLTOWT 
C 
      COMMON /ATOX2/ ATACET(7,11), ATACETBASE(10), 
     &               ATBASE(3,13), ATTRGT(13), ATWTOXY, 
     &               ATBENZ(6,11), ATBENZBASE(10), 
     &               ATBUTA(6,11), ATBUTABASE(10), 
     &               ATCENTER(14), ATCENTERA(14), 
     &               ATFORM(7,11), ATFORMBASE(10), 
     &               ATMTBE(4,10), 
     &               ATVOC(10,11), ATVOCBASE(10), 
     &               ATACROLEIN(28,2), 
     &               VOLTOWT(4), CMHC(10) 
