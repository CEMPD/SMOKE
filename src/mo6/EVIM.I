C 
C  Changes: (Last change first) 
C 
C 15 Sep 00 AIR Task 08: New include file. 
C 
C  The data is initialized in BDEVIM. 
C 
      REAL FAIL_IL, FAIL_AZFP, FAIL_AZGC 
      REAL EVNOIM, EVIM 
      REAL EVSAW 
      REAL GLLRDC 
C 
      COMMON/EVIM/ FAIL_IL(3,26), FAIL_AZFP(3,26), FAIL_AZGC(2,26), 
     *             EVNOIM(4,26,2), EVIM(4,26,2), 
     *             EVSAW(4,25,2), GLLRDC 
       
