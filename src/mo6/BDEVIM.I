C 
C  Changes: (Last change first) 
C 
C 15 Sep 00 AIR Task 08: New include.I 
C 
C The data is initialized in BDEVIM. 
C 
      REAL TEST_GC 
      INTEGER MYB_GC 
C 
      REAL TEST_FP 
      INTEGER MYB_FP 
C 
      COMMON/BDEVIM/TEST_GC(16),MYB_GC, 
     *              TEST_FP(18),MYB_FP 
