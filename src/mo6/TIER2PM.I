C 
C  Changes: (Last change first) 
C 
C  27 Apr 01 AIR Task 12: New INCLUDE 
C 
C   Common block array dictionary: 
C 
C   Name          Type              Description 
C  ------         ----  --------------------------------------------- 
C  T2PM_FRAC  R    Tier2 phase-in bin fractions for PM. 
C  T2PM_50K    R    PM 50K Certification Standards for LDGV, LDGT1, 
C                       LDGT2, LDGT3, LDGT4, LDDV, LDDT12 and LDDT34 
C                       each of 12 bins. 
C 
C The initialization of the data is found in BD_TIER2PM. 
C 
      REAL T2PM_FRAC 
      REAL T2PM_50K 
C 
      COMMON /T2PM/ T2PM_50K(8,12), 
     *              T2PM_FRAC(12,12,8) 
