C 
C  PART4 contains the  EPA-based Zero Mile Levels (ZML) 
c  and Deterioration Rates (DR). 
C 
C  Changes: (Last change first) 
C 
C  28 Aug 01: Air Task 23c: Implement DR1 and DR2 rates. 
C  24 Aug 01: AIR Task 23b: EPA-based ZMLs/DRs 
C  05 Jun 01: AIR Task 17: Added CAAROMFLG. 
C  05 Jun 01: AIR Task 17: Added base sulfur level array. 
C  23 May 01: AIR Task 17: Changed from 1975 to 1950 base and only 2 modes. 
C  20 Mar 01: AIR Task 11: New INCLUDE. 
C 
C  BSULFLVL(3)      - BSULFLVL (see below) 
C  PMDR1(71,2,28) - PMDR1(MY,ICAT,IV) 
C  PMDR2(71,2,28) - PMDR2(MY,ICAT,IV) 
C  PMZML(71,2,28) - PMZLM(MY,ICAT,IV) 
C 
C  MY is model year, 1950-=1, 2020+=71 
C 
C  ICAT is for catalyst equipment, as follows: 
C    1 = no catalyst, diesel 
C    2 = catalyst 
C 
C  Common Block Dictionary: 
C 
C   Name      Type                  Description 
C  ------     ----  ----------------------------------------------------- 
C /PART4/ 
C   BSULFLVL  R     Base sulfur level, in ppm, of PM ZML/DR values for 
C                   1=Gas CY1999-, 2=Gas CY2000+, 3=Diesel 
C   PMDR1     R     First PM Deterioration Rate (g/mi/10K mi) 
C   PMDR2     R     Second PM Deterioration Rate (g/mi/10K mi) 
C   PMDRAGE   I     Age at which PMDR1 ends and PMDR2 begins 
C   PMZML     R     PM Zero Mile Level (g/mi) 
C 
C The data is initialized in BD44. 
C 
      INTEGER PMDRAGE 
      REAL    BSULFLVL 
      REAL    PMDR1 
      REAL    PMDR2  
      REAL    PMZML 
C 
      COMMON /PART4/ PMZML(71,2,28),PMDR1(71,2,28), 
     *               PMDR2(71,2,28),PMDRAGE(28),BSULFLVL(3) 
