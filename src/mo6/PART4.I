C 
C  PART4 contains the  EPA-based Zero Mile Levels (ZML) 
c  and Deterioration Rates (DR) for diesels 
C 
C  Changes: (Last change first) 
C 
C  20 Dec 01: AIR Task 25: The associated diesel sulfur level is now by-model-year. 
C  16 Oct 01; AIR Task 25: Changed to handle only diesel values. 
C  28 Aug 01: Air Task 23c: Implement DR1 and DR2 rates. 
C  24 Aug 01: AIR Task 23b: EPA-based ZMLs/DRs 
C  05 Jun 01: AIR Task 17: Added CAAROMFLG. 
C  05 Jun 01: AIR Task 17: Added base sulfur level array. 
C  23 May 01: AIR Task 17: Changed from 1975 to 1950 base and only 2 modes. 
C  20 Mar 01: AIR Task 11: New INCLUDE. 
C 
C  PMDDR1(71,28) - PMDDR1(MY,IV) 
C  PMDDR2(71,28) - PMDDR2(MY,IV) 
C  PMDZML(71,28) - PMDZML(MY,IV) 
C  BSULFLVL(71)  - BSULFLVL(MY) 
C 
C  MY is model year, 1950-=1, 2020+=71 
C 
C  Common Block Dictionary: 
C 
C   Name      Type                  Description 
C  ------     ----  ----------------------------------------------------- 
C /PART4/ 
C   BSULFLVL  R     Base sulfur level, in ppm, of PM ZML/DR values for diesels 
C                   for each model year 
C   PMDDR1    R     First PM Deterioration Rate (g/mi/10K mi) for diesels 
C   PMDDR2    R     Second PM Deterioration Rate (g/mi/10K mi) for diesels 
C   PMDRAGE   I     Age at which PMDDR1 ends and PMDDR2 begins 
C   PMDZML    R     PM Zero Mile Level (g/mi) for diesels 
C 
C The data is initialized in BD44. 
C 
      INTEGER PMDRAGE 
      REAL    BSULFLVL 
      REAL    PMDDR1 
      REAL    PMDDR2  
      REAL    PMDZML 
C 
      COMMON /PART4/ PMDZML(71,28),PMDDR1(71,28), 
     *               PMDDR2(71,28),PMDRAGE(28),BSULFLVL(71) 
