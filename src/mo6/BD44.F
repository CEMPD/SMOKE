      BLOCK DATA BD44 
C 
C  EPA-based Zero Mile Levels (ZML) and Deterioration Rates (DR1 and DR2). 
C 
C  Changes: (Last change first) 
C 
C  20 Dec 01: AIR Task 25: The associated diesel sulfur level is now by-model-year. 
C  20 Oct 01: AIR Task 25: Gas Carbon (from PART5) and Total Diesel PM. 
C  28 Aug 01: Air Task 23c: Implemented DR1 and DR2 rates and PMDRAGE. 
C  24 Aug 01: AIR Task 23b: EPA-based ZMLs/DRs 
C  05 Jun 01: AIR Task 17: Added CAAROMFLG. 
C  05 Jun 01: AIR Task 17: Added base sulfur level array. 
C  23 May 01: AIR Task 17: Changed from 1975 to 1950 base and only 2 modes. 
C  20 Mar 01: AIR Task 11: New INCLUDE. 
C 
C  PMGZML(71,4,4) - PMGZML(MY,ICAT,IV) 
C 
C  BSULFLVL(3)     - BSULFLVL (see below) 
C  PMDDR1(71,28)   - PMDDR1(MY,IV) 
C  PMDDR2(71,28)   - PMDDR2(MY,IV) 
C  PMDRAGE(28)    - PMDRAGE(IV) 
C  PMDZML(71,28)   - PMZLM(MY,IV) 
C  PMGDR1(71,4,28) - PMGDR1(MY,ICAT,IV) 
C  PMGDR2(71,4,28) - PMGDR2(MY,ICAT,IV) 
C  PMGZML(71,4,28) - PMGZLM(MY,ICAT,IV) 
C  BSULFLVL(71)    - BSULFLVL(MY) 
C 
C  MY is model year, 1950-=1, 2020+=71 
C 
C  Common Block Dictionary: 
C 
C   Name      Type                  Description 
C  ------     ----  ----------------------------------------------------- 
C /EFPMC/: 
C   PMGZML  R    The carbon ZML emission factor, g/mi for gasoline 
C                vehicles referenced by model year, vehicle type 
C                and ICAT:  
C                  1 = Leaded Fuel (regardless of Cat or Air) 
C                  2 = Unleaded Fuel, Catalyst, No Air Pump 
C                  3 = Unleaded Fuel, Catalyst, Air Pump 
C                  4 = Unleaded Fuel, No Catalyst 
C   PMGDR1  R    The first carbon DR emission deterioration, g/mi/10k mi,  
C                for gasoline vehicles referenced by model year, 
C                vehicle type and ICAT:  
C                  1 = Leaded Fuel (regardless of Cat or Air) 
C                  2 = Unleaded Fuel, Catalyst, No Air Pump 
C                  3 = Unleaded Fuel, Catalyst, Air Pump 
C                  4 = Unleaded Fuel, No Catalyst 
C   PMGDR2  R    The second carbon DR emission deterioration, g/mi/10k mi,  
C                for gasoline vehicles referenced by model year, 
C                vehicle type and ICAT:  
C                  1 = Leaded Fuel (regardless of Cat or Air) 
C                  2 = Unleaded Fuel, Catalyst, No Air Pump 
C                  3 = Unleaded Fuel, Catalyst, Air Pump 
C                  4 = Unleaded Fuel, No Catalyst 
C   NCAT    I    Number of catalyst/fuel categories 
C 
C /PART4/ 
C   BSULFLVL  R     Base diesel sulfur level, in ppm, of PM ZML/DR for each model year 
C   PMDDR1    R     First PM Deterioration Rate (g/mi/10K mi) 
C   PMDDR2    R     Second PM Deterioration Rate (g/mi/10K mi) 
C   PMDRAGE  I     Age at which PMDDR1 ends and PMDDR2 begins 
C   PMDZML    R     PM Zero Mile Level (g/mi) 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EFGAS.I' 
      INCLUDE 'PART4.I' 
C 
      DATA PMGZML/7952*0./ 
      DATA PMGDR1/7952*0./ 
      DATA PMGDR2/7952*0./ 
      DATA NCAT/4/ 
C 
      DATA BSULFLVL/71*0./ 
      DATA PMDDR1/1988*0./ 
      DATA PMDDR2/1988*0./ 
      DATA PMDRAGE/28*0/ 
      DATA PMDZML/1988*0./ 
C 
      END 
       
