      BLOCK DATA BD44 
C 
C  EPA-based Zero Mile Levels (ZML) and Deterioration Rates (DR1 and DR2). 
C 
C  Changes: (Last change first) 
C 
C  28 Aug 01: Air Task 23c: Implemented DR1 and DR2 rates and PMDRAGE. 
C  24 Aug 01: AIR Task 23b: EPA-based ZMLs/DRs 
C  05 Jun 01: AIR Task 17: Added CAAROMFLG. 
C  05 Jun 01: AIR Task 17: Added base sulfur level array. 
C  23 May 01: AIR Task 17: Changed from 1975 to 1950 base and only 2 modes. 
C  20 Mar 01: AIR Task 11: New INCLUDE. 
C 
C  BSULFLVL(3)    - BSULFLVL (see below) 
C  PMDR1(71,2,28) - PMDR1(MY,ICAT,IV) 
C  PMDR2(71,2,28) - PMDR2(MY,ICAT,IV) 
C  PMDRAGE(28)    - PMDRAGE(IV) 
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
      IMPLICIT NONE 
C 
      INCLUDE 'PART4.I' 
C 
      DATA BSULFLVL/3*0./ 
      DATA PMDR1/3976*0./ 
      DATA PMDR2/3976*0./ 
      DATA PMDRAGE/28*0/ 
      DATA PMZML/3976*0./ 
C 
      END 
       
