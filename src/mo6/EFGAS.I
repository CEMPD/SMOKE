C  EFGAS.I 
C 
C 
C  EFGAS contains the  EPA-based Zero Mile Levels (ZML) 
C  and Deterioration Rates (DR) for gasoline vehicles 
C 
C Changes: 
C 
C 16 Oct 01 AIR Task 25: New MOBILE6.2 include. 
C 18 Apr 01 AIR Task 16: New MOBILE6.1 include. 
C 14 Nov 00 AIR Task 06: Added NCAT to common EFPMC. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  PMGDR1(71,4,28) - PMGDR1(MY,ICAT,IV) 
C  PMGDR2(71,4,28) - PMGDR2(MY,ICAT,IV) 
C  PMGZML(71,4,28) - PMGZLM(MY,ICAT,IV) 
C 
C  MY is model year, 1950-=1, 2020+=71 
C  ICAT is described below 
C 
C 
C  Common Block Dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /EFPMC/: 
C 
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
C  The data are initialized in BD44. 
C 
C 
      INTEGER NCAT 
      REAL    PMGZML 
      REAL    PMGDR1 
      REAL    PMGDR2 
C 
      COMMON /EFPMC/ PMGZML(71,4,28),PMGDR1(71,4,28), 
     *               PMGDR2(71,4,28), NCAT 
