C 
C  Changes: (Last change first) 
C 
C 24 May 01: AIR Task 17: Expanded PRTCHK to 3 required inputs. 
C 09 May 01: AIR Task 11/5: Updated to handle idle. 
C 30 Jan 01: AIR Task 11: New INCLUDE. 
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  /PART1/ 
C  D_PPM     R    diesel sulfur content, in ppm (required input) 
C  FUEL_S    R    EMFAC2000v2.02 fuel sulfur levels 
C  IVEHWL    I    Average Number of Vehicle Wheels for 28 Vehicle Classes 
C  MAXIPPM   I    Maximum number of particulate pollutants, as follows: 
C                 1=SO4, 2=OCARBON, 3=ECARBON, 4=GASPM, 5=Lead, 6=SO2, 7=MH3, 
C                 8=Brake, 9=Tire, 10=EMFAC_PM, 11=Idle 
C  PARTFLG   L    switch to run the particulate subroutine 
C  PRTCHK    I    checklist of required particulate run inputs (1=not input, 2=input) 
C                 (1=particle size, 2=diesel sulfur level, 3=ZML/DR files) 
C  PSZCUT    R    particle size cutoff, in micrometers (required input) 
C  BERSULF   R    Fuel sulfur level associated with the basic emission factor (BER). 
C 
C The data is initialized in BD43. 
C 
      INTEGER IVEHWL 
      INTEGER MAXIPPM 
      INTEGER PRTCHK 
      LOGICAL PARTFLG 
      REAL    BERSULF 
      REAL    D_PPM 
      REAL    FUEL_S 
      REAL    PSZCUT 
C 
      COMMON /PART1/ PRTCHK(3),PSZCUT,D_PPM,MAXIPPM, 
     *               IVEHWL(28),FUEL_S(3,29),PARTFLG, 
     *               BERSULF(25,28) 
