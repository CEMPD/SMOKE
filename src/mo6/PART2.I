C 
C  Changes: (Last change first) 
C 
C 15 Feb 01: AIR Task 11: New INCLUDE. 
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  /PART2/ 
C  BSWGTR   R    Base gasoline fuel sulfur content, in wt% 
C  CTLFRC   R    Catalyst fraction for a given IV and MY  
C  FECONM   R    Fuel economy, mpg 
C  FRACA    R    Fraction of vehicles with: 
C                 1: Ox No Air 
C                 2: Ox Air 
C                 3: 3w No Air 
C                 4: 3w Air 
C  GSWPTR   R    In-use gasoline fuel sulfur content, in wt% 
C  RMIS     R    Misfueling rate 
C  SULFCF   R    Fuel Sulfur Correction Factor 
C  TAMFRC   R    PM Tampering fraction 
C 
C The data is initialized in BD43. 
C 
C 
      REAL    BSWGTR 
      REAL    CTLFRC 
      REAL    FECONM 
      REAL    FRACA 
      REAL    GSWGTR 
      REAL    RMIS 
      REAL    SULFCF 
      REAL    TAMFRC 
       
C 
      COMMON /PART2/ RMIS,FRACA(4),CTLFRC,FECONM, 
     *               GSWGTR,BSWGTR,SULFCF,TAMFRC 
       
