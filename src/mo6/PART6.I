C 
C  PART6 contains the Ammonia (NH3) Basic Emission Rates 
C  and Sulfur Deterioration Rates 
C 
C  Changes: (Last change first) 
C 
C  09 Apr 01: AIR Task 11: New INCLUDE. 
C 
C  NH3BER(3,4,28) - NH3BER(ICAT,IROAD,IVEH) 
C  NH3SDR(3,4,28) - NH3SDR(ICAT,IROAD,IVEH) 
C 
C  ICAT is for catalyst, as follows: 
C    1 = no catalyst or diesel 
C    2 = oxidation catalyst 
C    3 = 3-way catalyst 
C 
C  IROAD and R are the roadway type, as follows: 
C    1 = freeways 
C    2 = arterial 
C    3 = local 
C    4 = ramp 
C 
C  NH3BER values are given in g/mi. 
C  NHRSDR values are given in g/mi/ppm. 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C /PART6/ 
C   NH3BER   R    Ammonia (NH3) Basic Emission Rate (g/mi) 
C   NH3SDR   R    Ammonia (NH3) Sulfur Deterioration Rate (g/mi/ppm) 
C 
C 
C The data is initialized in BD46. 
C 
      REAL NH3BER 
      REAL NH3SDR 
C 
      COMMON /PART6/ NH3BER(3,4,28),NH3SDR(3,4,28) 
