      BLOCK DATA BD46 
C 
C  Ammonia (NH3) Basic Emission Rates (g/mi) 
C  and Sulfur Deterioration Rates (g/mi/ppm). 
C 
C  Changes : 
C  28 Mar 01 AIR Task 11: Due to Lahey bug, had to replace NH3BER 
C            MODULE with GETNHEBER, GETNH3SDR, PMNH3BER.D and PMNH3SDR.D 
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
      IMPLICIT NONE 
C 
      INCLUDE 'PART6.I' 
C 
      DATA NH3BER/336*0./ 
      DATA NH3SDR/336*0./ 
C 
      END 
