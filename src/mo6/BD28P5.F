      BLOCK DATA BD28P5 
C 
C  BLOCK DATA Subprogram 28: /RLCOM4/ 
C 
C  BLOCK DATA 28 defines the road fuel economy rates for MC only. 
C 
C Changes: 
C 
C 13 Feb 81 AIR Task 11: Trimmed to only initialize MC. All other 
C           gasoline vehicle fuel economies are initialized in BD28. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block array subscripts: 
C  ROADFEMC(12)  -  ROADFEMC ( MDX ) 
C 
C  Common block dictionary: 
C 
C   Name   Type             Description 
C  ------  ----  -------------------------------------- 
C  /RLCOM4/ block data supplied road fuel economy rates 
C 
C  ROADFEMC   R    road fuel economy rates for motorcycles 
C 
C 
C  Notes: 
C 
C  The MDX Mapping to MY for ROADFEMC (MC only) is: 
C 
C   MDX     MY 
C    1      PRE -1977 
C    2-11   1977-1986 (i.e., MDX=MY-1976+1) 
C    12     POST 1987 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'RLCOM4.I' 
C 
      DATA ROADFEMC/ 
C  MC 
     5 65.937, 65.937, 64.120, 63.137, 61.070, 57.834, 
     5 57.394, 57.583, 56.782, 57.102, 56.794, 56.794/ 
C 
      END 
