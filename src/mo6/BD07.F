      BLOCK DATA BD07 
C 
C  BLOCK DATA Subprogram 7: /OPCNTL/ 
C 
C  BLOCK DATA 7 defines the variables that determine the output formats. 
C 
C  Changes : 
C 
C  8 Jun 00 @EPA-ddj,Bug120, Removed SSFLAG. 
C 
C  Common block array subscripts: 
C 
C  OPTABL(4)       -  OPTABL ( IVEXP ) 
C 
C  Common block dictionary: 
C 
C   Name   Type               Description 
C  ------  ----  -------------------------------------------- 
C  /OPCNTL/: 
C 
C  OPTABL   I    Controls output of the expanded tables for some vehicles. 
C                  1 = skip    2 = print   3 = print with I/M table 
C                   OPTABL(1) = LDGT 
C                   OPTABL(2) = HDGV 
C                   OPTABL(3) = HDDV 
C                   OPTABL(4) = BUS and CNG 
C 
      IMPLICIT NONE 
      INCLUDE 'OPCNTL.I' 
C 
C  /OPCNTL/: 
C 
      DATA OPTABL/4*1/ 
C 
      END 
