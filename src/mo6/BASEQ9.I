C 
C  Changes: (Last change first) 
C 
C  22 Jun 99 @DynTel-MLA 2-694  Removed P94REF and P94SEF. Moved VA015 
C            to SAVE07. 
C 
C     The data is initialized in BD05 
C 
      INTEGER MAXLYR 
      PARAMETER (MAXLYR = 32) 
      REAL LEVIMP, LDTAB 
      COMMON /BASEQ9/ LEVIMP(11,MAXLYR,5), LDTAB(2,2) 
