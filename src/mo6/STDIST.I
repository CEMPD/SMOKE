C 
C  Changes: (Last change first) 
C 
C 
C 10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C 23 May 00 @EPA-djb Added SDIST. 
C  6 May 99 @DynTel-MLA 2-650  Changed the number of hourly intervals in 
C           STARTDIST from 14 to 24. 
C         
C 
C The data is initialized in BD42. 
C 
      REAL STARTDIST,STARTSPERDAY 
      INTEGER SDIST 
      COMMON /STDIST/  STARTDIST(24,2),STARTSPERDAY(25,28,2),SDIST 
