C 
C  Changes: (Last change first) 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased 
C         size of MYGERU 
C   9 July 99 @Dyntel- JWR Increased the first dimension of MYGERU to 33 
C            to allow user input BERs from 1950 - present. 
C 
C   5 Nov 96 DynTel@MLA  Changed the fourth dimension of MYGERU and the  
C            second dimension of NUMERU from 8 to 24 to accomodate 7 new 
C            HDG, the 7 new HDD vehicle types, and the 2 new LDGT types. 
C 
C  2 Apr 97 - DynTel@MLA  Expanded arrays MYGERU and NUMERU to accomodate  
C             the new BUS vehicle categories. 
C 
C  12 May 98 @DynTel-RJD Added an extra dimension to MYGERU and NUMERU to 
C             distinguish between start and running emissions. 
C 
C  The initialization of the data could be found in BD03. 
C 
      INTEGER MYGERU,MAXERU,NUMERU 
      COMMON /BASEQ6/ MYGERU(35,2,3,28,2,2),MAXERU,NUMERU(3,28,2,2) 
