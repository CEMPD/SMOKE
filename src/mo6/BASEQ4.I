C 
C  Changes: (Last change first) 
C 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; Expanded array 
C                sizes further to 8400 
C     21 July 99@ Dyntel JWR Expanded array sizes further to 7920 
C                (=33*3*16*2  + 33*3*12*2*2) to allow all user 
C                BERs in a single run 
C     9 July 99  @Dyntel JWR Expanded array sizes from 100 to 250 
C                to allow more user input BERs for each run. 
C                Arrays affected are NEWPAR, BERNEW, NEWFIT, NEWOML. 
C  The data is initialized in BD03. 
C 
      INTEGER NEWPAR,NEWFIT,NEWCT,MAXCT 
      REAL BERNEW 
      CHARACTER *1 NEWOML 
      COMMON /BASEQ4/NEWPAR(5,8400),BERNEW(3,8400),NEWFIT(8400),NEWCT 
     *              ,MAXCT,NEWOML(8400) 
