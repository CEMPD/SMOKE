C 
C  Changes: (Last change first) 
C 
C 15 Nov 96 DynTel@MLA Expanded arrays MYGTCF and MYGTFD to accomodate the  
C           new Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C  7 Apr 97 - DynTel@MLA  Expanded array MYGTCF to accomodate the new BUS  
C             vehicle types. 
C 
C 
C  The data is initilized in BD17 
C 
      INTEGER MAXTCF, MYGTCF, MYGTFD, ISHIFT 
      COMMON /TEMPC2/ MAXTCF,MYGTCF(9,3,15),MYGTFD(5),ISHIFT                                          
