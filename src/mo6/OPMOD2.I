C 
C  Changes: (Last change first) 
C 
C 20 Nov 96 DynTel@MLA Expanded array MYGBFR to accomodate the new LDGT 
C           vehicle categories. 
C 
C 
C The data is initialized in BD19. 
C 
      INTEGER MYSTEP, MAXBFR, MYGBFR 
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,8) 
