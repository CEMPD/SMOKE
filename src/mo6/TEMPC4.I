C 
C  Changes: (Last change first) 
C 
C 12 Dec 96 DynTel@MLA Expanded array MYGCOO to accomodate the new LDGT 
C           vehicle categories. 
C 
C 
C The data is initialized in BD17 
C 
      INTEGER MYGCOO, LOWCO 
      REAL COOTFC, COOTFD 
      COMMON /TEMPC4/ MYGCOO(5),LOWCO,COOTFC(3),COOTFD(3)                                            
