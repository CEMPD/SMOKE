C 
C  Changes: (Last change first) 
C 
C 15 Nov 96 DynTel@MLA Expanded array CCSTD to accomodate the new LDGT 
C           vehicle categories. 
C 
C 
C The data is initialized in BD17. 
C 
      INTEGER NCC, MYCOL 
      REAL CCSTD 
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(5,2)                                                        
