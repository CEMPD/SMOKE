C       
C  Changes: (Last change first) 
C 
C   7 Nov 96 DynTel@MLA  Changed the dimension of LB1STS from 3 to 5 to 
C            accomodate the new LGDT vehicle types (LDGT1A and LDGT1B 
C            replace LDGT1, and LDGT2A and LDGT2B replace LDGT2). 
C 
C 
C   The data is initialized in BD02. 
C 
      INTEGER MAXSP1, LB1STS, LBLAST, IVB, IVA 
      COMMON /SPEED4/ MAXSP1,LB1STS(5),LBLAST,IVB,IVA                                                
