C  
C  Changes: (Last change first) 
C 
C   7 Nov 96 DynTel@MLA  Changed the last dimension of SIXT5 from 3 to 5 
C            to account for the new LDGT categories (LDGT1A and LDGT1B 
C            replace LDGT1, and LDGT2A and LDGT2B replace LDGT2). 
C 
C 
C   The data is initialized in BD02. 
C 
      REAL SIXT5 
      COMMON /SPEED8/ SIXT5(11,3,5) 
