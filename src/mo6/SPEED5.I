C  Changes: (Last change first) 
C 
C   7 Nov 96 DynTel@MLA  Changed the dimension of MAXHSP from 3 to 5 and 
C            the last dimension of MYHSPC from 3 to 5 to accomodate the 
C            new vehicle types LDGT1A and LDGT1B which replace LDGT1, and 
C            the new types LDGT2A and LDGT2B which replace LDGT2. 
C  
C       
C   The data is initialized in BD02. 
C      
      INTEGER MYHSPC, MAXHSP 
      COMMON /SPEED5/ MYHSPC(2,5),MAXHSP(5) 
