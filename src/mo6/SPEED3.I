C 
C  Changes: (Last change first) 
C 
C  2 Apr 97 - DynTel@MLA  Expanded array SCUNA1 to accomodate the new BUS  
C             vehicle categories. 
C   5 Nov 96 DynTel@MLA  Changed the last dimension of SCUNA1 from 8 to 24 
C            to acount for 7 new HDG vehicle types, the 7 new HDD types, 
C            and the two new LDGT types (LDGT1A AND LDGT1B in place of 
C            LDGT1, and LDGT2A and LDGT2B in place of LDGT2). 
C 
C 
C 
C 
C  The data is initialized in BDSPEED1. 
C 
      REAL SADJ, SCUNA1 
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,28)                                                            
