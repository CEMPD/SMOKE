C 
C  Changes: (Last change first) 
C 
C  2 Apr 97 - DynTel@MLA  Expanded array GPBSCO to accomodate the new BUS  
C             vehicle categories. 
C  5 Nov 96 - DynTel@MLA  Changed last dimension of GPBSCO from 4 to 18 
C             to account for 7 new HDG and 7 new HDD vehicle categories. 
C 
C 
C 
C 
C  The data is initialized in BD02. 
C 
      REAL GPBSCO 
      COMMON /SPEED7/ GPBSCO(3,3,22) 
