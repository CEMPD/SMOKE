C 
C  Changes: (Last change first) 
C 
C 21 Apr 98 @DynTel-RJD  Added an extra dimension to SALHCF to distinguish 
C            between start and running emissions. 
C  2 Apr 97 - DynTel@MLA  Expanded array SALHCF to accomodate the new BUS  
C             vehicle categories. 
C   5 Nov 96 DynTel@MLA  Changed the last dimension of SALHCF from 8 to 24 
C            to acount for 7 new HDG and 7 new HDD vehicle types and the 
C            new LGDT1A and LDGT1B types which replace LDGT1, and the new 
C            LDGT2A and LDGT2B types which replace LDGT2. 
C 
C 
C 
C  The data is initialized in BD02. 
C 
      REAL SALHCF, SCFIDL 
      LOGICAL IDLYES 
      COMMON /SPEED6/ SALHCF(25,3,28,2),SCFIDL,IDLYES                                                   
