C 
C  Changes: (Last change first) 
C 
C   7 Nov 96 DynTel@MLA  Changed the last dimensions of MYGSP1 and MYGSP2 
C            from 3 to 5 to accomodate new bvehicle types LDGT1A and 
C            LDGT1B (in place of LDGT1) and LDGT2A and LDGT2B (in place of 
C            LDGT2). 
C 
C 
C   The data is initialized in BD02. 
C 
      INTEGER MYGSP1, MYGSP2, MYGSMC 
      COMMON /SPEED2/ MYGSP1(18,2,5),MYGSMC(18,2,3),MYGSP2(3,5)                                      
