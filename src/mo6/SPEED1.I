C 
C  Changes: (Last change first) 
C 
C   7 Nov 96 DynTel@MLA  Changed the last dimensions of SK2C and SK4C 
C            from 3 to 5 to accomodate new vehicle types LDGT1A and 
C            LDGT1B (in place of LDGT1) and LDGT2A and LDGT2B (in place 
C            of LDGT2). 
C 
C 
C   The data is initialized in BD01. 
C 
      REAL SK1, SK2C, SK4C 
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,5),SK4C(2,14,3,5)                                      
