C 
C  Changes: (Last change first) 
C 
C  8 Apr 99 @DynTel-DDJ 2-696 Expanded array ALHRET to accomodate hourly intervals 
C  3 Dec 96 DynTel@MLA Expanded array ALHRET to accomodate the new LDGT 
C           vehicle categories. 
C 
C  The initialization of the data can be found in BD27. 
C 
      REAL ALHRET 
      COMMON /ALUOUT/ ALHRET(25,6,3,24) 
