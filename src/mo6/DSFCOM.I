C 
C  Changes: (Last change first) 
C 
C  25 Jan 2000 @EPA-djb Task 27 Replaced arrays with final values and 
C             restructured the array sizes. 
C  26 Nov 96 DynTel@MLA Expanded the arrays to accomodate the new Heavy Duty 
C           Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD23. 
C 
      INTEGER MAXDSF 
      REAL DSFRAC 
      COMMON /GSFCOM/ MAXDSF,DSFRAC(25,14) 
