C 
C  The data is initialized in BD30. 
C 
C  Changes: (Last change first) 
C 
C 
C  Name    Type  Description 
C  ------  ----  ---------------------------------------------------- 
C  RVPBAS   R    SCNAME's base RVP 
C  RVPICY   R    RVP to be used for not cert controlled my for SCNAME: 
C  RVPUWX   R    un-weathered RVPICY 
C 
C  Notes: 
C  
C  6 May 98 @DynTel-MLA 2-659  Added GBLRVP for Mobile6.  
C 
C 
      REAL RVPBAS,RVPICY,RVPUWX,GLRVPBAS 
      COMMON /CITRV1/ RVPBAS,RVPICY,RVPUWX,GLRVPBAS 
