C 
C  VDATA contains the variable that holds the value of the name of the version 
C  of the model. 
C 
C  Changes: (Last change first) 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  VERSION   C    holds the value of the name of the version of the model. 
C 
C 
C  Notes: 
C 
C  12-Dec-1994 @ CSC-rwk Req. 2-000 VDATA was created to make it easier to 
C  change the name of the version of the model.  This was done for quality 
C  control. 
C  15-Oct-1997 @ DynTel-gjr 2-000 Changed VERSON to VERSION.  
C 
C  The data is initialized in VNAME.I 
C  
      CHARACTER*80 VERSION 
      COMMON /VDATA/ VERSION 
