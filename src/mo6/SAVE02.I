C 
C  Changes: (Last change first) 
C 
C  15 Jan 97  Changed the last dimension of arrays VA003 and VA004 from       
C       16 to 18 to match the declarations of CRD4VA and CRD4VB in include 
C       files ICR4VA and ICR4VB, even though variables VA003 and VA004 are  
C       not currently used in Mobile6. 
C 
C  The data is initiliazed in BD12. 
C 
      REAL VA003,VA004 
      COMMON /SAVE02/ VA003(24,3,18),VA004(24,3,18) 
