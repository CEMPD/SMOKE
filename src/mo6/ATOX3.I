C 
C  Changes: (Last change first) 
C 
C 30 Jun 01: AIR Task 11: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      REAL GROUP 
C 
      REAL BENZRATIO 
      REAL BUTARATIO 
      REAL ACETRATIO 
      REAL FORMRATIO 
      REAL MTBERATIO 
      REAL VOCCOMPLEX 
C 
      COMMON /ATOX3/ BENZRATIO(10), BUTARATIO(10), ACETRATIO(10),  
     &               FORMRATIO(10), MTBERATIO(10), GROUP(10), 
     &               VOCCOMPLEX(10) 
