C 
C  Changes: (Last change first) 
C 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used 
C            for user input of ETP phase in 
C 
C  The data is initialized in BD11. 
C 
      INTEGER IETPYR,NETPYR 
      REAL ETP,ETPPIV 
      COMMON /ETPPAR/ ETP(7,3),IETPYR(4),ETPPIV(4),NETPYR 
