C 
C  YEARS4 sets four different calendar years 
C 
C  Changes: (Last change first) 
C 
C  Parameter Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ------------------------------------------------- 
C  IYSTART   I    4-digit representation of start year 
C  IY1960    I    4-digit representation of year 1960 
C  IY1994    I    4-digit representation of year 1994 
C  IYEND     I    4-digit representation of end year  
C 
C  Notes: 
C 
C 
C  The data is initilized in BD16. 
C  
      INTEGER IYSTART, IY1960, IY1994, IYEND 
      COMMON /YEARS4/ IYSTART,IY1960,IY1994,IYEND                                                    
