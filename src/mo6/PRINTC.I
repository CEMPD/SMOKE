C 
C  Common blocks: 
C 
C     /PRINTC/ - Holds printer control characters and other output 
C                formatting data. 
C 
C  The data is initialized in BD24. It may be reset when a labeled input 
C  record designating a specific platform is found in the input. 
C 
C  Changes: (Last change first) 
C 
C 
C 
C This include file is new for Mobile6. 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  /PRINTC/: 
C 
C  NEWPAGE   C   Form feed (page break) character. 
C 
C  The data is initialized in BD24. 
C 
      CHARACTER  NEWPAGE 
C 
      COMMON /PRINTC/ NEWPAGE 
 
