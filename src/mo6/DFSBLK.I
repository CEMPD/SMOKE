C  DFSBLK.I 
C 
C Changes: 
C 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /DFSBLK/: 
C 
C  DFS      R    Fraction of vehicles designed to run on leaded 
C                gasoline that use unleaded gasoline, i.e. 
C                discretionary fuel switching rates, 1-DFS is the 
C                fraction of vehicles designed to run on leaded 
C                gasoline that use leaded gasoline 
C  MYLEAD   I    Model year groups for discretionary fuel switching 
C  ABURN    R    Fraction of lead burned that is exhausted 
C                (denoted as a in particulate document) 
C  MAXLED   I    Number of model year groups for base lead factors 
C 
C 
C  The data are initialized in BD01P5. 
C 
C 
      INTEGER MAXLED 
      INTEGER MYLEAD 
      REAL    ABURN 
      REAL    DFS 
C  
      COMMON /DFSBLK/ DFS(4,5), MYLEAD(4,5), ABURN(3), MAXLED 
C 
