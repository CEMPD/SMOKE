C  IDLE.I 
C 
C Changes: 
C 
C 14 Nov 00 AIR Task 06: Added NIDLYR to common BSIDL2. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /BSIDL/: 
C 
C  BASIDL   R    base idle emission rates 
C 
C  /BSIDL2/: 
C 
C  IDLYRS   I    base idle emission years 
C  NIDLYR   I    number of idle years 
C 
C  /IDLCOM/: 
C 
C  IDLEPMFLG L   flag to print idle PM factors 
C 
C 
C  The data are initialized in BD08P5. 
C 
C 
      INTEGER IDLYRS 
      INTEGER NIDLYR  
      LOGICAL IDLEPMFLG 
      REAL    BASIDL 
C  
      COMMON /BSIDL/  BASIDL(4,5) 
      COMMON /BSIDL2/ IDLYRS(4), NIDLYR 
      COMMON /IDLCOM/ IDLEPMFLG 
C 
