      SUBROUTINE GETD_PPM(INPREC,INERR,RC) 
C 
C  GETD_PPM processes the DIESEL SULFUR labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C 30 Jan 01 AIR Task 11: New routine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    NONE 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /PART1/ D_PPM,PRTCHK 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'PART1.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        IOS 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
C 
C     Read and check the diesel sulfur data. 
C 
      READ(DATAFIELD,*,ERR=95,END=95,IOSTAT=IOS) D_PPM 
C 
      IF (D_PPM.LT.0.01 .OR. D_PPM.GT.5000.) THEN 
        CALL QUITER(D_PPM,0,181,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that the required diesel sulfur was entered for this particulate program. 
C 
      PRTCHK(2) = 2 
C 
      GOTO 99 
C 
C Error or EOF on any attempted read. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
C 
   99 RETURN 
      END 
