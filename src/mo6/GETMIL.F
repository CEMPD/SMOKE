      SUBROUTINE GETMIL(INPREC,INERR,RC) 
C 
C  GETMIL processes the OBD WORKING MIL labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK, QUITER 
C 
C  Changes: (Last change first) 
C 
C 27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C 15 Sep 00 AIR Task 08: Change MIL to WMIL for easier variable tracking 
C 12 May 99 @DynTel-MLA 2-694  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /LDGOBD/ WMIL 
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
      INCLUDE 'LDGOBD.I' 
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
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) WMIL 
C 
C     Check the data. 
C 
      IF(WMIL.LT.0.0.OR.WMIL.GT.100.) THEN 
        CALL QUITER(WMIL,1,522,INERR) 
        RC = 2 
      ENDIF 
C 
C     Convert the input from a percent to a fraction for use by 
C     the I/M OBD code. 
C 
      WMIL = WMIL / 100.0 
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
