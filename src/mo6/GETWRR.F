      SUBROUTINE GETWRR(INPREC,INERR,RC) 
C 
C  GETWRR processes the I/M WAIVER REDUCTION labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C 27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C 28 Aug 00 AIR Task 08: Expanded to handle EVAP I/M. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /CHKLST/ IMCHECK 
C    /LDGIMC/ WAVRDC 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  IMPGM      I   The I/M program number which uses this data. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  X          R   Temporary real to integer variable 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'LDGIMC.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        IMPGM 
      INTEGER        IOS 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
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
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program. 
C 
      IMCHECK(1,IMPGM) = 2 
      IMCHECK(27,IMPGM) = 2 
C 
C     Read the waiver reduction rate for the I/M program. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) WAVRDC(IMPGM) 
C 
C     Check the data. 
C 
      IF(WAVRDC(IMPGM).LT.0.0.OR.WAVRDC(IMPGM).GT.100.0) THEN 
        CALL QUITER(WAVRDC(IMPGM),0,34,INERR) 
        RC = 2 
      ENDIF 
C 
C     Scale the waiver reduction rate as expected by the code. 
C 
      WAVRDC(IMPGM) = WAVRDC(IMPGM) / 100.0 
C 
      GOTO 99 
C 
C     Error or EOF on any attempted read. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
C 
   99 RETURN 
      END 
