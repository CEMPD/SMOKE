      SUBROUTINE GETTEF(INPREC,INERR,RC) 
C 
C  GETTEF processes the I/M TTC EFFICIENCY labeled input record. 
C                                       
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C 27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C 15 Sep 00 AIR Task 08: Expanded number of I/M programs to handle EVAP 
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
C    /TTCINF/ TTCEFF 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  IMPGM      I   The I/M program number which uses this data. 
C  IOS        I   IOSTAT return from FORTRAN I/O routines. 
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
      INCLUDE 'TTCINF.I' 
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
C 
C     Read the effectiveness rate. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) TTCEFF(IMPGM) 
C 
C     Check the data. 
C 
      IF(TTCEFF(IMPGM).LT.0.0.OR.TTCEFF(IMPGM).GT.100.0) THEN 
        CALL QUITER(TTCEFF(IMPGM),0,511,INERR) 
        RC = 2 
      ENDIF 
C 
      TTCEFF(IMPGM) = TTCEFF(IMPGM) / 100.0 
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
