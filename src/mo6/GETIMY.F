      SUBROUTINE GETIMY(INPREC,INERR,RC) 
C 
C  GETIMY processes the I/M MODEL YEARS labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C  04 Jan 00 BAG Bug294 Changed i&m model year out of bounds errors to 
C            QUITER messages 12 and 13 
C  27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C  25 Aug 00 AIR Task 08: Allowed for up to 10 I/M and EVAP programs. 
C  16 Jul 99 @DynTel-bag 2-705  Limits on i&m model years changed to 1941 and 2050 
C 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /YEARS4/ IYEND,IYSTART 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /CHKLST/ IMCHECK 
C    /IMPAR1/ MODYR 
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
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'YEARS4.I' 
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
      READ(DATAFIELD,*,ERR=95,END=95) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program and check off 
C     the I/M MODEL YEARS data item for this I/M program. 
C 
      IMCHECK(1,IMPGM) = 2 
      IMCHECK(8,IMPGM) = 2 
C 
C     Read the first and last model years covered. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      MODYR(1,IMPGM)=INT(X) 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      MODYR(2,IMPGM)=INT(X) 
C 
C     Check the data. 
C 
      IF (MODYR(1,IMPGM).LT.1941.OR. 
     &    MODYR(1,IMPGM).GT.2050) THEN 
        CALL QUITER(0.0,0,12,INERR) 
        RC = 2 
      ENDIF 
C 
      IF (MODYR(2,IMPGM).LT.1941.OR. 
     &    MODYR(2,IMPGM).GT.2050) THEN 
        CALL QUITER(0.0,0,13,INERR) 
        RC = 3 
      ENDIF 
C 
      IF(MODYR(1,IMPGM).GT.MODYR(2,IMPGM)) THEN 
        CALL QUITER(0.,0,37,INERR) 
        RC = 4 
      ENDIF 
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
