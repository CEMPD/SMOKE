      SUBROUTINE GETORE(INPREC,INERR,RC)
C
C  GETORE processes the OBD REPAIRED EFS labeled input record.
C
C  Called by DSPTCH.
C
C  Calls NXTTOK, QUITER
C
C  Changes: (Last change first)
C
C 02 May 02 AIR Task 32: Changed QUITER call due to addition
C           of more OBDTYPEs in QUITER
C 27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads
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
C    /LDGOBD/ OBDREP
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
      READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) OBDREP
C
C     Check the data.
C
      IF(OBDREP.LT.0.0.OR.OBDREP.GT.100.) THEN
        CALL QUITER(OBDREP,5,522,INERR)
        RC = 2
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
