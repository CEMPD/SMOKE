      SUBROUTINE GETE200(INPREC,INERR,RC)
C
C  GETE200 processes the E200 labeled input record.
C
C  Called by DSPTCH.
C
C  Calls NXTTOK,QUITER
C
C  Changes: (Last change first)
C
C 02 Jul 02 AIR Task Bug Fix 27: Changed limits on AT inputs.
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
C    /ATOX1/ E200,TOXCHK
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
      INCLUDE 'ATOX1.I'
C
      CHARACTER*(*)  INPREC
      INTEGER        INERR
      INTEGER        RC
C
      CHARACTER*20   DATAFIELD
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
C     Read and check the E200 data. If out of range, warn the user that
C     the results are most likely erroneous
C
      READ(DATAFIELD,*,ERR=95,END=95,IOSTAT=IOS) E200
C
      IF (E200.LT.30. .OR. E200.GT.70.) CALL QUITER(E200,0,142,INERR)
C
C     Show that the required E200 was entered for this air toxic program.
C
      TOXCHK(1) = 2
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
