      SUBROUTINE GETBARO(INPREC,INERR,RC)
C
C  GETBARO processes the BAROMETRIC PRES labeled input record.
C
C  Called by DSPTCH.
C
C  Calls NXTTOK,QUITER
C
C  Changes: (Last change first)
C
C 11 Oct 02 @epa-bag bug 36: Extended values of limits on barometeric pressure
C 02 Jul 02 AIR Task 36: New Routine
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
C    /ALUHIN/ BARO, BARO_FLAG
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
      INCLUDE 'ALUHIN.I'
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
      BARO_FLAG=.FALSE.
C
C     Look for the first data field following the delimiter in column 20.
C     Initialize by setting PTR2 to 21 (the position of the blank following
C     the delimiter). PTR1 does not need to be initialized.
C
      PTR2 = 21
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
C
C     Read and check the barometric pressure. Error out if outside limits.
C
      READ(DATAFIELD,*,ERR=95,END=95,IOSTAT=IOS) BARO
C
      IF (BARO.LT.13.00 .OR. BARO.GT.33.00) THEN
        CALL QUITER(BARO,0,4,INERR)
        RC = 1
        GOTO 99
      END IF
      
      BARO_FLAG=.TRUE.
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
