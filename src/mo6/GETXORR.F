      SUBROUTINE GETXORR(INPREC,INERR,RC)
C
C  GETXORR processes the OBD EXHAUST RESPONSE (RATE) labeled input record.
C
C  Called by DSPTCH.
C
C  Calls NXTTOK, QUITER
C
C  Changes: (Last change first)
C
C 23 Sep 02 AIR Task 35b: prevent the use of OBD EXH RESPONSE whenever an
C           I/M program is present.
C 27 Jun 02 @EPA-BG QUITER message numbers
C 04 Apr 02 AIR Task 32: OBD Exhaust Response modifications (New M6 Routine).
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
C    /LDGOBD/ OBDXM,OBDXF,OBDXFLG,OBDXLIN
C
C  Local variable / array dictionary:
C
C   Name    Type                      Description
C  ------   ----  ---------------------------------------------------------
C  DATAFIELD  C   One blank-delimited data item from the input record.
C  I          I   Loop control variable.
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
      INTEGER        I
      INTEGER        IOS
      INTEGER        PTR1
      INTEGER        PTR2
C
C .....................................................................
C
      RC = 0
      OBDXFLG=.FALSE.
C
C     Do not allow the OBD EXH RESPONSE command to work whenever an
C     I/M program is input.
C
      IF(OBDIMFLG) THEN
        CALL QUITER(0.0,0,118,INERR)
        GOTO 99
      END IF
C
C     Look for the first data field following the delimiter in column 20.
C     Initialize by setting PTR2 to 21 (the position of the blank following
C     the delimiter). PTR1 does not need to be initialized.
C
      PTR2 = 21
C
C     Obtain the 2 exhaust inflection points, in 1000s of miles
C
      DO I=1,2
C
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
        READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) OBDXM(I)
C
        IF(OBDXM(I).LT.0.0 .OR. OBDXM(I).GT.200.) THEN   ! 0 to 200,000 miles
          CALL QUITER(OBDXM(I),2,588,INERR)
          RC = 2
        ENDIF
C
C       Convert the input from 1000s of miles to 10,000s of miles.
C
        OBDXM(I) = OBDXM(I) / 10.0
C
      END DO
C
C     Ensure that 2nd inflection point occurs after the 1st
C
      IF(OBDXM(2).LE.OBDXM(1)) THEN
        CALL QUITER(OBDXM(2),2,587,INERR)
        RC = 2
      END IF
C
C     Obtain the 3 exhaust response rates
      
      DO I=1,3
C
C       Read the data.
C
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
        READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) OBDXF(I)
C
        IF(OBDXF(I).LT.0.0 .OR. OBDXF(I).GT.100.) THEN
          CALL QUITER(OBDXF(I),2,522,INERR)
          RC = 2
        ENDIF
C
C       Convert the input from a percent to a fraction.
C
        OBDXF(I) = OBDXF(I) / 100.0
C
      END DO
C
C     Allow OBD EXH RESPONSE only if all inputs are correct
C
      IF(RC.EQ.0.AND.INERR.EQ.0) OBDXFLG=.TRUE.
C
C     Obtain the Linear switch keyword, if available (default is STEP)
C
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
      OBDXLIN=.FALSE.  ! Step
      IF(TRIM(DATAFIELD).EQ.'LINEAR') OBDXLIN=.TRUE. ! Linear
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
