      SUBROUTINE GETVORR(INPREC,INERR,RC)
C
C  GETVORR processes the OBD EVAP RESPONSE (RATE) labeled input record.
C
C  Called by DSPTCH.
C
C  Calls NXTTOK, QUITER
C
C  Changes: (Last change first)
C
C 27 Jun 02 @EPA-BG QUITER message numbers
C 04 Apr 02 AIR Task 32: OBD Evap Response modifications (New M6 Routine).
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
C    /LDGOBD/ OBDVA,OBDVF,OBDVFLG,OBDVLIN
C    /MAXIMA/ MAXYRS
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
      INCLUDE 'MAXIMA.I'
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
      OBDVFLG=.FALSE.
C
C     Look for the first data field following the delimiter in column 20.
C     Initialize by setting PTR2 to 21 (the position of the blank following
C     the delimiter). PTR1 does not need to be initialized.
C
      PTR2 = 21
C
C     Obtain the 2 evap inflection points, in years
C
      DO I=1,2
C
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
        READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) OBDVA(I)
C
        IF(OBDVA(I).LT.1 .OR. OBDVA(I).GT.MAXYRS) THEN   ! 1 to 25 years
          CALL QUITER(FLOAT(OBDVA(I)),3,588,INERR)
          RC = 2
        ENDIF
C
      END DO
C
C     Ensure that 2nd inflection point occurs after the 1st
C
      IF(OBDVA(2).LE.OBDVA(1)) THEN
        CALL QUITER(FLOAT(OBDVA(2)),3,587,INERR)
        RC = 2
      END IF
C
C     Obtain the 3 evap response rates
      
      DO I=1,3
C
C       Read the data.
C
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
        READ (DATAFIELD,*,ERR=95,IOSTAT=IOS) OBDVF(I)
C
        IF(OBDVF(I).LT.0.0 .OR. OBDVF(I).GT.100.) THEN
          CALL QUITER(OBDVF(I),3,522,INERR)
          RC = 2
        ENDIF
C
C       Convert the input from a percent to a fraction.
C
        OBDVF(I) = OBDVF(I) / 100.0
C
      END DO
C
C     Allow OBD EVAP RESPONSE only if all inputs are correct
C
      IF(RC.EQ.0.AND.INERR.EQ.0) OBDVFLG=.TRUE.
C
C     Obtain the Linear switch keyword, if available (default is STEP)
C
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2)
      OBDVLIN=.FALSE.  ! Step
      IF(TRIM(DATAFIELD).EQ.'LINEAR') OBDVLIN=.TRUE. ! Linear
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
