      SUBROUTINE GETHSGLL(INPREC,INERR,RC)
C
C  GETFPPR reads in the hot soak liquid leak rates
C
C  Called by DSPTCH.
C
C  Calls QUITER.
C
C  Changes: (Last change first)
C
C  02 May 02 @EPA-bag user input of hot soak liquid leak rates
C
C  Input on call:
C
C    argument list: INPREC
C    common blocks:
C
C  Output on return:
C
C   argument list: INERR,RC
C    common block:
C    /IOUCOM/ IOUREP
C
C    modules
C    LEAKYDAT : HSUSERGL_INPT, HSGLL_USER
C
C  Local variable / array dictionary:
C
C  Name    Type              Description
C  ------  ----  ----------------------------------------------------
C  SUM      I    Summation
C  I        I    Do Loop Index
C
      USE LEAKYDAT, ONLY : HSUSERGL_INPT, HSGLL_USER
C
      IMPLICIT NONE
      INCLUDE 'IOUCOM.I'
C
      CHARACTER INPREC*(*)
      INTEGER   INERR
      INTEGER   RC
      INTEGER   I
      INTEGER, SAVE :: IYES = 1
C
      RC = 0
C
C  Obtain the 4 Fixed P/P Rates: Pass/Pass, Pass Pressure/Fail Purge,
C  Fail Pressure/Pass Purge, Fail/Fail.
C
      WRITE(IOUREP,10)
  10  FORMAT(/,2x,
     &         'Reading User Supplied Hot Soak Liquid Leak Rates.')
C
      READ(INPREC(21:),*,ERR=95,END=96) HSGLL_USER
C
C  Check that the rates are greater than 0. If they don't, warn the user.
C
C
      DO I=1,3
        IF(HSGLL_USER(I).LT.0) THEN
          CALL QUITER(HSGLL_USER(I),I,524,INERR)
          GOTO 95
        ENDIF
      ENDDO
C
      HSUSERGL_INPT=IYES
      CALL QUITER(0.0,0,700,INERR)
      WRITE(IOUREP,100) HSGLL_USER(1),HSGLL_USER(2),HSGLL_USER(3)
  100 FORMAT(/,2x,
     &         'User Supplied Hot Soak Liquid Leak Rates:'/
     &         10X,'Carbureted      :',1X,F6.2,1X,'gram/hour'/
     &         10X,'Throttle-body FI:',1X,F6.2,1X,'gram/hour'/
     &         10X,'Ported FI       :',1X,F6.2,1X,'gram/hour'/)
      GOTO 99
C
C  Error or EOF on any attempted read.
C
   95 RC = 1
      GOTO 99
   96 RC = -1
C
   98 CALL QUITER(0.,0,525,INERR)
C
   99 RETURN
      END
