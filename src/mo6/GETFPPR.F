      SUBROUTINE GETFPPR(INPREC,INERR,RC) 
C 
C  GETFPPR reads in the user-supplied Fixed P/P Rates. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  18 Jul 2000 AIR Task 06 New Routine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /IOUCOM/ IOUGEN 
C 
C  Output on return: 
C 
C   argument list: INERR,RC 
C    common block: 
C    /EVAPACT/ FPPRUSER, FPPRFLAG 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ---------------------------------------------------- 
C  SUM      I    Summation 
C  I        I    Do Loop Index 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'EVAPACT.I' 
C 
      CHARACTER INPREC*(*) 
      INTEGER   INERR 
      INTEGER   RC 
      INTEGER   I 
      REAL      SUM 
C      
      RC = 0 
C 
C  Obtain the 4 Fixed P/P Rates: Pass/Pass, Pass Pressure/Fail Purge, 
C  Fail Pressure/Pass Purge, Fail/Fail. 
C 
      WRITE(IOUREP,10) 
  10  FORMAT(/,2x, 
     &         'Reading User Supplied Fixed P/P Rates.') 
C 
      READ(INPREC(21:),*,ERR=95,END=96) FPPRUSER 
C 
C  Check that the 4 rates sum to 1.0. If they don't, warn the user. 
C 
      SUM=0. 
      DO I=1,4 
        SUM=SUM+FPPRUSER(I) 
      END DO 
C 
      IF(ABS(SUM-1.0).GT.0.0001) THEN 
        INERR=INERR+1 
        CALL QUITER(SUM,0,552,INERR) 
        GOTO 95 
      ENDIF 
C 
      FPPRFLAG = 1 
C 
      GOTO 99 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 1 
      GOTO 98 
   96 RC = -1 
C 
   98 CALL QUITER(0.,0,553,INERR) 
C 
   99 RETURN 
      END 
