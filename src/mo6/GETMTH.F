      SUBROUTINE GETMTH(INPREC,INERR,RC) 
C 
C  GETMTH reads the evaluation month for a scenario. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  25 May 2000 @EPA-elg  Changed Fixed format input to Free Format. 
C  16 OCT 1998 @DynTel-ZL 2-000 Explicit Typing  
C  Input on call: 
C 
C    argument list: INPREC 
C 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C 
C    common blocks: 
C    /EVAL/ MEVAL 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DATA     C    Holds the data value from the input record. 
C 
C  Notes: 
C 
C  27 Mar 98 @DynTel-MLA 2-659  New routine for Mobile6 (the transitional 
C            labeled inputs version). This subroutine sets the RFGFLG. 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAL.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     IOS 
C 
      RC = 0 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) MEVAL 
C 
      IF(MEVAL.LE.0) MEVAL=1 
      IF(MEVAL.NE.1.AND.MEVAL.NE.7) THEN  
        CALL QUITER(0.,MEVAL,145,INERR) 
        RC = 2 
      ENDIF 
C      
      GOTO 99 
C 
C  Error or EOF on any attempted read. 
C 
   97 RC = IOS 
      INERR = INERR + 1 
      GOTO 99 
C 
   98 RC = -1 
      INERR = INERR + 1 
      GOTO 99 
C 
   99 RETURN 
      END 
