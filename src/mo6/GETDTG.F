      SUBROUTINE GETDTG(NSCEN,INERR,RC) 
C 
C  GETDTG reads in the detergent gasoline parameters record. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  31 May 00 @EPA-elg Changed Read to be Free Format. 
C  08 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /IOUCOM/ IOUGEN 
C 
C  Output on return: 
C 
C    common block: 
C    /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDEFF,FIDFRC,IVDEFF,IVDFRC 
C    /FLAGS1/ DTGFLG 
C 
C    INERR - Running tally of errors. 
C    RC    - The return code. 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ---------------------------------------------------- 
C  ERRVAL   R    used to store erroneous value found in detergent gasoline input 
C  IERRFL   I    used to store an integer that points to the name of the 
C                detergent gasoline input variable that is in error 
C 
C  Notes: 
C 
C  GETDTG is new for Mobile6. This code used to reside in RDCITY. 
C     
      IMPLICIT NONE 
      INCLUDE 'DTGGAS.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER   I 
      INTEGER   IERRFL 
      INTEGER   INERR 
      INTEGER   NSCEN          
      INTEGER   RC 
C 
      REAL      ERRVAL 
C 
      RC = 0 
      IF (DTGFLG.EQ.0) THEN 
        RC = 1 
        GOTO 99 
      ENDIF 
      DTGFLG = 2 
C 
C     set defaults 
C 
      DGPISY    = 1995 
      DGPIEY    = 1998 
      DGPHSF    = 0.50 
      FIDEFF(1) = 0.025 
      FIDEFF(2) = 0.0165 
      FIDEFF(3) = 0.0 
      FIDFRC    = 0.019 
      IVDEFF(1) = 0.054 
      IVDEFF(2) = 0.065 
      IVDEFF(3) = 0.118 
      IVDFRC    = 0.18 
C 
      READ (IOUGEN,*,ERR=95,END=96) 
     *     DGPIEY,DGPHSF,FIDFRC,(FIDEFF(I),I=1,3),IVDFRC, 
     *     (IVDEFF(I),I=1,3) 
C  130 FORMAT (I4,9(1X,F5.2)) 
      IF (DGPIEY.LE.1995) CALL QUITER(0.,DGPIEY,169,INERR) 
C 
      DGPHSF = DGPHSF / 100. 
      FIDFRC = FIDFRC / 100. 
      IVDFRC = IVDFRC / 100. 
      DO 81 I=1,3 
        FIDEFF(I) = FIDEFF(I) / 100. 
        IVDEFF(I) = IVDEFF(I) / 100. 
   81 CONTINUE 
C 
      ERRVAL = 0.0 
      IF (DGPHSF.LT.0.0.OR.DGPHSF.GT.1.0) THEN 
          ERRVAL = DGPHSF 
          IERRFL = 1 
      ELSEIF (FIDFRC.LT.0.0.OR.FIDFRC.GT.1.0) THEN 
          ERRVAL = FIDFRC 
          IERRFL = 2 
      ELSEIF (IVDFRC.LT.0.0.OR.IVDFRC.GT.1.0) THEN 
          ERRVAL = IVDFRC 
          IERRFL = 3 
      ELSEIF (FIDEFF(1).LT.0.0.OR.FIDEFF(1).GT.1.0) THEN 
          ERRVAL = FIDEFF(1) 
          IERRFL = 4 
      ELSEIF (FIDEFF(2).LT.0.0.OR.FIDEFF(2).GT.1.0) THEN 
          ERRVAL = FIDEFF(2) 
          IERRFL = 5 
      ELSEIF (FIDEFF(3).LT.0.0.OR.FIDEFF(3).GT.1.0) THEN 
          ERRVAL = FIDEFF(3) 
          IERRFL = 6 
      ELSEIF (IVDEFF(1).LT.0.0.OR.IVDEFF(1).GT.1.0) THEN 
          ERRVAL = IVDEFF(1) 
          IERRFL = 7 
      ELSEIF (IVDEFF(2).LT.0.0.OR.IVDEFF(2).GT.1.0) THEN 
          ERRVAL = IVDEFF(2) 
          IERRFL = 8 
      ELSEIF (IVDEFF(3).LT.0.0.OR.IVDEFF(3).GT.1.0) THEN 
          ERRVAL = IVDEFF(3) 
          IERRFL = 9 
      ENDIF 
      IF (ERRVAL.NE.0.0) CALL QUITER(ERRVAL,IERRFL,70,INERR) 
C 
      IF (NSCEN.EQ.0) THEN 
        GLDTGFLG = 2 
        GLDGPSY = 1995 
        GLDGPEY = DGPIEY 
        GLDGPHSF = DGPHSF 
        GLFIDEF(1) = FIDEFF(1) 
        GLFIDEF(2) = FIDEFF(2) 
        GLFIDEF(3) = FIDEFF(3) 
        GLFIDFC = FIDFRC 
        GLIVDEF(1) = IVDEFF(1) 
        GLIVDEF(2) = IVDEFF(2) 
        GLIVDEF(3) = IVDEFF(3) 
        GLIVDFC = IVDFRC 
      ENDIF 
C 
      GOTO 99 
C 
C  Error or EOF on any attempted read => take alternate 
C  return 1 => run aborts. 
C 
   95 RC = 1 
      GOTO 98 
   96 RC = -1 
C 
   98 CALL QUITER(0.,0,69,INERR) 
C 
   99 RETURN 
      END 
