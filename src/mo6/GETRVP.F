      SUBROUTINE GETRVP(INPREC,NSCEN,INERR,RC) 
C 
C  GETRVP reads in the fuel RVP value for a scenario. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  25 May 2000 @EPA-elg  Changed Fixed format input to Free Format. 
C  22 Oct 1998 @DynTel-ZK 2-000 Explict Typing 
C  Input on call: 
C 
C    common blocks: 
C 
C    /CITRV1/  RVPBAS 
C    /RVPNAT/  RVPLIM 
C 
C  Output on return: 
C 
C    common block: 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ---------------------------------------------------- 
C 
C  Notes: 
C 
C  GETRVP is new for Mobile6. This code used to reside in RDCITY. 
      IMPLICIT NONE  
      INCLUDE 'CITRV1.I' 
      INCLUDE 'RVPNAT.I' 
C 
      CHARACTER INPREC*(*) 
      INTEGER   INERR 
      INTEGER   RC 
      INTEGER   NSCEN 
      RC = 0 
       
C  Read and verify the fuel RVP. 
C 
      READ(INPREC(21:),*,ERR=95,END=96) RVPBAS 
C 
C  Check RVPs and re-assign if necessary. 
C 
      IF(RVPBAS.LT.RVPLIM(1)) CALL QUITER(RVPLIM(1),1,90,INERR) 
      IF(RVPBAS.LT.RVPLIM(1)) RVPBAS=RVPLIM(1) 
C 
      IF(RVPBAS.GT.RVPLIM(2)) CALL QUITER(RVPLIM(2),1,90,INERR) 
      IF(RVPBAS.GT.RVPLIM(2)) RVPBAS=RVPLIM(2) 
C 
      IF (NSCEN.EQ.0) THEN 
        GLRVPBAS = RVPBAS 
      ENDIF 
C 
      GOTO 99 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 1 
      INERR = INERR + 1 
      GOTO 99 
   96 RC = -1 
      INERR = INERR + 1 
      GOTO 99 
C 
   99 RETURN 
      END 
