      SUBROUTINE GETWIN(INPREC,NSCEN,INERR,RC) 
C 
C  GETWIN sets the WINFLG.  Sets the program to winter or summer. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER, UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C  29 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    argument list: INPREC 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C    common blocks: 
C    /RFORM1/ WINFLG 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DATA     C    Holds the data value from the input record. 
C  X        R    Temporary real to integer variable 
C 
C  Notes: 
C 
C  27 Mar 98 @DynTel-MLA 2-659  New routine for Mobile6 (the transitional 
C            labeled inputs version). This subroutine sets the TEMFLG. 
C  24 Jun 98 @DynTel-bag 2-000 Changed QUITER argument RVALUE (29 --> 8)  
C            to reflect changes made in the NAMFLG array. 
C             
      IMPLICIT NONE 
      INCLUDE 'RFORM1.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     NSCEN 
      INTEGER     RC 
C 
      CHARACTER   DATA*5 
      INTEGER I   OS 
      REAL        X 
C 
      RC = 0 
C 
      CALL UCCOMP(INPREC(21:), DATA) 
      READ(DATA,*,IOSTAT=IOS,ERR=97,END=98) X 
      WINFLG=INT(X) 
C      
      IF (WINFLG.LT.0.OR.WINFLG.GT.2) THEN 
        CALL QUITER(8.,WINFLG,1,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
      IF (NSCEN.EQ.0) THEN 
        GLWIN = WINFLG 
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
C 
   99 RETURN 
      END 
