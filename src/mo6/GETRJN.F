      SUBROUTINE GETRJN(INPREC,INERR,RC) 
C 
C  GETRJN reads the region flag for a scenario. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER, UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  28 Dec 00 @EPA-bag bug294 QUITER message for IREJN out of bounds is now 
C            message #24 
C  31 May 00 @EPA-elg Changed Read to be Free Format. 
C  22 Jun 99 @DynTel-MLA 2-694  Region now must be either 1 or 2. Regions 
C            3 and 4 are no longer valid. 
C  22 Oct 1998 @DynTel-ZK 2-000 Explicit Typing  
C  Input on call: 
C 
C    argument list: INPREC 
C 
C    common blocks: 
C    /REGION/ FEET 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C 
C    common blocks: 
C    /REGION/ IREJN,ALT 
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
C  24 Jun 98 @DynTel-bag 2-000 Changed QUITER argument RVALUE ( --> 1)  
C            to reflect the changes made in the NAMFLG array. 
C 
      IMPLICIT NONE 
      INCLUDE 'REGION.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
C 
      INTEGER IOS 
       
      RC = 0 
C 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) IREJN 
C 
C      
      IF (IREJN.LT.1.OR.IREJN.GT.2) THEN 
        CALL QUITER(1.,IREJN,24,INERR) 
        IREJN = 1 
        RC = 1 
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
