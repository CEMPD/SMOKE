      SUBROUTINE GETATWAV(INPREC,INERR,RC) 
C 
C  GETATWAV reads the Air Toxic RVP Waiver flag for a Run. 
C  (This flag is exactly the same as the RVP WAIVER flag 
C   found in the OXYGENATEF FUELS command.) 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER, UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  03 Jul 01 AIR Task 15, new subroutine 
C 
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
C    /OXY1/ IGASHW,GLIGAS 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'OXY1.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
C 
      INTEGER IOS 
       
      RC = 0 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) IGASHW 
C 
      GLIGAS = IGASHW  !GETATWAV is a Run Level Command 
C      
      IF (IGASHW.LT.1.OR.IGASHW.GT.2) THEN 
        CALL QUITER(1.,IGASHW,575,INERR) 
        IGASHW = 1 
        GLIGAS = 1 
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
