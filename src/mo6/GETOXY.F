      SUBROUTINE GETOXY(INPREC,NSCEN,INERR,RC) 
C 
C  GETOXY reads in the oxygenated fuels data. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  15 Jan 2002 AIR Task 15: Changed maximum alcohol oxy limit from 3.5% to 3.7% 
C  16 Nov 2000 @EPA-bag bug302 copy SHRMKT(1) to GLSHARE(1) global/run level 
C              save variables 
C  31 May 2000 @EPA-elg  Fixed READ Format Changed to Free READ format. 
C  20 Oct 1998 @DynTel-ZK 2-000 Explicit Typing  
C 
C  Input on call: 
C 
C    argument list: NSCEN 
C    common blocks: 
C    /IOUCOM/ IOUGEN 
C 
C  Output on return: 
C 
C   argument list: INERR,INPREC,RC 
C    common block: 
C    /FLAGS1/ GLOXFLAG,OXYFLG 
C    /OXY1/   GLIGAS,GLOXCNT,GLSHARE,IGASHW,OXYCNT,SHRMKT 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ---------------------------------------------------- 
C  MKTSHR   I 
C  NDP      I 
C  RVAL     R 
C  NVAL     I 
C 
C  Notes: 
C 
C  GETOXY is new for Mobile6. This code used to reside in RDCITY. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'OXY1.I' 
C 
      CHARACTER INPREC*(*) 
      INTEGER   INERR 
      INTEGER   NSCEN 
      INTEGER   RC 
      INTEGER   MKTSHR, NDP, IS, NVAL 
      REAL      RVAL  
C      
      RC = 0 
C 
C  Obtain and verify the oxy fuel inputs. These might be global values 
C  or values that are local to the current scenario. 
C 
      OXYFLG = 2 
C 
      READ(INPREC(21:),*,ERR=95,END=96) SHRMKT(2),SHRMKT(3), 
     *        OXYCNT(1),OXYCNT(2),IGASHW 
C 
C     Check the parameters for oxygenated fuels.  
C 
      IF(SHRMKT(2).LT.0..OR.SHRMKT(2).GT.1.0) 
     *    CALL QUITER(SHRMKT(2),0,2,INERR) 
      IF(SHRMKT(3).LT.0..OR.SHRMKT(3).GT.1.0) 
     *    CALL QUITER(SHRMKT(3),0,3,INERR) 
C 
C  Turn SHRMKT(IS) into NDP significant digits 
C 
      MKTSHR=0 
      NDP=3 
      DO 80 IS=2,3 
         RVAL=SHRMKT(IS)*(10**NDP) 
         NVAL=INT(RVAL) 
         MKTSHR=MKTSHR+NVAL 
         RVAL=1.0*RVAL/(10**NDP) 
         SHRMKT(IS)=RVAL 
   80 CONTINUE 
      IF(MKTSHR.GT.1000) CALL QUITER(SHRMKT(2)+SHRMKT(3),0,4,INERR) 
      SHRMKT(1)=1.0-SHRMKT(2)-SHRMKT(3) 
      SHRMKT(1)=(NINT(SHRMKT(1)*1000))/1000. 
C 
      IF(OXYCNT(1).LT.0..OR.OXYCNT(1).GT..027001) 
     *    CALL QUITER(OXYCNT(1),0,5,INERR) 
      IF(OXYCNT(2).LT.0..OR.OXYCNT(2).GT..037001) 
     *    CALL QUITER(OXYCNT(2),0,6,INERR) 
C 
      IF(IGASHW.LT.1.OR.IGASHW.GT.2) CALL QUITER(0.,IGASHW,7,INERR) 
C 
      IF(SHRMKT(2).GT.0..AND.OXYCNT(1).EQ.0. .OR. 
     *   SHRMKT(2).EQ.0..AND.OXYCNT(1).GT.0.) 
     *   CALL QUITER(0.,0,8,INERR) 
      IF(SHRMKT(3).GT.0..AND.OXYCNT(2).EQ.0. .OR. 
     *   SHRMKT(3).EQ.0..AND.OXYCNT(2).GT.0. .OR. 
     *   SHRMKT(3).EQ.0..AND.IGASHW.EQ.2) 
     *   CALL QUITER(0.,0,9,INERR) 
C 
      IF (NSCEN.EQ.0) THEN 
        GLOXFLAG = 2 
        GLSHARE(1) = SHRMKT(1) 
        GLSHARE(2) = SHRMKT(2) 
        GLSHARE(3) = SHRMKT(3) 
        GLOXCNT(1) = OXYCNT(1) 
        GLOXCNT(2) = OXYCNT(2) 
        GLIGAS = IGASHW 
      ENDIF 
C 
      GOTO 99 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 1 
      GOTO 98 
   96 RC = -1 
C 
   98 CALL QUITER(0.,0,137,INERR) 
C 
   99 RETURN 
      END 
