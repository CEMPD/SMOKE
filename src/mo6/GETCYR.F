      SUBROUTINE GETCYR(INPREC,ICY,INERR,RC) 
C 
C  GETCYR reads the calendar year record for a scenario. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C 19 Oct 00 @EPA-bag Added data checks for database model years 
C  25 May 00 @EPA-elg  Changed Fixed format input to Free Format. 
C  07 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: /YEARS4/ IYEND 
C 
C  Output on return: 
C 
C    argument list: ICY, INERR, RC 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
C  27 Mar 98 @DynTel-MLA 2-659  New routine for Mobile6 (the transitional 
C            labeled inputs version). This subroutine reads the calendar 
C            year for a scenario. 
C 
      USE DATABASE, ONLY : DBSELAGE,DB_MDLYR,DBMODELYR 
      IMPLICIT NONE 
      INCLUDE 'YEARS4.I' 
      INCLUDE 'IOUCOM.I' 
C 
C     Declare parameter list 
      INTEGER,   INTENT(OUT) :: ICY 
      INTEGER, INTENT(INOUT):: INERR 
      INTEGER,   INTENT(OUT):: RC 
      CHARACTER, INTENT(IN) :: INPREC*(*) 
C 
C     Declare local variables 
      INTEGER IOS 
      CHARACTER*18, SAVE :: MWARN = 'Data Base Warning:' 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
      CHARACTER*5, SAVE         :: LOWER='Lower' 
      CHARACTER*5, SAVE         :: UPPER='Upper' 
      INTEGER, SAVE :: LOWER_LMT=1 
      INTEGER, SAVE :: UPPER_LMT=2 
      INTEGER, SAVE :: IYES=1 
      REAL          :: R_ICY 
      INTEGER       :: DIFF 
      INTEGER       :: ITEMPAGE 
C 
      RC = 0 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) R_ICY 
      ICY=INT(R_ICY) 
C 
      IF (ICY.LT.IYSTART.OR.ICY.GT.IYEND) THEN 
        CALL QUITER(0.0,ICY,27,INERR) 
        RC = 1 
      ENDIF 
C 
      IF(DB_MDLYR==IYES) THEN 
C 
        DBSELAGE(UPPER_LMT) = DBMODELYR(UPPER_LMT)    ! for multiple scenarios get back model year 
        DBSELAGE(LOWER_LMT) = DBMODELYR(LOWER_LMT)    ! for multiple scenarios get back model year 
        ITEMPAGE = DBSELAGE(LOWER_LMT) 
        DIFF=ICY - DBSELAGE(UPPER_LMT) 
C 
        DBSELAGE(LOWER_LMT) = MAX(0,DIFF) 
        IF(DIFF<0) THEN 
          WRITE(IOUOUT,201) MWARN 
          WRITE(IOUOUT,202) DBSELAGE(UPPER_LMT),IOOFBD,ICY-24,ICY,UPPER, 
     &    ICY 
        ENDIF 
C 
        DIFF=ICY - ITEMPAGE 
        DBSELAGE(UPPER_LMT) = MIN(24,DIFF) 
        IF(DIFF>24) THEN 
          WRITE(IOUOUT,201) MWARN 
          WRITE(IOUOUT,202) ITEMPAGE,IOOFBD,ICY-24,ICY,LOWER,ICY-24 
        ENDIF 
C 
      ENDIF 
C 
  201 FORMAT(/T3,A18) 
  202 FORMAT(15X,I4,5A4,'vehicle model years (',I4,' to ',I4,').'/ 
     &       15X,A5,' limit set to ',I4,'.'/) 
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
