      SUBROUTINE DBSELA(INPREC,INERR,RC) 
C 
C  DBSELA selects the AGE/MODEL YEAR INTERVALS that will appear in the database output. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELAGE,DB_MDLYR 
C    /YEARS4/    IYSTART, IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /IOUCOM/ IOUERR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C 
C  Notes: 
C 
      USE DATABASE, ONLY : DBSELAGE,DB_MDLYR,DBMODELYR 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'YEARS4.I' 
C 
C     Declare Paramter List 
      CHARACTER*(*), INTENT(IN) :: INPREC 
      INTEGER, INTENT(INOUT)    :: INERR 
      INTEGER, INTENT(OUT)      :: RC 
C 
C     Declare Local Variables 
      CHARACTER*18, SAVE :: MWARN = 'Data Base Warning:' 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
      CHARACTER*5, SAVE         :: LOWER='Lower' 
      CHARACTER*5, SAVE         :: UPPER='Upper' 
      CHARACTER*2, SAVE         :: ZERO=' 0' 
      CHARACTER*2, SAVE         :: TWENTYFOUR='24' 
      INTEGER, SAVE             :: UPPER_LMT=2 
      INTEGER, SAVE             :: LOWER_LMT=1 
      INTEGER, PARAMETER        :: IYES=1 
      INTEGER, PARAMETER        :: INO=0 
      INTEGER                   :: IOS 
      REAL                      :: TMP 
C 
C ..................................................................... 
C 
      RC = 0 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) DBSELAGE(LOWER_LMT), 
     &                                             DBSELAGE(UPPER_LMT) 
C 
      IF(DBSELAGE(LOWER_LMT) > DBSELAGE(UPPER_LMT)) THEN 
        TMP=DBSELAGE(LOWER_LMT) 
        DBSELAGE(LOWER_LMT)=DBSELAGE(UPPER_LMT) 
        DBSELAGE(UPPER_LMT)=TMP 
      ENDIF 
C 
      SELECT CASE(DB_MDLYR) 
C 
        CASE(INO) 
          IF(DBSELAGE(LOWER_LMT) < 0.OR. 
     &       DBSELAGE(LOWER_LMT) > 24) THEN 
            WRITE(IOUOUT,201) MWARN 
            WRITE(IOUOUT,202) DBSELAGE(LOWER_LMT),IOOFBD,LOWER,ZERO 
            DBSELAGE(LOWER_LMT) = 0 
          ENDIF 
          IF(DBSELAGE(UPPER_LMT) > 24) THEN 
            WRITE(IOUOUT,201) MWARN 
            WRITE(IOUOUT,202) DBSELAGE(UPPER_LMT),IOOFBD,UPPER, 
     &                        TWENTYFOUR 
            DBSELAGE(UPPER_LMT) = 24 
          ENDIF 
          GOTO 99 
C 
        CASE(IYES) 
          DBMODELYR(UPPER_LMT) = DBSELAGE(UPPER_LMT) 
          DBMODELYR(LOWER_LMT) = DBSELAGE(LOWER_LMT) 
          IF(DBSELAGE(UPPER_LMT) - DBSELAGE(LOWER_LMT) > 25) THEN 
            WRITE(IOUOUT,201) MWARN 
            WRITE(IOUOUT,204) DBSELAGE(LOWER_LMT),DBSELAGE(UPPER_LMT) 
          ENDIF 
C 
          IF(DBSELAGE(LOWER_LMT) < (IYSTART-24) .OR. 
     &       DBSELAGE(LOWER_LMT) > IYEND) THEN 
            WRITE(IOUOUT,201) MWARN 
            WRITE(IOUOUT,203) DBSELAGE(LOWER_LMT),IOOFBD, 
     &                        (IYSTART-24),IYEND 
            DBSELAGE(LOWER_LMT) = (IYSTART-24) 
          ENDIF 
          IF(DBSELAGE(UPPER_LMT) > IYEND) THEN 
            WRITE(IOUOUT,201) MWARN 
            WRITE(IOUOUT,203) DBSELAGE(UPPER_LMT),IOOFBD, 
     &                        IYSTART,IYEND 
            DBSELAGE(UPPER_LMT) = IYEND 
          ENDIF 
C 
          DBMODELYR(UPPER_LMT) = DBSELAGE(UPPER_LMT)    ! save model year info. for multiple scenarios 
          DBMODELYR(LOWER_LMT) = DBSELAGE(LOWER_LMT)    ! save model year info. for multiple scenarios 
      END SELECT 
c 
      GOTO 99 ! fatal read errors should have been caught 
              ! in READ statement; skip to RETURN 
C 
  201 FORMAT(T3,A18) 
  202 FORMAT(15X,I2,5A4,'vehicle age (0 to 24).'/ 
     &       15X,A5,' limit set to ',A2,'.'/) 
  203 FORMAT(15X,I4,5A4,'min/max model years (',I4,' to ',I4,')'/) 
  204 FORMAT(15X,'User entered database model year range'/ 
     &       15X,'(',I4,' to ',I4,') is greater than 25 years.'/) 
C 
C  Error or EOF on any attempted read. 
C 
   97 WRITE(IOUOUT,215) INPREC(1:19) 
      WRITE(IOUERR,215) INPREC(1:19) 
  215 FORMAT( 
     *'*** Error:  Please check data entered in ',A19,' command') 
      RC = IOS 
      INERR = INERR + 1 
      GOTO 99 
C 
   98 RC = -1 
      INERR = INERR + 1 
C 
   99 RETURN 
      END 
