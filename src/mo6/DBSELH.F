      SUBROUTINE DBSELH(INPREC,INERR,RC) 
C 
C  DBSELE selects the HOURLY INTERVALS that will appear in the database output. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELHR 
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
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C 
C  Notes: 
C 
      USE DATABASE, ONLY : DBSELHR 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
C     Declare parameter list 
      CHARACTER*(*), INTENT(IN) :: INPREC 
      INTEGER, INTENT(INOUT)    :: INERR 
      INTEGER, INTENT(OUT)      :: RC 
C 
C     Declare local variables 
      CHARACTER*18, SAVE :: MWARN = 'Data Base Warning:' 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
      CHARACTER*5, SAVE         :: LOWER='Lower' 
      CHARACTER*5, SAVE         :: UPPER='Upper' 
      CHARACTER*2, SAVE         :: ONE=' 1' 
      CHARACTER*2, SAVE         :: TWENTYFOUR='24' 
      INTEGER, SAVE             :: UPPER_LMT=2 
      INTEGER, SAVE             :: LOWER_LMT=1 
      INTEGER                   :: IOS 
      REAL                      :: TMP 
C 
C ..................................................................... 
C 
      RC = 0 
C 
      READ(INPREC(21:),*,IOSTAT=IOS,ERR=97,END=98) DBSELHR(LOWER_LMT), 
     &                                             DBSELHR(UPPER_LMT) 
C 
      IF(DBSELHR(LOWER_LMT) > DBSELHR(UPPER_LMT)) THEN 
        TMP=DBSELHR(LOWER_LMT) 
        DBSELHR(LOWER_LMT)=DBSELHR(UPPER_LMT) 
        DBSELHR(UPPER_LMT)=TMP 
      ENDIF 
C 
      IF (DBSELHR(LOWER_LMT) < 1.OR.DBSELHR(LOWER_LMT) > 24) THEN 
        WRITE(IOUOUT,201) MWARN 
        WRITE(IOUOUT,202) DBSELHR(LOWER_LMT),IOOFBD,LOWER,ONE 
        DBSELHR(LOWER_LMT) = 1 
      ENDIF 
      IF (DBSELHR(UPPER_LMT) > 24) THEN 
        WRITE(IOUOUT,201) MWARN 
        WRITE(IOUOUT,202) DBSELHR(UPPER_LMT),IOOFBD,UPPER, 
     &                        TWENTYFOUR 
        DBSELHR(UPPER_LMT) = 24 
      ENDIF 
C 
  201 FORMAT(T3,A18) 
  202 FORMAT(15X,I2,5A4,'hourly range (1 to 24).'/ 
     &       15X,A5,' limit set to ',A2,'.'/) 
C 
      GOTO 99 ! fatal read errors should have been caught 
              ! in READ statement; skip to RETURN 
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
