      SUBROUTINE GETSPD(INPREC,INERR,RC) 
C 
C    GETSPD reads in the user supplied vehicle 
C    type START PER DAY array. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, NXTREC, UCCOMP, QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C 10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C           Also includes reading each set of day type (weekday, weekend) and 
C           vehicle class one at a time. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  30 May 00 @EP-DJB revised for start distribution input 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  31 Aug 99 @DynTel-ddj 2-699, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /STDIST/ STARTSPERDAY 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'STDIST.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     LABB 
C 
      LOGICAL     EXTDATA 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      CHARACTER   M6REC*150 
      INTEGER     SECTION 
      INTEGER     SOURCE 
      CHARACTER   STDLAB*40 
      INTEGER     IV,IWEEK,JDX,L 
      REAL        DUM(25), DUMM 
C 
      RC = 0 
      SPDFLG = SPDFLG + 1 
      IF (SPDFLG.GT.1) THEN 
C 
C     Error: Starts/day distribution has already been read in. 
C     
        WRITE(IOUERR,296) 
  296 FORMAT(/'*** ERROR: Starts/day Distribution', 
     &        'has already been read for this run.', 
     &       /'*** Only one Starts/day Distribution ', 
     &        'can be entered per run.',/, 
     &        '*** The distribution entered here will be ignored.  ') 
 
        RC = 3 
        INERR = INERR + 1 
      ELSE 
C 
C     Find out whether we are reading from an external data file. 
C 
        IF (INPREC(20:20).EQ.'@') THEN 
          EXTDATA = .FALSE. 
          SOURCE = IOUGEN 
          M6REC = INPREC 
        ELSE 
C 
C       Read the data from an external data file. 
C 
          EXTDATA = .TRUE. 
          SOURCE = IOUALT 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
          CALL CHKLAB('STARTS PER DAY',STDLAB,LABB,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
          LABEL = STDLAB(1:1) 
          LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs          CALL UCCOMP(INPREC(21:), FILENAME) 
          CALL LJCOMP(INPREC(21:), FILENAME)

          OPEN(SOURCE,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *         IOSTAT=RC,ERR=80) 
          IF (RC.NE.0) THEN 
C 
C     Error opening the external data file. 
C 
   80       WRITE(IOUERR,280) FILENAME 
  280       FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
            INERR = INERR + 1 
          ELSE 
            WRITE(IOUREP,200) FILENAME 
  200       FORMAT(/'* Reading start Starts/day distribution from the', 
     &              ' following external',/'* data file: ',A) 
C 
C       Read the required header line from the file. 
C 
            CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
            IF (RC.NE.0) THEN 
C 
C     Error reading the header of the external data file. 
C 
              SELECT CASE (RC) 
                CASE (:-1) 
                  WRITE(IOUERR,285) FILENAME 
  285             FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
                  RC = 1 
                  INERR = INERR + 1 
                CASE (1:) 
                  WRITE(IOUERR,286) FILENAME 
  286             FORMAT(/'*** ERROR reading the header of the ', 
     &                  'following external data file:',/'*** ',A) 
                CASE DEFAULT 
              END SELECT 
            ELSE 
C 
C       Check for a valid labeled input record. 
C 
              CALL UCCOMP(M6REC(1:19),LABEL) 
              CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
              IF (LABTYPE.NE.LABB) THEN 
C 
C     Invalid label on the header record. 
C 
                WRITE(IOUERR,290) M6REC 
  290           FORMAT(/'*** ERROR: the following record is not ', 
     &          'allowed in the header of the external data file:', 
     &         /'*** ',A) 
                RC = 2 
                INERR = INERR + 1 
              ELSE 
C 
C     Read the data field from the input record. 
C     Find the first data line. 
C 
                DO L=1,500  !up to 500 user input lines of data. 
                  CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
                  IF (RC.LT.0) THEN 
                    RC=0 
                    EXIT    !no more data to read. 
                  ENDIF 
                  BACKSPACE(SOURCE) 
C 
C      Read in user supplied data for starts per day. 
C 
                  DO JDX=1,25 
                    DUM(JDX)=0.0 
                  END DO 
                  READ (IOUALT,*,IOSTAT=RC) 
     *                 IWEEK, IV, (DUM(JDX),JDX=1,25) 
C 
C     Check the entered values for starts per day. 
C     Only values between 0.0 and 100.0 are allowed. 
C 
                  DO JDX=1,25 
                    IF(DUM(JDX).LT.0.0 .OR. DUM(JDX).GT.100.0) THEN 
                      DUMM=DUM(JDX) 
                      CALL QUITER(DUMM,0,183,INERR) 
                      EXIT 
                    ENDIF 
                  END DO 
C 
C     Check the input parameters and store the results. 
C 
                  IF(IV.GE.1 .AND. IV.LE.28) THEN 
                    IF(IWEEK.GE.1 .AND. IWEEK.LE.2 ) THEN 
                      DO JDX = 1,25 
                         STARTSPERDAY(JDX,IV,IWEEK) = DUM(JDX) 
                      END DO 
C 
C     One of the input paramters is out of range. 
C 
                    ELSE 
                      CALL QUITER(0.0,IWEEK,186,INERR) 
                    ENDIF 
                  ELSE 
                    CALL QUITER(0.0,IV,25,INERR) 
                  ENDIF 
                  IF(RC.LT.0) EXIT 
                END DO  !reading user inputs 
                IF (EXTDATA) CLOSE(IOUALT) 
              END IF 
            END IF 
          END IF 
        END IF 
      END IF 
      RETURN 
      END 
