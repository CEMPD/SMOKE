      SUBROUTINE GETHSACT(INPREC,INERR,RC) 
C 
C 
C  GETHSACT reads user supplied hot soak activity fractions for 14 
C  hourly periods (6-7AM to 7-8PM) and 1-60 minutes. 
C 
C 
C 
C  Called by: DSPTCH. 
C 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  31 jan 01 bag          long form of command (HOT SOAK ACTIVITY) used 
C                         to check the input file header 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  11 Oct 00 bag          removed superfluous GOTO 98 
C  21 Sep 00 bag bug232   routine did not close data file after sucessfull read 
C  12 Jul 00 AIR Task 06, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C    /IOUCOM/ IOUREP 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /EVAPACT/ HSACTUSER, HSACTFLAG 
C 
C  Common block array Subscripts: 
C 
C 
C 
C  Local variable Subscripts: 
C 
C    DUMMY       DUMMY(MY, IV) 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store user supplied hot soak  
C                   activity fractions prior to validation of input. 
C   SUM        R    Sum of 60 by-miniute hot activity fractions for a 
C                   given hourly period. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAPACT.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
C  
      INTEGER, EXTERNAL :: JUNIT
C
      INTEGER    INERR 
      CHARACTER  INPREC*(*) 
      INTEGER    RC 
C 
      CHARACTER  FILENAME*80 
      CHARACTER  LABEL*40 
      INTEGER    LABTYPE 
      INTEGER    LDGIMP 
      CHARACTER  M6REC*150 
      INTEGER    SECTION 
      CHARACTER  STDLAB*40 
      INTEGER    IMIN,IH,I 
      REAL       DUMMY(14) 
      REAL       SUM 
      LOGICAL :: FOPEN 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('HOT SOAK ACTIVITY',STDLAB,LDGIMP,SECTION) 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file. 
C 

ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)
 
      IOUALT = JUNIT()

      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) TRIM(FILENAME) 
  120 FORMAT(/20X,' User supplied hot soak activity: ',A) 
C 
C     Read the required header line from the file. 
C 
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC.NE.0) THEN 
        GOTO 85 
      ENDIF 
C 
C     Check for a valid labeled input record. 
C 
      CALL UCCOMP(M6REC(1:19),LABEL) 
      CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
      IF (LABTYPE.NE.LDGIMP) THEN 
        GOTO 90 
      ENDIF 
C 
C     Initialize the user supplied hot soak array to zero. 
C 
      DO IH=1,14 
        DO IMIN=1,60 
          HSACTUSER(IMIN,IH)=0.0 
        END DO 
      END DO 
C       
C     Read the file containing the hot soak activity fractions. Read one data 
C     line (by minute over 14 hour periods) at a time, skipping any blank lines 
C     or comments between the data blocks. 
C 
   75 WRITE(IOUREP,300) 
  300 FORMAT(/,2x, 
     &         'Reading User Supplied Hot Soak Activity Fractions') 
C 
C       Find the first data line in the next data block. 
C 
      DO I = 1,60 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read the hot soak activty fractions for the designated minute 
C       across all 14 hourly periods. Perform QC on the data. 
C 
        READ (IOUALT,*,IOSTAT=RC) IMIN, DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
        IF (IMIN.LT.1 .OR. IMIN.GT.60) THEN 
        INERR=INERR+1 
        CALL QUITER(0.,IMIN,548,INERR) 
        RC=1 
        EXIT 
        ENDIF 
C 
        DO IH=1,14 
          HSACTUSER(IMIN,IH)=DUMMY(IH) 
        END DO 
C 
      END DO 
C 
C     Verify that the total for each hourly period, over all 60 minutes, 
C     sums to 1.0 +/- 0.00005. If not, warn user where the error occurred. 
C 
      DO IH=1,14 
        SUM=0.0 
        DO IMIN=1,60 
          SUM=SUM+HSACTUSER(IMIN,IH) 
        END DO 
        IF(ABS(SUM-1.0).GT.0.0001) THEN 
        INERR=INERR+1 
        CALL QUITER(SUM,IH,549,INERR) 
        RC=1 
        EXIT 
        ENDIF 
      END DO 
C 
      IF(RC.EQ.0) HSACTFLAG = 1 
C 
      GOTO 99 
C 
C     Close the datafile and write the appropriate warning message. 
C 
  80  WRITE(IOUOUT,180,ERR=98) TRIM(FILENAME) 
      WRITE(IOUERR,180,ERR=98) TRIM(FILENAME) 
 180  FORMAT('*** Error opening external data file ',A80) 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) TRIM(FILENAME) 
  185     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186) TRIM(FILENAME) 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 2 
C 
C 
  98  INERR = INERR + 1 
  99  INQUIRE(IOUALT,OPENED=FOPEN) 
      IF (FOPEN) CLOSE(IOUALT) 
      RETURN 
      END 
