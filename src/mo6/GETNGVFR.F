      SUBROUTINE GETNGVFR(INPREC,INERR,RC) 
C 
C 
C  GETNGVFR reads user supplied NGV by model year phase-in fractions for 
C  vehicle types 1 through 28. 
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
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  23 Feb 00 @EPA-elg 1-021, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C 
C  Common block array Subscripts: 
C 
C 
C 
C  Local variable Subscripts: 
C 
C    NGVPHIN       NGVPHIN(MY, IV) 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   NGVPHIN    R    Temporary array to store user supplied Tier2 bin 
C                   fractions prior to validation of input. 
C   BIN_TOTAL  R    Sum of 12 Tier2 bin fractions for a given year. 
C   BIN_DIFF   R    Difference between BIN_TOTAL and 1.0. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'P94INF.I' 
      INCLUDE 'NGV.I' 
C 
      INTEGER    INERR 
      CHARACTER  INPREC*(*) 
      INTEGER    RC 
C 
      CHARACTER  FILENAME*80 
      INTEGER    IV,IMY 
      INTEGER    IYR 
      CHARACTER  LABEL*40 
      INTEGER    LABTYPE 
      INTEGER    LDGIMP 
      CHARACTER  M6REC*150 
      INTEGER    NGVMAX 
      INTEGER    SECTION 
      CHARACTER  STDLAB*40 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      NGVMAX = 57 
      RC = 0 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('NGV FRACTION',STDLAB,LDGIMP,SECTION) 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file. 
C 

ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)

      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) TRIM(FILENAME) 
  120 FORMAT(/'* Reading NGV IMPLEMENTATION SCHEDULE ', 
     &        'from the following external',/'* data file: ',A) 
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
C     Read the file containing the phase-in schedule. Read one data 
C     block (for one vehicle type) at a time, skipping any blank lines 
C     or comments between the data blocks. 
C 
C       Find the first data line in the next data block. 
C 
        WRITE(IOUREP,300) 
  300   FORMAT(/,2x, 
     &         'Reading User Supplied NGV phase-in fractions') 
C 
      DO IV = 1,28 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) EXIT 
        BACKSPACE(IOUALT) 
C 
C       Read the implementation schedule for the current vehicle type. 
C 
        READ (IOUALT,*,IOSTAT=RC) 
     &         (NGVPHIN(IMY,IV),IMY=1,NGVMAX) 
C 
        IF(RC.LT.0) EXIT 
C 
      END DO 
C 
C     Implementation array items must be between 0.0 and 1.0. The phase- 
C     in values for all of the standards categories must sum to 1.0 in 
C     each calendar year. 
C 
C     IYR of 1 equals the 1994 model year.  This is implemented later in NGV_EF. 
C 
C 
      DO IV = 1,28 
        DO IYR = 1, NGVMAX 
           NGVPHIN(IYR,IV) = NGVPHIN(IYR,IV) / 100. 
           IF (NGVPHIN(IYR,IV).LT.0.0 .OR. 
     *         NGVPHIN(IYR,IV).GT.1.0) THEN 
              CALL QUITER(NGVPHIN(IYR,IV),0,535,INERR) 
              RETURN 
           END IF 
        END DO 
      END DO 
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
      GOTO 97 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 2 
      GOTO 97 
C 
  97  CLOSE(IOUALT) 
  98  INERR = INERR + 1 
  99  RETURN 
      END 
