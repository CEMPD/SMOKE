      SUBROUTINE GETNGVEF(INPREC,INERR,RC) 
C 
C 
C  GETNGVEF reads user supplied NGV by model year emission factors for 
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
C  02 Oct 00 AIR Task 07: Bug Fix 268, revised error handling to 
C            properly CLOSE file.     
C  15 Jun 00 @EPA-bag bug150, calls to quiter for error #'s 
C             532,533, and 534 (vehicle type out of bounds) 
C             had IVALUE and RVALUE interchanged 
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
C    DUMMY       DUMMY(MY, IV) 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store user supplied Tier2 bin 
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
      CHARACTER  LABEL*40 
      INTEGER    LABTYPE 
      INTEGER    LDGIMP 
      CHARACTER  M6REC*150 
      INTEGER    SECTION 
      CHARACTER  STDLAB*40 
      INTEGER    IVC,IPC,ISRC,I 
      REAL       DUMZ,DUMD,DUMH 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('NGV EF',STDLAB,LDGIMP,SECTION) 
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
  120 FORMAT(/'* Reading NGV Emission Factors ', 
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
      DO I = 1,500 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read the implementation schedule for the current vehicle type. 
C 
        READ (IOUALT,*,IOSTAT=RC) 
     &          IVC, IPC, ISRC, DUMZ, DUMD, DUMH 
C 
        IF(RC.LT.0) EXIT 
C 
        IF (IVC.GT.28 .OR. IVC.LT.1) THEN 
          CALL QUITER(0.,IVC,532,INERR) 
        ELSEIF (IPC.GT.3  .OR. IPC.LT.1) THEN 
          CALL QUITER(0.,IVC,533,INERR) 
        ELSEIF (ISRC.GT.2 .OR. ISRC.LT.1) THEN 
          CALL QUITER(0.,IVC,534,INERR) 
        ENDIF 
C 
        NGVZML(IVC,IPC,ISRC)  = DUMZ 
        NGVDET(IVC,IPC,ISRC)  = DUMD 
        NGVHIGH(IVC,IPC,ISRC) = DUMH 
C 
      END DO 
C 
   75 WRITE(IOUREP,300) 
  300 FORMAT(/,2x, 
     &         'Reading User Supplied NGV Emission Factors') 
C 
      RC = 0 
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
  99  CLOSE(IOUALT) 
      RETURN 
      END 
