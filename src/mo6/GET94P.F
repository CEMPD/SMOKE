      SUBROUTINE GET94P(INPREC,INERR,RC) 
C 
C  GET94P reads the LDG 94+ implementation array and LDG 94+ emissions 
C  coefficients for user-supplied emission factors. 
C 
C  Called by DSPTCH 
C 
C  Calls  ABS (FORTRAN library function), QUITER.  
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  25 Sep 00 bug257 Added IZFLG to account for 100% ZEV vehicles, and removed 
C            Fixed Format input requirement. 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  24 Apr 00 AIR Task 02: Fixed comments regarding parameter list. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  07 Jul 99 @DynTel-MLA 2-694  Standard Mobile6 external data file 
C            header is now required; read data in blocks by IV; allow 
C            blank lines and comments in the data file. 
C  22 Jun 99 @DynTel-MLA 2-694  Removed 94+ start and running emission 
C            factor inputs, which need to be rewritten for Mobile6. 
C  11 Mar 99 @DynTel-MLA 2-701  Changed size of P94NAM (file name) to 80. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C  03 Mar 99 @DynTel-MLA 2-659 Replaced call to FEXIST with call to OPEN. 
C  06 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C    argument list: INPREC, INERR, RC 
C    constant: MAXLYR 
C    common blocks: 
C    /BASEQ9/ MAXLYR 
C    /IOUCOM/ IOUREP,IOUOUT 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /BASEQ9/ LEVIMP   
C    /FLAGS1/ USR94P 
C    /P94INF/ P94AR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  FILENAME C    Name of external data file. 
C  ISTD     I    Loop control variable. 
C  IV       I    Loop control variable. 
C  IYR      I    Loop control variable. 
C  LABEL    C    The label field from an input record. 
C  LABTYPE  I    Token number corresponding to an input label. 
C  LDGIMP   I    Token number for input label '94+ LDG IMPLEMENTATION.' 
C  M6REC    C    Record from the external data file. 
C  SECTION  I    Dummy parameter for CHKLAB. 
C  SUM      R    Sum of phase-in values for each calendar year. 
C  STDLAB   C    Standard form of input label '94+ LDG IMPLEMENTATION.' 
C 
C  Notes: 
C 
C  April-20-1998 @Dyntel-ID GET94P was created to read and check the data 
C  in data files entered by user to phase-in LDG 94+ schemes.   
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'P94INF.I' 
C 
      INTEGER, EXTERNAL :: JUNIT
C
      INTEGER    INERR 
      CHARACTER  INPREC*(*) 
      INTEGER    RC 
C 
      CHARACTER  FILENAME*80 
      INTEGER    ISTD 
      INTEGER    IV 
      INTEGER    IYR 
      CHARACTER  LABEL*40 
      INTEGER    LABTYPE 
      CHARACTER  LASTREC*150 
      INTEGER    LDGIMP 
      CHARACTER  M6REC*150 
      INTEGER    SECTION 
      REAL       SUM 
      CHARACTER  STDLAB*40 
      INTEGER    IZFLG 
            
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('94+ LDG IMPLEMENT',STDLAB,LDGIMP,SECTION) 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file. 
C 

ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)
      
ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )
      M6FILE = TRIM( M6DIR ) // '/' // FILENAME
      
      IOUALT = JUNIT()
      
      OPEN(IOUALT,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80)
      WRITE(IOUREP,120) TRIM(FILENAME) 
  120 FORMAT(/'* Reading 94+ LEV IMPLEMENTATION SCHEDULE from the ', 
     &        'following external',/'* data file: ',A) 
C 
C     Store name of the file in P94AR array to print names  
C     of 94P+ phase-in data files. 
C 
      P94AR(1) = FILENAME 
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
      DO IV = 1, 5 
C 
C       Find the first data line in the next data block. 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) GOTO 95 
        BACKSPACE(IOUALT) 
C 
C       Read the implementation schedule for the current vehicle type. 
C 
        DO IYR = 1, MAXLYR 
          READ (IOUALT,*,IOSTAT=RC,ERR=95,END=95) 
     &              (LEVIMP(ISTD,IYR,IV),ISTD=1,11) 
        END DO 
C 
      END DO 
C 
      CLOSE(IOUALT) 
C 
C     Implementation array items must be between 0.0 and 1.0. The phase- 
C     in values for all of the standards categories must sum to 1.0 in 
C     each calendar year. 
C 
      IZFLG = 0 
      DO IV = 1, 5 
        DO IYR = 1, MAXLYR 
          SUM = 0.0 
          DO ISTD = 1, 11 
            IF (LEVIMP(ISTD,IYR,IV).LT.0.0.OR. 
     &          LEVIMP(ISTD,IYR,IV).GT.1.0) THEN 
              CALL QUITER(LEVIMP(ISTD,IYR,IV),0,175,INERR) 
            ELSE 
              SUM = SUM + LEVIMP(ISTD,IYR,IV) 
            ENDIF 
            IF(LEVIMP(11,IYR,IV).EQ.1.0) THEN 
               IZFLG = 1 
            ENDIF 
          END DO   
          IF (ABS(SUM - 1.0).GT.0.0001)THEN 
            CALL QUITER(SUM,0,176,INERR)                   
          ENDIF 
        END DO 
        IF(IZFLG.EQ.1) THEN 
          WRITE(IOUERR,187) IV 
  187     FORMAT(/' *** WARNING: ZEV Phase-in is set to 100%', 
     *            'for vehicle type: ',I2,/) 
          IZFLG = 0 
        ENDIF 
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
   95 IF (RC.LT.0) THEN 
        WRITE(IOUERR,195,ERR=98) TRIM(FILENAME) 
  195   FORMAT('*** Error reading  external data file ',A80, 
     &        /'*** (Unexpected end of file.)') 
      ELSE 
        BACKSPACE(IOUALT) 
        READ(IOUALT,196,ERR=97,END=97) LASTREC 
  196   FORMAT(A) 
        WRITE(IOUERR,197,ERR=97) TRIM(FILENAME), TRIM(LASTREC) 
  197   FORMAT('*** Error reading external data file ',A80, 
     &        /'*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
C 
  97  CLOSE(IOUALT) 
  98  INERR = INERR + 1 
  99  RETURN 
      END 
