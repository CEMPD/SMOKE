      SUBROUTINE GETPMDDR1(FILENAME,INERR,RC) 
C 
C  GETPMDDR1 reads the first set of PM Deterioration Rates (g/10,000 miles) 
C  from the user-supplied file. 
C 
C  Changes: (Last change first) 
C 
C  28 Aug 01 Air Task 23c: Implemented DR1 and DR2 rates. 
C  24 Aug 01 AIR Task 23b: Implement EPA-based PM ZMLs/DRs  
C  25 May 01 AIR Task 17: Changed from fixed to user-supplied input file. 
C  23 May 01 AIR Task 17: Changed from 1975 to 1950 base and only 2 modes. 
C  09 Apr 01 AIR Task 11: New Subroutone 
C 
C  Called by: GETPMEF. 
C 
C  Calls CHKLAB 
C 
C  Input on call: 
C 
C    argument list: INERR,FILENAME 
C 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C 
C  Common block array Subscripts: 
c 
C  PMDDR1(71,28) - PMDDR1(MY,IV) 
C 
C  MY is Model Year Group, 1=1950-, 71=2020+ 
C 
C  IV is the MOBILE6 Vehicle Class 
C 
C  Local variable Subscripts: 
C 
C    DUMMY       DUMMY( MY ) 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store 71 model years (1950- through 
C                   2020+ of the EPA-based First Deterioration Rates 
C                   in g/10000 miles. 
C                   Equation Coefficients. 
C   I          I    Do Loop index 
C   L          I    Lines read in 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART4.I' 
C 
      INTEGER, EXTERNAL :: JUNIT 
C 
      CHARACTER*(*) FILENAME 
      INTEGER       INERR 
      INTEGER       RC 
C 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       IV 
      INTEGER       L 
      INTEGER       LABTYPE 
      INTEGER       LDGIMP 
      INTEGER       MY 
      INTEGER       SECTION 
      REAL          DUMMY(71) 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Open the external data file for PM DR. 
C       
      IOUALT = JUNIT() 
      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading the First PM Deterioration Rates '/ 
     &        '* from the external data file ',A) 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('PMDDR1',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
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
C     Read the file containing the first PM deterioration rates. 
C     Read data a line a time, skipping any blank lines or comments. 
C 
      L = 0  ! Initialize lines read 
C 
      DO I = 1,32767 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read and validate the deterioration rates. Data must be 
C       for diesel vehicles 
C 
        READ (IOUALT,*,IOSTAT=RC) IV,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
        IF(IV.GT.MAXVEH.OR.IV.LT.1) GOTO 92 
        IF(VVDSL(IV).NE.1) GOTO 92 
C 
        DO MY=1,71 
          PMDDR1(MY,IV) = DUMMY(MY) 
        END DO 
C 
        L = L + 1  ! Update lines read in 
C 
      END DO 
C 
C     Make sure all DRs are read in since EXCEL does not terminate 
C     the file properly. 
C 
      IF(L.EQ.13) THEN 
        RC = 0 
        GOTO 99 
      ELSE 
        WRITE(IOUOUT,170,ERR=98) TRIM(FILENAME), L 
        WRITE(IOUERR,170,ERR=98) TRIM(FILENAME), L 
 170    FORMAT('*** ERROR:  invalid number of data lines in the'/ 
     &         '    external PM DR1 data file ',A,'. The number'/ 
     &         '    of lines read in was ',I3,', out of 13. Please'/ 
     &         '    correct the file, especially checking for the'/ 
     &         '    proper termination of the last line.') 
        RC = 2 
        GOTO 98 
      END IF 
C 
C     Close the datafile and write the appropriate warning message. 
C 
  80  WRITE(IOUOUT,180,ERR=98) FILENAME 
      WRITE(IOUERR,180,ERR=98) FILENAME 
 180  FORMAT('*** Error opening external PM DR1 data file ',A) 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            'the external PM DR1 data file ',A) 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            'external PM DR1 data file: ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) FILENAME, TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in'/ 
     &        'the header of the external PM DR1 data file ',A/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   92 WRITE(IOUERR,192) FILENAME, TRIM(M6REC) 
  192 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        'indexes or data in the PM DR1 external data file ',A/ 
     &        1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
