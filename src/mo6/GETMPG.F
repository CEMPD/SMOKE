      SUBROUTINE GETMPG(INPREC,INERR,RC) 
C 
C  GETFEMPG reads the gasoline and diesel fuel economies,  
C  miles per gallon, from the user-supplied file. 
C 
C  Changes: (Last change first) 
C 
C  06 Aug 02 @EPA-bag bug 22: Renamed label in external file
C             from ROADFE to MPG ESTIMATES
C  13 Mar 02 AIR Task 27: Replaced MPG.I with RLCOM.I. 
C  06 Nov 01 AIR Task 27: New Subroutone 
C 
C  Called by: DSPTCH. 
C 
C  Calls CHKLAB, UCCOMP 
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
C 
C  ROADFE(99,28,2) - ROADFE (MY,IV,IMYU) 
C 
C  MY is Model Year Group, 1=1952-, 101=2050 
C  IV is the MOBILE6 Vehicle Class 
C  IM6U is the MOBILE6 or user-supplied toggle 
C 
C  Local variable Subscripts: 
C 
C    DUMMY       DUMMY( MY ) 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store 99 model years (1952- through 
C                   2050 of the gas and diesel fuel economies (mi/gal). 
C   I          I    Do Loop index 
C   L          I    Lines read in 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'RLCOM1.I' 
C 
      INTEGER, EXTERNAL :: JUNIT 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER     DATAFIELD*80 
      CHARACTER     FILENAME*80 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       IV 
      INTEGER       L 
      INTEGER       LABTYPE 
      INTEGER       LDGIMP 
      INTEGER       MY 
      INTEGER       PTR1 
      INTEGER       PTR2 
      INTEGER       SECTION 
      REAL          DUMMY(99) 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=50) M6DIR 
      CHARACTER(LEN=80) M6FILE 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Set MPG flag to false to use default M6 values 
C 
      MPGFLAG = .FALSE. 
C 
C    Look for the first data field (gas carbon ZML filename) following the delimiter 
C    in column 20. Initialize by setting PTR2 to 21 (the position of the blank 
C    following the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME) 
C 
C     Open the external data file for fuel economy. 
C 
ccs   All M6 input files live in SMK_M6PATH 
      CALL GETENV( "SMK_M6PATH", M6DIR ) 
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME 
       
      IOUALT = JUNIT() 
      OPEN(IOUALT,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading Gas and Diesel Fuel Economies'/ 
     &        '* from the external data file ',A) 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('MPG ESTIMATES',STDLAB,LDGIMP,SECTION)
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
C     Read the file containing the Gas and Diesel Fuel Economies (mi/gal) 
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
C  Read and validate the fuel economy levels. 
C 
        READ (IOUALT,*,IOSTAT=RC) IV,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
C  Verify the indexing input 
C 
        IF (IV.GT.MAXVEH.OR.IV.LT.1) GOTO 95 
C 
C  Store the fuel economy values in the second half of the ROADFE array. 
C  The first half contains the default MOBILE6 values. 
C 
        DO MY=1,99 
          ROADFE(MY,IV,2) = DUMMY(MY) 
        END DO 
C 
        L = L + 1  ! Update lines read in 
C 
      END DO 
C 
C     Make sure all fuel economy lines are read in since EXCEL does not terminate 
C     the file properly. Set MPGFLAG to true to use user-supplied values. 
C 
      IF(L.EQ.MAXVEH) THEN 
        RC = 0 
        MPGFLAG=.TRUE. 
        GOTO 99 
      ELSE 
        WRITE(IOUOUT,170,ERR=98) TRIM(FILENAME), L 
        WRITE(IOUERR,170,ERR=98) TRIM(FILENAME), L 
 170    FORMAT('*** ERROR:  invalid number of data lines in the'/ 
     &         '    external Fuel Economy data file ',A,'. The number'/ 
     &         '    of lines read in was ',I3,', out of 28. Please'/ 
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
 180  FORMAT('*** Error opening external Fuel Economy data file ',A) 
      RC = 2 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            '    the external Fuel Economy data file: ',A) 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            '    external Fuel Economy data file ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) FILENAME, TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in the'/ 
     &        '    header of the external Fuel Economy data file ',A/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   95 WRITE(IOUERR,195) FILENAME, TRIM(M6REC) 
  195 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    indexes or data in the external Fuel Economy'/ 
     &        '    data file ',A 
     &       /1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
