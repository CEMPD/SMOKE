      SUBROUTINE GETPMZML(FILENAME,INERR,RC) 
C 
C  GETPMSCF reads the base BER sulfur levels, the PM DR Ages, and the 
C  PM Zero Mile Levels (g/mi) from the user-supplied file. 
C 
C  Changes: (Last change first) 
C 
C  28 Aug 01 AIR Task 23c: Implement PMDRAGE. 
C  24 Aug 01 AIR Task 23b: Implement EPA-based PM ZMLs/DRs  
C  07 Jun 01 AIR Task 17: Added input of base sulfur levels and CA diesel aromatic 
C                         adjustment flag. 
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
C 
C  PMZML(71,2,28) - PMZML (MY,ICAT,IV) 
C  PMDRAGE(28)    - PMDRAGE(IV) 
C 
C  MY is Model Year Group, 1=1950-, 71=2020+ 
C 
C  ICAT is for catalyst equipment, as follows: 
C    1 = no catalyst, diesel 
C    2 = catalyst 
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
C                   2020+ of the EMFAC2000v2.02-based Zero Mile Levels (g/mi). 
C   I          I    Do Loop index 
C   L          I    Lines read in 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART4.I' 
C 
      CHARACTER*(*) FILENAME 
      INTEGER       INERR 
      INTEGER       RC 
C 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       ICAT 
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
C     Open the external data file for PM ZML. 
C 
      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading PM Zero Mile Levels '/ 
     &        '* from the external data file ',A) 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('PMZML',STDLAB,LDGIMP,SECTION) 
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
C     Find, read and verify the base sulfur levels for the 
C     PM ZMLs/DRs. The input order is: Gas CY1999-, Gas CY2000+, 
C     and Diesel. 
C 
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC.LT.0) GOTO 91 
      BACKSPACE(IOUALT) 
C 
      READ (IOUALT,*,IOSTAT=RC) BSULFLVL 
C 
      DO I=1,3 
      IF(BSULFLVL(I).LT.   1. .OR. 
     *   BSULFLVL(I).GT.5000.) GOTO 91 
      END DO 
C 
C     Find, read and verify the age when PMDR1 ends and PMDR2 begins. 
C     Order is MOBILE6 IV. 
C 
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC.LT.0) GOTO 91 
      BACKSPACE(IOUALT) 
C 
      READ (IOUALT,*,IOSTAT=RC) PMDRAGE 
C 
      DO IV=1,MAXVEH 
      IF(PMDRAGE(IV).LT. 0 .OR. 
     *   PMDRAGE(IV).GT.25) GOTO 93 
      END DO 
C       
C     Read the file containing the PM Zero Mile Levels (g/mi) 
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
C       Read and validate the zero mile levels. 
C 
        READ (IOUALT,*,IOSTAT=RC) IV,ICAT,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
        IF (IV.GT.MAXVEH.OR.IV.LT.1 .OR. 
     *      ICAT.GT.2.OR.ICAT.LT.1) GOTO 92 
C 
        DO MY=1,71 
          PMZML(MY,ICAT,IV) = DUMMY(MY) 
        END DO 
C 
        L = L + 1  ! Update lines read in 
C 
      END DO 
C 
C     Make sure all ZMLs are read in since EXCEL does not terminate 
C     the file properly. 
C 
      IF(L.EQ.56) THEN 
        RC = 0 
        GOTO 99 
      ELSE 
        WRITE(IOUOUT,170,ERR=98) TRIM(FILENAME), L 
        WRITE(IOUERR,170,ERR=98) TRIM(FILENAME), L 
 170    FORMAT('*** ERROR:  invalid number of data lines in the'/ 
     &         '    external PM ZML data file ',A,'. The number'/ 
     &         '    of lines read in was ',I3,', out of 56. Please'/ 
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
 180  FORMAT('*** Error opening external PM ZML data file ',A) 
      RC = 2 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            '    the external PM ZML data file: ',A) 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            '    external PM ZML data file ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) FILENAME, TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in'/ 
     &        '    the header of the external PM ZML data file ',A/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Missing or invalid sulfur data. 
C 
   91 WRITE(IOUERR,191) FILENAME, TRIM(M6REC) 
  191 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    sulfur data in the external PM ZML data file ',A/ 
     &       1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   92 WRITE(IOUERR,192) FILENAME, TRIM(M6REC) 
  192 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    indexes or data in the external PM ZML data file ',A 
     &       /1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   93 WRITE(IOUERR,192) FILENAME, TRIM(M6REC) 
  193 FORMAT(/'*** ERROR: the following record contains invalid PM'/ 
     &        '    DR age indexes in the external PM ZML data file ',A 
     &       /1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
