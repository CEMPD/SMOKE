      SUBROUTINE GETREG(INPREC,INERR,RC) 
C 
C  GETREG checks for and reads in and screens user supplied registration 
C  distributions. 
C 
C  Called by DSPTCH 
C 
C  Calls CHKLAB, UCCOMP, NXTREC, QUITER. 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  25 May 2000 @EPA-elg Inputs are now Free Format rather than Fixed Format. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C   2 Mar 99 @DynTel-MLA 2-659  Cleaned up comments and error messages. 
C  21 OCT 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: INERR,INPREC,RC 
C    common blocks: 
C    /FLAGS2/ MYMRFG 
C    /IOUCOM/ IOUGEN,IOUALT 
C    /MAXIMA/ MAXYRS 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /MYRSAV/ AMAR,JULMYR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CHKMYR    R   myr mix check: sum of the myr fractions = 1.0 +/- 1.E-6 
C                (a normalized myr distribution has been entered for each IV) 
C  EXTDATA   L   Is an external data file open? 
C  FILENAME  C   Name of the external data file. 
C  I         I   Loop control variable. 
C  IV        I   Vehicle type. 
C  JDX       I   Model year. 
C  LABEL     C   The label field from a labeled input record. 
C  LABTYPE   I   Integer code for the label type. 
C  LASTREC   C   Previous record read from the input. 
C  M6REC     C   Record from the external data file. 
C  REGDIST   I   Numeric code for the REGISTRATION DISTRIB input label. 
C  SECTION   I   Dummy argument for CHKLAB. 
C  SOURCE    I   Logical unit number for the input file. 
C  STDLAB    C   Dummy argument for CHKLAB. 
C 
C  Notes: 
C 
C   2 Mar 99 @DynTel-MLA 2-659  Cleaned up error messages. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. Altered code to read 
C            Mobile6 transition input formats. 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C 
C  27 Mar 98 @DynTel-MLA 2-659  Updated for Mobile6 (the transitional 
C            labeled inputs version). This subroutine now sets MYMRFG. 
C            This subroutine was derived from the obsolete GETMRN, which 
C            was divided into the new GETMAR and GETREG subroutines. The 
C            consistency checks in GETMRN were moved to the new CHKMRN 
C            subroutine. 
C 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'NAMES1.I' 
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
      INTEGER    IVC,I,J 
      REAL       CHKMYR, DUM1(25) 
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C 
C     MYMRFG should be set to 2 if the user has entered only new mileage 
C     accumulation rates, to 3 if the user has entered only registration 
C     distributions, and to 4 if the user has entered both new mileage 
C     accumulation rates and registration distributions. 
C 
      IF (MYMRFG.EQ.3.OR.MYMRFG.EQ.4) THEN 
        GOTO 95 
      ELSE IF (MYMRFG.EQ.2) THEN 
        MYMRFG = 4 
      ELSE 
        MYMRFG = 3 
      ENDIF 
 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('REGIST DISTRIB',STDLAB,LDGIMP,SECTION) 
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
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME
      
      IOUALT = JUNIT()

      OPEN(IOUALT,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) TRIM(FILENAME) 
  120 FORMAT(/'* Reading Registration Distributions ', 
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
C  Get vehicle registration distributions. 
C 
C  Combined diesel/gas sales vehicle class groupings are: 
C 
C         IVC 
C        ----- 
C     1 : LDV 
C     2 : LDT1 
C     3 : LDT2 
C     4 : LDT3 
C     5 : LDT4 
C     6 : HDV2B 
C     7 : HDV3 
C     8 : HDV4 
C     9 : HDV5 
C    10 : HDV6 
C    11 : HDV7 
C    12 : HDV8A 
C    13 : HDV8B 
C    14 : HDBS (school bus) 
C    15 : HDBT (transit bus) 
C    16 : MC 
C 
C 
C 
      DO 31 I = 1,500 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read the implementation schedule for the current vehicle type. 
C 
C 
        READ (IOUALT,*,IOSTAT=RC) 
     &          IVC, (DUM1(J),J=1,25) 
C 
        IF(IVC.GE.1 .AND. IVC.LE.16) THEN 
            DO J = 1,25 
              JULMYR(J,IVC) = DUM1(J) 
            END DO 
 
        ELSE 
          WRITE(IOUERR,*) '*** ERROR: Illegal input for Mileage ', 
     *                    'Accumulation Rates', IVC 
          INERR = INERR + 1 
        ENDIF 
 
        IF(RC.LT.0) EXIT 
C 
        CHKMYR=0.0 
        DO 30 J=1,25 
          CHKMYR=CHKMYR + JULMYR(J,IVC) 
   30   CONTINUE 
C 
        IF(CHKMYR-.00001.GT.1..OR.CHKMYR+.00001.LT.1.) 
     *    CALL QUITER(CHKMYR,0,49,INERR) 
C 
   31 END DO 
C 
C 
      DO IVC=1,16 
        DO J=1,MAXYRS 
          IF(JULMYR(J,IVC).LT.0.0) 
     &        CALL QUITER(JULMYR(J,IVC),0,20,INERR) 
        END DO 
      END DO 
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
      GOTO 98 
C 
C 
   95 WRITE(IOUERR,195) 
  195 FORMAT(/'*** ERROR: Alternate registration distributions ', 
     &        'have already been read in this run.', 
     &       /'*** Only one set of registration distributions ', 
     &        'can be entered in each run.'/) 
      RC = 3   
C 
  98  INERR = INERR + 1 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
