      SUBROUTINE GETXIM(INPREC,INERR,RC) 
C 
C  GETXIM processes the I/M EXTERNAL FILE labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,GETCMO,GETCOV,GETCMO,GETCPC,GETCIM,GETCPL,GETCUT,GETEFR, 
C        GETFAL,GETIMV,GETIMY,GETIPR,GETGRA,GETLIM,GETLOS,GETNVT,GETRCC, 
C        GETRNC,GETRSO,GETRTP,GETRTV,GETSGC,GETTEF,GETTTC,GETWAV,GETXMP, 
C        GETWRR,GETFBEN,LJCOMP,NXTTOK,QUITER,UCOMP 
C 
C  Changes: (Last change first) 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C 15 Sep 00 AIR Task 08: Added EVAP I/M FRAC BEN. 
C 13 Oct 99 @Dyntel-JWR Amended PROCLABEL section to read M6REC. Also changed 
C           the TTC label to relflect change in the default. Commented out the 
C           non-periodic COIM and RSD labels. 
C 
C 28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C           module subroutine CHKLAB was removed from that module 
C           and converted to an ordinary external procedure. 
C 12 May 99 @DynTel-MLA 2-694  Added cases for OBD I/M and modified the 
C           call to CHKLAB because the calling sequence for that subroutine 
C           was changed. 
C  3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUGEN,IOUREP 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /CHKLST/ IMCHECK 
C    /IMPAR6/ IFILET 
C    /IMPAR7/ TNAME 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  BUFFER1    C   Buffer used to compose output messages. 
C  BUFFER2    C   Buffer used to compose output messages. 
C  FILENAME   C   Name of external I/M description data file. 
C  LABEL      C   The label field from a labeled input record. 
C  LABTYPE    I   A numeric code for the label type. 
C  M6REC      C   A record from the external I/M description data file. 
C  SECTION    I   Dummy argument for CHKLAB. 
C  STDLAB     C   Dummy argument for CHKLAB. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER, EXTERNAL :: JUNIT   
C
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*150  BUFFER1 
      CHARACTER*150  BUFFER2 
      CHARACTER*80   FILENAME 
      CHARACTER*19   LABEL 
      INTEGER        LABTYPE 
      CHARACTER*150  M6REC 
      INTEGER        SECTION 
      CHARACTER      STDLAB*40 
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=50) M6DIR
      CHARACTER(LEN=80) M6FILE
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Get the external data file name from the input 
C     record and open the file. 
C 

ccs      CALL UCCOMP(INPREC(22:), FILENAME) 
      CALL LJCOMP(INPREC(22:), FILENAME)

ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME
      
      IOUXIM = JUNIT()
      OPEN(IOUXIM,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,110) FILENAME 
  110 FORMAT(/'* Reading I/M program description records from the', 
     &        ' following external',/'* data file: ',A) 
C 
C     Read and process records from the external file. 
C 
   20 CALL NXTREC(IOUXIM,IOUREP,M6REC,INERR,RC) 
      IF (RC.LT.0) THEN 
        GOTO 95 
      ELSE IF (RC.GT.0) THEN 
        GOTO 85 
      ENDIF 
C 
C     Check for a valid labeled input record. If the label type is 
C     unrecognized, print an error message and continue with the next 
C     input record. 
C 
      CALL UCCOMP(M6REC(1:19),LABEL) 
   25 CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
C     The next assignment statement is present solely to suppress 
C     compiler generated warnings. It has no functional importance. 
      RC = SECTION 
C 
      IF (LABTYPE.EQ.0) THEN 
        WRITE(IOUERR,130) M6REC 
  130   FORMAT('*** ERROR: The following record does not contain', 
     &         ' a valid input label:',/A) 
        INERR = INERR + 1 
        GOTO 20 
      ENDIF 
C 
C     Process the input record. (NOTE: the following list 
C     must be maintained in parallel with the list in DSPTCH.) 
C 
      RC = 0 
C 
      PROCLAB: SELECT CASE (STDLAB) 
C 
        CASE ('I/M PROGRAM', 
     &        'I/M PGM') 
          CALL GETIPR(M6REC,INERR,RC) 
        CASE ('I/M MODEL YEARS') 
          CALL GETIMY(M6REC,INERR,RC) 
        CASE ('I/M EFFECTIVENESS') 
          CALL GETEFR(M6REC,INERR,RC) 
        CASE ('I/M CREDIT FILE') 
C         *** Save file name in common IMPAR7. 

ccs          CALL UCCOMP(M6REC(22:),TNAME) 
          CALL LJCOMP(M6REC(22:),TNAME)   

          IFILET = 2 
        CASE ('I/M VEHICLES') 
          CALL GETIMV(M6REC,INERR,RC) 
        CASE ('I/M STRINGENCY') 
          CALL GETSGC(M6REC,INERR,RC) 
        CASE ('I/M COMPLIANCE') 
          CALL GETCPL(M6REC,INERR,RC) 
        CASE ('I/M WAIVER RATES') 
          CALL GETWAV(M6REC,INERR,RC) 
        CASE ('I/M WAIVER REDUCTION') 
          CALL GETWRR(M6REC,INERR,RC) 
        CASE ('I/M CUTPOINTS') 
          CALL GETCUT(M6REC,INERR,RC) 
        CASE ('I/M EXEMPTION AGE') 
          CALL GETXMP(M6REC,INERR,RC) 
        CASE ('I/M GRACE PERIOD') 
          CALL GETGRA(M6REC,INERR,RC) 
        CASE ('NO I/M TTC CREDITS') 
          CALL GETTTC(M6REC,INERR,RC) 
        CASE ('I/M TTC EFFECTIVENESS') 
          CALL GETTEF(M6REC,INERR,RC) 
        CASE ('EVAP I/M FRAC BEN') 
          CALL GETFBEN(INPREC,INERR,RC) 
        CASE ('I/M DESCRIPTION FILE') 
          WRITE(IOUERR,160) FILENAME 
  160     FORMAT('*** ERROR: I/M DESCRIPTION FILE labeled input', 
     &           ' records can not be used', 
     &          /' *** in an external I/M descripiton data file.', 
     &           ' The record was found in', 
     &          /'*** file ',A) 
          INERR = INERR + 1 
        CASE DEFAULT       
C          *** The label is not an I/M record type. 
          WRITE(IOUERR,170) M6REC 
  170     FORMAT('*** ERROR: The following record can not be', 
     &           ' used in an I/M descripiton data file:', 
     &          /'*** ',A) 
          INERR = INERR + 1 
C 
      END SELECT PROCLAB 
C 
      IF (RC.NE.0) THEN 
        WRITE(BUFFER1,190,IOSTAT=RC,ERR=95) RC,LABTYPE 
  190   FORMAT('*** ERROR ',I11,' processing the following input', 
     &         ' record (LABTYPE=',I3,') from') 
        CALL LJCOMP(BUFFER1,BUFFER2) 
        WRITE (IOUERR,194,IOSTAT=RC,ERR=95) BUFFER2,FILENAME,M6REC 
  194   FORMAT(A,/'*** external file ',A,/'*** ',A/) 
      ENDIF 
C 
C     Get the next records from the input file. 
C 
      GOTO 20 
C 
C     Error opening the external data file. 
C 
   80 WRITE(IOUERR,180) FILENAME 
  180 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     Error reading the external data file. 
C 
   85 CLOSE(IOUXIM) 
      WRITE(BUFFER1,185) RC,FILENAME 
  185 FORMAT('*** ERROR ',I9,' reading I/M description data', 
     &       ' file',A80) 
      CALL LJCOMP(BUFFER1,BUFFER2) 
      WRITE (IOUERR,186) BUFFER2 
  186 FORMAT(/A) 
      INERR = INERR + 1 
      GOTO 99 
C 
   95 CLOSE(IOUXIM) 
      IF (RC.LT.0) THEN 
        RC = 0 
      ENDIF 
C 
   99 RETURN 
      END 
