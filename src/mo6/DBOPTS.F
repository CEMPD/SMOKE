      SUBROUTINE DBOPTS(INPREC,INERR,RC) 
C 
C  DBOPTS reads database selection records (SELECT EMISSIONS, FACILITIES, 
C  POLLUTANTS, or VEHICLES) from an external file. 
C  MODIFIED TO ALSO PROCESS DAILY OUTPUT AND WITH FIELDNAMES LABELS 
C 
C  Called by PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection 
C 19 Oct 00 @EPA-djb Task W2, Added AGGREGATED OUTPUT command. 
C 02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C 13 Sep 00 @EPA-elg bug228 changed error handling conditions. 
C 05 MAY 00 @EPA TASK E22: ADD DAILY DATABASE OUTPUT OPTIONS 
C 11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C           so that Lahey Standard for keyboard/monitor I/O is followed. 
C 10 Apr 00 AIR Task 02: Changed DATABASE POLLUTANTS to POLLUTANTS 
C 10 Feb 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELVEH 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /IOUCOM/ IOUERR 
C    DATABASE MODULE: 
C      DBDAILY AND DBHEAD 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C 
C  Notes: 
C 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER, EXTERNAL :: JUNIT
C
      INTEGER,INTENT(IN OUT)           ::  INERR 
      CHARACTER(LEN=*),INTENT(IN)      ::  INPREC 
      INTEGER,INTENT(OUT)              ::  RC 
C 
      CHARACTER (LEN=150)  ::  BUFFER1 
      CHARACTER (LEN=150)  ::  BUFFER2 
      INTEGER              ::  DBOPTIONS 
      CHARACTER (LEN=80)   ::  FILENAME 
      CHARACTER (LEN=19)   ::  LABEL 
      INTEGER              ::  LABTYPE 
      CHARACTER (LEN=150)  ::  M6REC 
      INTEGER              ::  SECTION 
      CHARACTER (LEN=40)   ::  STDLAB 
      INTEGER, SAVE        ::  IYES=1 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=50) M6DIR
      CHARACTER(LEN=80) M6FILE
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Find the token number of the label that must appear in the input 
C     label field of the first data record in the header of the external 
C     database options selection file. 
C 
      LABEL = 'DBASE OPT' 
      CALL CHKLAB(LABEL,STDLAB,DBOPTIONS,SECTION) 
C 
C     Open the external data file. 
C 

ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)

ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME

      IOUALT = JUNIT()

      OPEN(IOUALT,FILE=M6FILE,ACTION='READ',IOSTAT=RC,ERR=85) 
      WRITE(IOUOUT,120) FILENAME 
  120 FORMAT(/'* Reading database selection criteria from the ', 
     &        'following external',/'* data file: ',A) 
C 
C     Read the required header line from the file. 
C 
      CALL NXTREC(IOUALT,IOUOUT,M6REC,INERR,RC) 
      IF (RC.NE.0) THEN 
        GOTO 95 
      ENDIF 
C 
C     Check for a valid labeled input record. 
C 
      CALL UCCOMP(M6REC(1:19),LABEL) 
      CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
C     ! The next line is present only to suppress compiler-generated 
C     ! warning messages. It has no functional importance. 
      IF (LABTYPE /= DBOPTIONS) GOTO 90 
C 
      DO  ! Forever 
C 
         CALL NXTREC(IOUALT, IOUOUT, M6REC, INERR, RC) 
         IF (RC < 0) EXIT 
         IF (RC > 0) GOTO 95 
C 
         CALL UCCOMP(M6REC(1:19), LABEL) 
         CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
C 
         SELECT CASE (STDLAB) 
            CASE ('DATABASE EMISSIONS', 
     &            'DBASE EMISSIONS') 
               CALL DBSELE(M6REC,INERR,RC) 
            CASE ('DATABASE FACILITIES', 
     &            'DBASE FACILITIES') 
               CALL DBSELF(M6REC,INERR,RC) 
            CASE ('POLLUTANTS', 
     &            'POL') 
               CALL DBSELP(M6REC,INERR,RC) 
            CASE ('DATABASE VEHICLES', 
     &            'DBASE VEHICLES') 
               CALL DBSELV(M6REC,INERR,RC) 
            CASE ('DATABASE MODEL YEARS', 
     &            'DBASE MODEL YEARS') 
               DB_MDLYR=IYES 
               CALL DBSELA(M6REC,INERR,RC) 
            CASE ('DATABASE VEHICLE AGE', 
     &            'DBASE VEHICLE AGE') 
               CALL DBSELA(M6REC,INERR,RC) 
            CASE ('DATABASE HOURS', 
     &            'DBASE HOURS') 
               CALL DBSELH(M6REC,INERR,RC) 
            CASE ('DAILY OUTPUT') 
               DBDAILY = .TRUE. 
               DBAGGR  = .FALSE. 
            CASE ('AGGREGATED OUTPUT') 
               DBDAILY = .FALSE. 
               DBAGGR  = .TRUE. 
            CASE ('WITH FIELDNAMES') 
               DBHEAD = .TRUE. 
            CASE DEFAULT 
               WRITE(IOUOUT,130) M6REC 
               WRITE(IOUERR,130) M6REC 
  130          FORMAT (' *** The following record is not recognized', 
     &                 ' in a database selection file:',/' *** ',A) 
               RC = 1 
               GOTO 99 
         END SELECT 
C 
      END DO 
 
      IF (RC < 0) THEN 
         RC = 0 
      ENDIF 
C 
      CLOSE(IOUALT) 
      GOTO 99 
C 
C     Error opening the external data file. 
C 
   85 WRITE(IOUOUT,185) FILENAME 
      WRITE(IOUERR,185) FILENAME 
  185 FORMAT(/' *** ERROR: unable to open the following external', 
     &        ' data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUOUT,190) M6REC 
      WRITE(IOUERR,190) M6REC 
  190 FORMAT(/' *** ERROR: the following record is not allowed in', 
     &        ' the header of the external data file:', 
     &      /' *** ',A) 
      RC = 2   
      INERR = INERR + 1 
      GOTO 99 
C 
   95 IF (RC > 0) THEN 
         WRITE(BUFFER1,195) RC 
  195    FORMAT(/'*** ERROR ',I9,' reading the following record', 
     &        ' from the database options file:') 
         CALL LJCOMP(BUFFER1,BUFFER2) 
         WRITE(IOUOUT,"(' ',A,/' *** ',A)") BUFFER2,M6REC 
         WRITE(IOUERR,"(' ',A,/' *** ',A)") BUFFER2,M6REC 
      ELSE IF (RC < 0) THEN 
         WRITE(IOUOUT,196) FILENAME 
         WRITE(IOUERR,196) FILENAME 
  196    FORMAT(' *** ERROR  No header was found in the following', 
     &          ' external database options file:', 
     &         /' *** ',A) 
      ENDIF 
      INERR = INERR + 1 
C 
   99 RETURN 
      END 
