      SUBROUTINE GETSOAK(INPREC,INERR,RC) 
C 
C    GETSOAK reads in the user supplied START SOAK DISTRIBUTION 
C    array. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  13 Oct 00 AIR Task 02: Bug Fix 273: Removed SOAKFLAG code to allow for 
C            multi-scenario SOAK DISTRIBUTION commands. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  31 Aug 99 @DynTel-ddj 2-699, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS2/ IMFLAG 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS5/ TAMFLG 
C    /SKDIST/ SOAKDIST 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C FILENAME  C    name of external data file holding soak numbers, for SOAKDIST 
C LABEL     C    dummy variable to hold first character of standard label; used 
C                to suppress compiler warning 
C LASTREC   C    holds string which causes read error 
C M6REC     C    holds mobile6 command string from external data file for checks 
C                against standard label 
C STDLAB    C    holds standard returned from CHKLAB/ INILAB, i.e., lablist 
C SOAKFAT   I    label number 
C EXTDATA   L    used to check label record for external data file 
C LABTYPE   I    equated to SECTION 
C SECTION   I    returned from CHKLAB/lablist; indicates which section, run level(=1), 
C                scenario level(=3), run or scenario level(=2), or header level(=0) 
C                the label should be in 
C SOURCE    I    FORTRAN unit number assigned to external data file (=IOUGEN) 
C I,J,K     I    indices for reading SOAKDIST 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'SKDIST.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     SOAKFAT 
C 
      LOGICAL     EXTDATA 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
      INTEGER     SECTION 
      INTEGER     SOURCE 
      CHARACTER   STDLAB*40 
      INTEGER     I,J,K 
C 
      RC = 0 
C 
C     Find out whether we are reading from an external data file. 
C 
      IF (INPREC(20:20).EQ.'@') THEN 
        EXTDATA = .FALSE. 
        SOURCE = IOUGEN 
        M6REC = INPREC 
      ELSE 
C 
C     Read the data from an external data file. 
C 
        EXTDATA = .TRUE. 
        SOURCE = IOUALT 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
        CALL CHKLAB('SOAK DISTRIBUTION',STDLAB,SOAKFAT,SECTION) 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C     Open the external data file. 
C 

ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), FILENAME)

        OPEN(SOURCE,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *       IOSTAT=RC,ERR=80) 
        IF (RC.NE.0) THEN 
C 
C     Error opening the external data file. 
C 
   80     WRITE(IOUERR,280) FILENAME 
  280     FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
          INERR = INERR + 1 
        ELSE 
          WRITE(IOUREP,200) FILENAME 
  200     FORMAT(/'* Reading start SOAK distribution from the ', 
     &          'following external',/'* data file: ',A) 
C 
C       Read the required header line from the file. 
C 
          CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
          IF (RC.NE.0) THEN 
C 
C     Error reading the header of the external data file. 
C 
            SELECT CASE (RC) 
              CASE (:-1) 
                WRITE(IOUERR,285) FILENAME 
  285           FORMAT(/'*** ERROR: no header record was found in ', 
     &          'the following external data file:',/'***',A) 
                RC = 1 
                INERR = INERR + 1 
              CASE (1:) 
                WRITE(IOUERR,286) FILENAME 
  286           FORMAT(/'*** ERROR reading the header of the ', 
     &                'following external data file:',/'*** ',A) 
              CASE DEFAULT 
            END SELECT 
          ELSE 
C 
C       Check for a valid labeled input record. 
C 
            CALL UCCOMP(M6REC(1:19),LABEL) 
            CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
            IF (LABTYPE.NE.SOAKFAT) THEN 
C 
C     Invalid label on the header record. 
C 
              WRITE(IOUERR,290) M6REC 
  290         FORMAT(/'*** ERROR: the following record is not ', 
     &        'allowed in the header of the external data file:', 
     &         /'*** ',A) 
              RC = 2 
              INERR = INERR + 1 
            ELSE 
C 
C     Read the data field from the input record. 
C     Find the first data line. 
C 
              CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
              IF (RC.NE.0) THEN 
   98           CALL QUITER(0.,0,126,INERR) 
                IF (RC.LT.0) THEN 
                  WRITE(IOUERR,298) 
  298             FORMAT('*** (Unexpected end of file.)') 
                ELSE 
                  BACKSPACE(SOURCE) 
                  READ(SOURCE,199,ERR=98,END=98) LASTREC 
  199             FORMAT(A) 
                  WRITE(IOUERR,299) LASTREC 
  299             FORMAT('*** The following record caused the error:', 
     &                  /'*** ',A) 
                ENDIF 
                RETURN 
              ELSE 
                BACKSPACE(SOURCE) 
C 
C      READ IN DATA. 
C 
                READ(SOURCE, *, END=98, IOSTAT=RC) 
     &              (((SOAKDIST(I,J,K), I=1,24), J=1,70), K=1,2) 
C 
                IF (EXTDATA) CLOSE(IOUALT) 
C 
              END IF 
            END IF 
          END IF 
        END IF 
      END IF 
      RETURN 
      END 
