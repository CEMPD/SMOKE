      SUBROUTINE GETSCF(INPREC,INERR,RC) 
C 
C  GETSCF reads in, validates and stores SCF for speed 
C  distribution at different hours from an external data file. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C 11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C           so that Lahey Standard for keyboard/monitor I/O is followed. 
C 28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C           module subroutine CHKLAB was removed from that module 
C           and converted to an ordinary external procedure. 
C 11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C           calling sequence for that subroutine was changed. 
C  1 March 1999 @DynTel-HQ 2-698 Read in speed correction coefficient 
C                          from user input 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS1/ SCFFLG 
C    /SPEED9/ SCFF SCFAC 
C 
C  Local variable / array dictionary: 
C 
C  Name    Type  Description 
C  ------- ----  --------------------------------------------------- 
C  SECTION  I    Dummy argument for CHKLAB. 
C  STDLAB   C    Dummy argument for CHKLAB. 
C 
C  Notes: 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'SPEED9.I' 
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     SCFFAT 
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
      INTEGER  I, J, K 
C 
      RC = 0 
      SCFUSER=2 
      IF (SCFFLG.EQ.1) THEN 
        SCFFLG = 2 
      ELSE 
        GOTO 96 
      ENDIF 
C 
C     Find out whether we are reading from an external data file. 
C 
      IF (INPREC(20:20).EQ.'@') THEN 
        EXTDATA = .FALSE. 
        SOURCE = IOUGEN 
        M6REC = INPREC 
      ELSE 
C 
C       Read the data from an external data file. 
C 
        EXTDATA = .TRUE. 
        SOURCE = IOUALT 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
        CALL CHKLAB('SCF BY FWAY AND AC',STDLAB,SCFFAT,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), FILENAME)

        OPEN(SOURCE,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     &      IOSTAT=RC,ERR=80) 
        WRITE(IOUREP,200) FILENAME 
  200   FORMAT(/'* Reading non-default SCF BY FWAY AND AC from the ', 
     &          'following external',/'* data file: ',A) 
C 
C       Read the required header line from the file. 
C 
        CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) THEN 
          GOTO 85 
        ENDIF 
C 
C       Check for a valid labeled input record. 
C 
        CALL UCCOMP(M6REC(1:19),LABEL) 
        CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
        IF (LABTYPE.NE.SCFFAT) THEN 
          GOTO 90 
        ENDIF 
C 
      ENDIF 
C 
C     Read the data field from the input record. 
C     Find the first data line. 
C 
      CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
      IF (RC.NE.0) GOTO 98 
      BACKSPACE(SOURCE) 
C 
C 
C     WE READ IN DATA HERE. 
C   I ---- 1 to 24 Hour 
C 
  100 FORMAT(12F5.2) 
      DO 34 K=1, 14 
      READ(SOURCE, 100, END=98, IOSTAT=RC)  
     * ((SCFF(J,K,I), J=1,3),I=1,4) 
   34 CONTINUE 
      DO 36 K=1, 14 
      READ(SOURCE, 100, END=98, IOSTAT=RC)  
     * ((SCFAC(J,K,I), J=1,3),I=1,4) 
   36 CONTINUE 
 
C 
C 
      IF (EXTDATA) CLOSE(IOUALT) 
C 
   40 GOTO 99 
C 
C     Error opening the external data file. 
C 
   80 WRITE(IOUERR,280) FILENAME 
  280 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,285) FILENAME 
  285     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
          RC = 1 
          INERR = INERR + 1 
        CASE (1:) 
          WRITE(IOUERR,286) FILENAME 
  286     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 99 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,290) M6REC 
  290 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 2   
      INERR = INERR + 1 
      GOTO 99 
C 
C     Error: alternate speed VMT have been 
C            read in. 
C     
   96 WRITE(IOUERR,296) 
  296 FORMAT(/'*** ERROR: Alternate SFCF and SFCAC', 
     &        'have already been read in this run.', 
     &       /'*** Only one set of SFCF and AFCAC', 
     &        'can be entered in one run.'/) 
CCC  &      /'*** The factors entered here will replace the ', 
CCC  &       'rates that were entered previously.'/) 
      RC = 3   
      INERR = INERR + 1 
      GOTO 99 
C 
   98 CALL QUITER(0.,0,126,INERR) 
      IF (RC.LT.0) THEN 
        WRITE(IOUERR,298) 
  298   FORMAT('*** (Unexpected end of file.)') 
      ELSE 
        BACKSPACE(SOURCE) 
        READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199   FORMAT(A) 
        WRITE(IOUERR,299) LASTREC 
  299   FORMAT('*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
C 
   99 RETURN 
      END 
