      SUBROUTINE GETT2EVP(INPREC,INERR,RC) 
C 
C 
C  GETT2EVP reads user supplied Tier2 EVAPORATIVE phase-in fractions. 
C  Values read in will replace the default values located in the arry T2_EVAP. 
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
C 
C  23 Apr 01 AIR Task 11: Redirected error messages to IOUERR. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  15 Sep 00 @EPA-elg, New Subroutine fixes Bug 204. 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C    /IOUCOM/    IOUALT, IOUERR. IOUREP 
C    /TIER2/     T2_EVAP 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IOUCOM/   IOUALT, IOUERR, IOUREP 
C    /TIER2/    T2_EVAP 
C 
C  Common block array Subscripts: 
C 
C    T2_EVAP           T2_EVAP(MY,IV) 
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
C   DUMMY      R    Temporary array to store user supplied Tier2 EVAP 
C                   fractions prior to validation of input. 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'TIER2.I' 
C 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     SOURCE 
      INTEGER     I 
      INTEGER     IV 
      INTEGER     NAGE 
      INTEGER     LABTYPE 
      INTEGER     LDGIMP 
      INTEGER     NCHAR 
      INTEGER     SECTION 
C 
      REAL        DUMMY(12,5) 
C 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*40 
      CHARACTER   INPREC*(*) 
      CHARACTER   M6REC*150 
      CHARACTER   LASTREC*150 
      CHARACTER   STDLAB*40 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('T2 EVAP PHASE-IN',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Find out whether we are reading from an external data file. 
C 
      READ(INPREC(21:),*,IOSTAT=RC,END = 10,ERR=97) FILENAME 
C 
C       Open the external data file. 
C 
  10  IF ( RC.EQ.0 ) THEN 
         SOURCE = IOUALT 
         OPEN(SOURCE,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *        IOSTAT=RC,ERR=80) 
C 
C     Read the required header line from the file. 
C 
         CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
         IF (RC.EQ.0) THEN 
C 
C     Check for a valid labeled input record. 
C 
           CALL UCCOMP(M6REC(1:19),LABEL) 
           CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
           IF (LABTYPE.EQ.LDGIMP) THEN 
             READ(SOURCE,*,ERR=98,END=98,IOSTAT=RC) DUMMY 
             CLOSE(SOURCE) 
             WRITE(IOUREP,300) 
  300        FORMAT(/,2x, 
     &         'Reading User Supplied Tier2 EVAP phase-in fractions',/) 
C 
C     Check to verify that each EVAP fraction is between zero and one. 
C     If not, return to DSPTCH without updating T2_EVAP. 
C 
             DO IV = 1, 5 
                 DO NAGE = 1, 12 
                     IF (DUMMY(NAGE,IV).LT.0.0) THEN 
                       CALL QUITER(DUMMY(NAGE,IV),0,557,INERR) 
                       RETURN 
                     END IF 
                     IF (DUMMY(NAGE,IV).GT.1.0) THEN 
                       CALL QUITER(DUMMY(NAGE,IV),0,557,INERR) 
                       RETURN 
                     END IF 
                     T2_EVAP(NAGE,IV) = DUMMY(NAGE,IV) 
                 END DO 
             END DO 
C 
C 
C      Echo Print FILENAME for user supplied T2_BIN_FRAC 
C 
            WRITE(IOUREP,301) FILENAME 
  301       FORMAT('     Data read from file: ',A80,/ ) 
C 
C 
           ELSE   !  LABTYPE.NE.LDGIMP 
             WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
             RC = 2 
           END IF 
         ELSE 
           SELECT CASE (RC) 
           CASE (:-1) 
             WRITE(IOUERR,185) TRIM(FILENAME) 
  185        FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
             RC = 1 
           CASE (1:) 
             WRITE(IOUERR,186) TRIM(FILENAME) 
  186        FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
           END SELECT 
         END IF 
      ELSE ! RC not equal to zero 
C 
C     Error opening the external data file. 
C 
   80    WRITE(IOUERR,180) FILENAME 
  180    FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
   97    RC = MOD(RC,256) 
         IF (RC.EQ.71) THEN 
           WRITE(IOUERR,181) 
  181      FORMAT('*** (The file name or path could not be found.)') 
           DO I = LEN(FILENAME),1,-1 
             NCHAR = I 
             IF (FILENAME(I:I).NE.' ') EXIT 
           END DO 
           IF (NCHAR.GT.50) THEN 
             WRITE(IOUERR,182) 
  182        FORMAT('*** (The file name must not be longer than 50 ', 
     &        'characters.)') 
             INERR = INERR + 1 
           END IF 
         ELSE 
   98      IF (RC.NE.0) THEN 
              WRITE(IOUERR,183) 
  183         FORMAT('*** Error opening external file for User ', 
     &               'Supplied Tier2 EVAP phase-in fractions.') 
              IF (RC.LT.0) THEN 
                 WRITE(IOUERR,198) 
  198            FORMAT('*** (Unexpected end of file.)') 
              ELSE 
                 BACKSPACE(SOURCE) 
                 READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199            FORMAT(A) 
                 WRITE(IOUERR,200) LASTREC 
  200            FORMAT('*** The following record caused the error:', 
     &              /'*** ',A) 
              END IF 
           END IF 
         ENDIF 
      END IF 
C 
   99 CLOSE(IOUALT) 
      RETURN 
      END 
