      SUBROUTINE GETSDIST(INPREC,INERR,RC) 
C 
C    GETSDIST reads in the user supplied START DISTRIBUTION by hour 
C    array. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, NXTREC, UCCOMP, QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  29 Jun 00 @EPA-elg Initialize DISTSUM array to 0.0. 
C  23 May 00 @EPA-djb created to read engine start (by hour) distribution. 
C 
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
C    /STDIST/ STARTDIST,SDIST 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DISTUSM   R   Sum of hourly engine start distribution by weekday/weekend. 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'STDIST.I' 
C 
      INTEGER, EXTERNAL :: JUNIT   
C 
      CHARACTER, INTENT(IN) :: INPREC*(*) 
      INTEGER,INTENT(INOUT) :: INERR 
      INTEGER,  INTENT(OUT) ::  RC 
C 
      REAL, DIMENSION(2) :: DISTSUM 
C 
      INTEGER     I,J 
      INTEGER     LABTYPE 
      INTEGER     LABB 
      INTEGER     SECTION 
      INTEGER     SOURCE 
C 
      LOGICAL     EXTDATA 
C 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
      CHARACTER   STDLAB*40 
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
      RC = 0 
      SDIST = SDIST + 1 
      IF (SDIST.GT.1) THEN 
C 
C     Error: start distribution has already been read in. 
C     
        WRITE(IOUERR,296) 
  296 FORMAT(/'*** ERROR: Hourly Start Distribution', 
     &        'has already been read for this run.', 
     &       /'*** Only one Starts/hour Distribution ', 
     &        'can be entered per run.',/, 
     &        '*** The distribution entered here will be ignored.  ') 
 
        RC = 3 
        INERR = INERR + 1 
      ELSE 
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
          SOURCE = JUNIT() 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
          CALL CHKLAB('START DIST',STDLAB,LABB,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
          LABEL = STDLAB(1:1) 
          LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs          CALL UCCOMP(INPREC(21:), FILENAME) 
          CALL LJCOMP(INPREC(21:), FILENAME)

ccs       All M6 input files live in SMK_M6PATH
          CALL GETENV( "SMK_M6PATH", M6DIR )
          M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME

          OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *         IOSTAT=RC,ERR=80) 
          IF (RC.NE.0) THEN 
C 
C     Error opening the external data file. 
C 
   80       WRITE(IOUERR,280) FILENAME 
  280       FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
            INERR = INERR + 1 
          ELSE 
            WRITE(IOUREP,200) FILENAME 
  200       FORMAT(/'* Reading hourly start distribution from the', 
     &              ' following external',/'* data file: ',A) 
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
  285             FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
                  RC = 1 
                  INERR = INERR + 1 
                CASE (1:) 
                  WRITE(IOUERR,286) FILENAME 
  286             FORMAT(/'*** ERROR reading the header of the ', 
     &                  'following external data file:',/'*** ',A) 
                CASE DEFAULT 

              END SELECT 
            ELSE 
C 
C       Check for a valid labeled input record. 
C 
              CALL UCCOMP(M6REC(1:19),LABEL) 
              CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
              IF (LABTYPE.NE.LABB) THEN 
C 
C     Invalid label on the header record. 
C 
                WRITE(IOUERR,290) M6REC 
  290           FORMAT(/'*** ERROR: the following record is not ', 
     &          'allowed in the header of the external data file:', 
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
   98             CALL QUITER(0.,0,126,INERR) 
                  IF (RC.LT.0) THEN 
                    WRITE(IOUERR,298) 
  298               FORMAT('*** (Unexpected end of file.)') 
                  ELSE 
                    BACKSPACE(SOURCE) 
                    READ(SOURCE,199,ERR=98,END=98) LASTREC 
  199               FORMAT(A) 
                    WRITE(IOUERR,299) LASTREC 
  299               FORMAT('*** The following record caused the error:', 
     &                    /'*** ',A) 
                  ENDIF 
                  RETURN 
                ELSE 
                  BACKSPACE(SOURCE) 
C 
C      READ IN DATA. 
C 
                  READ(SOURCE, *, END=98, IOSTAT=RC) 
     &               ((STARTDIST(I,J), I=1,24), J=1,2) 
C 
C      NORMALIZE.  Sum of each 24 hours (weekday/weekend) should be one. 
C 
                  DISTSUM(1) = 0.0 
                  DISTSUM(2) = 0.0 
                  DO I=1,24 
                    DO J=1,2 
                    DISTSUM(J) = DISTSUM(J) + STARTDIST(I,J) 
                    END DO 
                  END DO 
                  DO I=1,24 
                    DO J=1,2 
                    STARTDIST(I,J) = STARTDIST(I,J) / DISTSUM(J) 
                    END DO 
                  END DO 
C 
                  IF (EXTDATA) CLOSE(SOURCE) 
C 
                END IF 
              END IF 
            END IF 
          END IF 
        END IF 
      END IF 
      RETURN 
      END 
