      SUBROUTINE GETSVMT(INPREC,NSCEN,INERR,RC) 
C 
C  GETSVMT reads in, validates and stores VMT for speed 
C  distribution at different hours from an external data file. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C 22 Nov 00 @EPA-bag bug301 added GSVMT to enable scenario and run level 
C           user input 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C 31 May 00 @EPA-elg  Changed input to Free Format and Normalized VMT fractions. 
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
C    /SPEED9/ GLSVMT,SVMT 
C 
C  Local variable / array dictionary: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'SPEED9.I' 
C 
      INTEGER, EXTERNAL :: JUNIT   
C 
C     Declare parameter list. 
C 
      CHARACTER,  INTENT(IN) :: INPREC*(*) 
      INTEGER,    INTENT(IN) :: NSCEN 
      INTEGER, INTENT(INOUT) :: INERR 
      INTEGER,   INTENT(OUT) :: RC 
C 
C     Declare local variables. 
C 
      INTEGER     LABB 
      LOGICAL     EXTDATA 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
      INTEGER     SECTION 
      INTEGER     SOURCE 
      CHARACTER   STDLAB*40 
      REAL        TSUM 
      INTEGER     I, IFAC, IHR, ISPD 
      REAL        DUM(2,14,24), DUM1(14) 
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
 
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
C       Read the data from an external data file. 
C 
        EXTDATA = .TRUE. 
        SOURCE = JUNIT() 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
        CALL CHKLAB('SPEED VMT',STDLAB,LABB,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), M6FILE)

ccs     Usually M6 input files live in SMK_M6PATH; however, we put
ccs     the full path of the files in the M6 input file

        OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *         IOSTAT=RC,ERR=80) 
        IF (RC.NE.0) THEN 
C 
C     Error opening the external data file. 
C 
   80     WRITE(IOUERR,280) M6FILE 
  280     FORMAT(/'*** ERROR: unable to open the following external ', 
     &      'data file:',/'*** ',A) 
          INERR = INERR + 1 
        ELSE 
          WRITE(IOUREP,200) M6FILE 
  200     FORMAT(/'* Reading Hourly, Roadway, and Speed VMT dist.', 
     &              ' from the following external',/'* data file: ',A) 
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
                WRITE(IOUERR,285) M6FILE 
  285           FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
                RC = 1 
                INERR = INERR + 1 
              CASE (1:) 
                WRITE(IOUERR,286) M6FILE 
  286           FORMAT(/'*** ERROR reading the header of the ', 
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
C     For a given HOUR and ROADWAY Type the VMT distribution 
C     of the 14 Speeds should sum to 1.00 
C 
                DO 34 I = 1, 48 
                  READ(SOURCE, *, END=98, IOSTAT=RC) 
     *            IFAC, IHR, (DUM1(ISPD), ISPD=1,14) 
C 
                  IF(IFAC.EQ.1 .OR. IFAC.EQ.2) THEN 
                    IF(IHR.GE.1 .AND. IHR.LE.24) THEN 
                      DO ISPD = 1,14 
                        DUM(IFAC,ISPD,IHR) = DUM1(ISPD) 
                      END DO 
                    ELSE 
          WRITE(IOUERR,*) '*** ERROR: Illegal input for SVMT Hour', 
     *                IHR 
                      INERR = INERR + 1 
                    ENDIF 
                  ELSE 
          WRITE(IOUERR,*) '*** ERROR: Illegal input for SVMT Roadway', 
     *              IFAC 
                    INERR = INERR + 1 
                  ENDIF 
   34           CONTINUE 
C 
C 
                DO 35 IFAC = 1,2 
                  DO 36 IHR = 1,24 
                    TSUM = 0.0 
                    DO 37 ISPD = 1,14 
                      TSUM = TSUM + DUM(IFAC,ISPD,IHR) 
   37               END DO 
C 
                    DO 38 ISPD = 1,14 
                      SVMT(IFAC,ISPD,IHR) = DUM(IFAC,ISPD,IHR) / TSUM 
   38               END DO 
   36             END DO 
   35           END DO 
C 
C 
                IF(NSCEN==0) GLSVMT=SVMT ! for scenario/run level user input 
                IF (EXTDATA) CLOSE(SOURCE) 
C 
              END IF 
            END IF 
          END IF 
        END IF 
      END IF 
      RETURN 
      END 
