      SUBROUTINE GETHVMT(INPREC,NSCEN,INERR,RC) 
C 
C    GETHVMT reads in the user supplied VMT distribution 
C    across the day by hour. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C 22 Nov 00 @EPA-bag bug301 added GHVMT to enable scenario and run level 
C           user input 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  21 Sep 00 @EPA-elg  bug156 changed vmt by hour, fac to vmt by hour 
C  29 Jun 00 @EPA-elg  Added code to create a very small non zero value for HVMT if 
C                      it is entered as a zero by the user. 
C  30 May 2000 @EPA-elg, Subroutine GETHVMT was completely re-written 
C                        but still performs basically the same function as the 
C                        previous version of the same name. 
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
C    /SPEED9/ GLHVMT,HVMT 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C 
C  Notes: 
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
      INTEGER     J,ICC 
      REAL        TSUM 
      
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
        CALL CHKLAB('VMT BY HOUR',STDLAB,LABB,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), FILENAME)

ccs     All M6 input files live in SMK_M6PATH
        CALL GETENV( "SMK_M6PATH", M6DIR )
        M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME

        OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *         IOSTAT=RC,ERR=80) 
        IF (RC.NE.0) THEN 
C 
C     Error opening the external data file. 
C 
   80     WRITE(IOUERR,280) FILENAME 
  280     FORMAT(/'*** ERROR: unable to open the following external ', 
     &      'data file:',/'*** ',A) 
          INERR = INERR + 1 
        ELSE 
          WRITE(IOUREP,200) FILENAME 
  200     FORMAT(/'* Reading Hourly VMT distribution from', 
     &            ' the following external',/'* data file: ',A) 
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
            IF (LABTYPE.NE.LABB) THEN 
C 
C     Invalid label on the header record. 
C 
              WRITE(IOUERR,290) M6REC 
  290         FORMAT(/'*** ERROR: the following record is not ', 
     &        'allowed in the header of the external data file:', 
     &        /'*** ',A) 
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
                ICC = 0 
                READ(SOURCE, *, END=98, IOSTAT=RC) 
     &             (HVMT(J), J=1,24) 
C 
                TSUM = 0.0 
                DO J = 1,24 
                  IF(HVMT(J).LE.0.0) THEN 
                    HVMT(J) = 0.00001 
                    ICC = ICC + 1 
                  ENDIF 
                  TSUM = TSUM + HVMT(J) 
                END DO 
C 
                IF(ICC.GT.0) THEN 
                  WRITE(IOUERR,291) ICC 
  291             FORMAT(/'*** WARNING *** The user has entered ', 
     *             I2,' zero values for HVMT.',/, 
     *             'These values have been set to 0.00001 and ', 
     *             'the entire distribution re-normalized') 
                ENDIF 
C 
                DO J = 1,24 
                  HVMT(J) = HVMT(J) / TSUM 
                END DO 
C 
                IF(NSCEN==0) GLHVMT=HVMT ! for scenario/run level user input 
C 
                IF (EXTDATA) CLOSE(SOURCE) 
C 
              END IF 
            END IF 
          END IF 
        END IF 
      END IF 
      RETURN 
      END 
