      SUBROUTINE GETFVMT(INPREC,NSCEN,INERR,RC) 
C 
C    GETFVMT reads in the user supplied VMT distribution across 
C    roadway type and hour. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C 16 Feb 01 @EPA-elg Changed GETFVMT input routine to allow input of FVMT as a 
C                    function of vehicle class. 
C 22 Nov 00 @EPA-bag bug301 added GFVMT to enable scenario and run level 
C           user input 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  30 May 2000 @EPA-elg, Major revision for input of VMT's 
C                        with free format 
C 10 Apr 00 AIR Task 01: replaced outask with iouin/iouout 
C 28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C           module subroutine CHKLAB was removed from that module 
C           and converted to an ordinary external procedure. 
C 11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C           calling sequence for that subroutine was changed. 
C  1 March 1999 @DynTel-HQ 2-698 Read in speed correction coefficient 
C                          from user input 
C 
C 
C  Input on call: 
C 
C    argument list: INPREC,NSCEN 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /SPEED9/ GLFVMT,FVMT 
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
C     Declare parameter list. 
C 
      CHARACTER,  INTENT(IN) :: INPREC*(*) 
      INTEGER,    INTENT(IN) :: NSCEN 
      INTEGER, INTENT(INOUT) :: INERR 
      INTEGER,   INTENT(OUT) :: RC 
C 
C     Declare local variables. 
C 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*40 
      INTEGER     LABTYPE 
      CHARACTER   M6REC*150 
      INTEGER     SECTION 
      CHARACTER   STDLAB*40 
      INTEGER     I,J,K 
      REAL        TSUM 
      INTEGER     LDGIMP 
      INTEGER     IVC 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('VMT BY FAC',STDLAB,LDGIMP,SECTION) 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file. 
C 

ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)

      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) TRIM(FILENAME) 
  120     FORMAT(/'* Reading Hourly Roadway VMT distribution from', 
     &              ' the following external',/'* data file: ',A) 
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
      DO K = 1,28 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) EXIT 
        BACKSPACE(IOUALT) 
C 
        READ (IOUALT,*,IOSTAT=RC) IVC 
        IF(RC.LT.0) EXIT 
        IF (IVC.GT.28 .OR. IVC.LT.1) THEN 
          CALL QUITER(0.,IVC,174,INERR) 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C 
C       Read the implementation schedule for the current vehicle type. 
C 
        READ (IOUALT,*,IOSTAT=RC) 
     &          IVC, ((FVMT(I,J,IVC), I=1,4), J=1,24) 
C 
        IF(RC.LT.0) EXIT 
C 
        DO J = 1,24 
          TSUM = 0.0 
          DO I = 1,4 
            TSUM = TSUM + FVMT(I,J,IVC) 
          END DO 
          IF(ABS(1.-TSUM)>0.0025) THEN 
            DO I = 1,4 
              FVMT(I,J,IVC) = FVMT(I,J,IVC) / TSUM 
            END DO 
          ENDIF 
        END DO 
C 
        IF(NSCEN==0) GLFVMT=FVMT ! for scenario/run level user input 
C 
      END DO 
C 
   75 WRITE(IOUREP,300) 
  300 FORMAT(/,2x, 
     &         'Reading User Supplied ROADWAY VMT Factors') 
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
C 
C 
  98  INERR = INERR + 1 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
