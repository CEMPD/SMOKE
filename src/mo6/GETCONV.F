      SUBROUTINE GETCONV(INPREC,INERR,RC) 
C 
C  GETCONV reads in, validates and stores conversion factors 
C  from an external data file. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, NXTREC, QUITER, UCCOMP. 
C 
C  Changes : (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C 26 Sep 00 @EPA-EG      Added HDCFAC Flag Structure. 
C 05 Jun 00 @EPA-EG      Revised for user can input values in free format 
C 11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C           so that Lahey Standard for keyboard/monitor I/O is followed. 
C 28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C           module subroutine CHKLAB was removed from that module 
C           and converted to an ordinary external procedure. 
C 11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C           calling sequence for that subroutine was changed. 
C  3 Mar 99 @DynTel-MLA 2-659 Updated comments, var. declarations, etc. 
C 15 Oct 1998 @DynTel-BK 2-000 Changes in write statements. Made  
C             format statement label conform to 2?? for write  
C             statements and 1??) for read statements.  Also, 
C             removed from comment statement the variables not used  
C             in the HDCCOM.I file.  
C  7 Oct 1998 @DynTel-HQ 2-678 Removed inclusion of FLAGS2.I 
C  6 Oct 1998 @DynTel-ZK 2-000 Explicit Typing 
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
C    /FLAGS1/ CONVFLG 
C    /HDCCOM/ HDCFAC 
C 
C  Local variable / array dictionary: 
C 
C  Name       Type Description 
C  --------   ---- --------------------------------------------------- 
C  CONVFACTOR  I   Heavy duty conversion factor read from the input. 
C  EXTDATA     L   Has the external data file been opened? 
C  FILENAME    C   Name of the external data file. 
C  I           I   Loop control variable. 
C  J           I   Loop control variable. 
C  K           I   Loop control variable. 
C  LABEL       C   The label on the labeled input record. 
C  LABTYPE     I   The numeric code for the input label type. 
C  LASTREC     C   The previous record read from the input. 
C  M6REC       C   The current input record 
C  SECTION     I   Dummy argument for CHKLAB. 
C  SOURCE      I   Logical unit number for the input file. 
C  STDLAB      C   Dummy argument for CHKLAB. 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'HDCCOM.I' 
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
      INTEGER    IVC2 
      REAL       DUM1(56) 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
 
      IF (CONVFLG.EQ.1) THEN 
        CONVFLG = 2 
      ELSE 
        GOTO 95 
      ENDIF 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('HD CONV FACTORS',STDLAB,LDGIMP,SECTION) 
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
  120 FORMAT(/'* Reading Heavy-duty Conversion Factors ', 
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
C     Initialize the user supplied conversion factor 
C     array to the default values because the user is not 
C     required to enter a complete set of conversion factors. 
C 
      DO J = 1,56 
        DO IVC2 = 1,19 
          HDCFAC2(J,IVC2) = HDCFAC(J,IVC2) 
        END DO 
      END DO 
C 
      DO 31 I = 1,500 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read the conversion factor for the current vehicle type. 
C 
        READ (IOUALT,*,IOSTAT=RC) 
     &          IVC, (DUM1(J),J=1,56) 
C 
        HDCFAC_FLG = 1 
        IF(IVC.GE.1 .AND. IVC.LE.19) THEN 
            DO J = 1,56 
              HDCFAC2(J,IVC) = DUM1(J) 
            END DO 
 
        ELSE 
          WRITE(IOUERR,*) '*** ERROR: Illegal input for Heavy-duty ', 
     *                    'Conversion Factors', IVC 
          INERR = INERR + 1 
        ENDIF 
 
        IF(RC.LT.0) EXIT 
C 
   31 END DO 
C 
  210 FORMAT(/'*** WARNING: The Heavy Duty Conversion Factor', 
     &        ' Is Out Of Bounds(0.00-5.00)', 
     &      /'*** ',F10.3) 
      DO 38 I=1, 19 
        DO 37 J=1,56 
          IF (HDCFAC(J, I).LT.0.000.OR.HDCFAC(J, I).GT.5.000) THEN 
              WRITE (IOUERR,210) HDCFAC(J, I) 
          ENDIF 
   37   CONTINUE 
   38 CONTINUE 
C 
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
