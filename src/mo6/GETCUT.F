      SUBROUTINE GETCUT(INPREC,INERR,RC) 
C 
C  GETCUT processes the I/M CUTPOINTS labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C  15 Sep 00 AIR Task 08: Expanded number of I/M programs to handle EVAP 
C  31 May 2000 @EPA-elg  Fixed READ Format Changed to Free READ format. 
C  22 May 00 @EPA-elg  Added code to restrict I/M cutpoints to final levels. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C  10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUGEN,IOUREP 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /IMPAR4/ CUTHC,CUTCO,CUTNO 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  IMPGM      I   The I/M program number which uses this data. 
C  MY         I   Implicit DO loop index (model year). 
C  IV         I   Implicit DO loop index (vehicle type). 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  SECTION    I   Dummy argument for CHKLAB. 
C  STDLAB     C   Dummy argument for CHKLAB. 
C  X          R   Temporary real to integer variable 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IOUCOM.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      INTEGER        CUTPOINTS 
      CHARACTER*80   DATAFIELD 
      LOGICAL        EXTDATA 
      CHARACTER*80   FILENAME 
      INTEGER        IMPGM 
      INTEGER        IMY 
      INTEGER        IOS 
      CHARACTER*19   LABEL 
      INTEGER        LABTYPE 
      INTEGER        MY 
      CHARACTER*150  M6REC 
      INTEGER        SECTION 
      INTEGER        SOURCE 
      CHARACTER*40   STDLAB 
      INTEGER        IVCUT 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
C 
C ..................................................................... 
C 
      RC = 0 
      EXTDATA = .FALSE. 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,IOSTAT=IOS,ERR=95) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program and check off 
C     the I/M CUTPOINTS data item for this I/M program. 
C 
      IMCHECK(1,IMPGM) = 2 
      IMCHECK(7,IMPGM) = 2 
C 
C     Read the cutpoints for the I/M program. 
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
        SOURCE = IOUALT 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
        CALL CHKLAB('I/M CUTPOINTS',STDLAB,CUTPOINTS,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       Open the external data file. 
C 
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 

ccs        CALL UCCOMP(DATAFIELD,FILENAME) 
        CALL LJCOMP(DATAFIELD,FILENAME)

        OPEN(SOURCE,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *       IOSTAT=RC,ERR=80) 
        EXTDATA = .TRUE. 
        WRITE(IOUREP,130) FILENAME 
  130   FORMAT(/'* Reading non-default I/M CUTPOINTS from the ', 
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
        IF (LABTYPE.NE.CUTPOINTS) THEN 
          GOTO 90 
        ENDIF 
C 
      ENDIF 
C 
C     Read 4 sets of cutpoints, one for each vehicle category. 
C 
      DO IVCUT = 1,4 
C 
C       Find the next data line. 
C     
        CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) GOTO 95 
        BACKSPACE(SOURCE) 
C 
C       Read the block of data for the current vehicle type. 
C 
        READ(SOURCE,*,IOSTAT=IOS,ERR=95) 
     &          (CUTHC(MY,IVCUT,IMPGM),MY=1,25), 
     &          (CUTCO(MY,IVCUT,IMPGM),MY=1,25), 
     &          (CUTNO(MY,IVCUT,IMPGM),MY=1,25) 
CCC  140   FORMAT((10(F7.3,1X),/10(F7.3,1X),/5(F7.3,1X))) 
C 
C 
        DO IMY = 1,25 
         IF(CUTHC(IMY,IVCUT,IMPGM) .LT. 0.80) THEN 
           WRITE(IOUERR,145) CUTHC(IMY,IVCUT,IMPGM),IMY,IVCUT 
 145       FORMAT('User Supplied HC IM240 Cutpoint of ',F6.3, 
     *            ' for year ',I3,' and Vehicle Type ',I3, 
     *            ' is Too Low and has been Reset to 0.80 g/mi ') 
           CUTHC(IMY,IVCUT,IMPGM) = 0.80 
         ENDIF 
C 
         IF(CUTCO(IMY,IVCUT,IMPGM) .LT. 15.0) THEN 
           WRITE(IOUERR,146) CUTCO(IMY,IVCUT,IMPGM),IMY,IVCUT 
 146       FORMAT('User Supplied CO IM240 Cutpoint of ',F6.3, 
     *            ' for year ',I3,' and Vehicle Type ',I3, 
     *            ' is Too Low and has been Reset to 15.0 g/mi ') 
           CUTCO(IMY,IVCUT,IMPGM) = 15.0 
         ENDIF 
C 
         IF(CUTNO(IMY,IVCUT,IMPGM) .LT. 2.0) THEN 
           WRITE(IOUERR,147) CUTNO(IMY,IVCUT,IMPGM),IMY,IVCUT 
 147       FORMAT('User Supplied NOx IM240 Cutpoint of ',F6.3, 
     *            ' for year ',I3,' and Vehicle Type ',I3, 
     *            ' is Too Low and has been Reset to 2.00 g/mi ') 
           CUTNO(IMY,IVCUT,IMPGM) = 2.0 
         ENDIF 
C 
        END DO 
C 
      END DO 
C 
      IF (EXTDATA) THEN 
        CLOSE(SOURCE) 
      ENDIF 
C 
      GOTO 99 
C 
C     Error opening the external data file. 
C 
   80 WRITE(IOUERR,180) FILENAME 
  180 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'*** ',A) 
          RC = 1 
          INERR = INERR + 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
          RC = 2 
        CASE DEFAULT 
      END SELECT 
      GOTO 99 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) M6REC 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 3   
      INERR = INERR + 1 
      GOTO 99 
C 
C Error or EOF on any attempted read. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
      IF (EXTDATA) THEN 
        CLOSE(SOURCE) 
      ENDIF 
C 
   99 RETURN 
      END 
