      SUBROUTINE WTLRNL(IWK,INPREC,NSCEN,INERR,RC) 
C 
C  WTLRNL reads in the 6 user-supplied %VMT trip length values. 
C  WTLRNL previously read in 3 weighting values (never used) for TLVMT. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C 28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C           module subroutine CHKLAB was removed from that module 
C           and converted to an ordinary external procedure. 
C 11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C           calling sequence for that subroutine was changed. 
C  2  Mar 99 @DynTel-BAG 2-664  Expanded TLVMTU array for hourly 
C            running losses (weekly and hourly indices); deleted RNLFLG. 
C  15 January 1999 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IWK,NSCEN,INERR 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUGEN 
C    /RUNLS1/ MAXTL 
C 
C  Output on return: 
C 
C    argument list: INERR,RC,INPREC 
C    common blocks: 
C    /RUNLS1/ GLVMTU,TLVMTU 
C 
C  Local array subscripts 
C 
C  TLVMTU(6)  -  TLVMTU( ITL ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  SECTION  I    Dummy argument for CHKLAB. 
C  STDLAB   C    Dummy argument for CHKLAB. 
C  TLSUM    R    sum of the trip length VMT fractions 
C 
C  Notes: 
C 
C  WTLRNL was activated beginning with Version D02 of MOBILE4.1 
C 
C  WTLRNL is not used in Release or In-House MOBILE4. 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C  Dec. 16 97 Dyntel@MMS 2-645 renamed SPDFLG to RNLFLG, 
C             reduced value of RNLFLG by 1. 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'RUNLS1.I' 
C  
      INTEGER, EXTERNAL :: JUNIT    
C 
C     declare parameter list 
C 
      INTEGER, INTENT(INOUT) :: INERR 
      INTEGER, INTENT(IN) :: NSCEN 
      INTEGER, INTENT(IN) :: IWK 
      INTEGER, INTENT(OUT):: RC 
      CHARACTER INPREC*(*) 
C 
C     declare local variables 
C 
      CHARACTER(LEN=80) :: FILENAME 
      CHARACTER(LEN=19) :: LABEL 
      CHARACTER(LEN=40) :: STDLAB 
      CHARACTER(LEN=150) :: LASTREC 
      CHARACTER(LEN=150) :: M6REC 
      INTEGER   IH 
      INTEGER   ITL 
      INTEGER   LABTYPE 
      INTEGER, PARAMETER :: MAXTL_IH=14 
      INTEGER   SECTION 
      INTEGER   SOURCE 
      INTEGER   WKDA_DIST 
      INTEGER   WKND_DIST 
      LOGICAL   EXTDATA 
      REAL      TLSUM 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=50) M6DIR 
      CHARACTER(LEN=80) M6FILE 
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
        SOURCE = JUNIT() 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
        IF( IWK .EQ. 1 ) CALL 
     *      CHKLAB('WEEK DAY TRIP LENGTH DISTRIBUTIONS',STDLAB, 
     *       WKDA_DIST,SECTION) 
        IF( IWK .EQ. 2 ) CALL 
     *      CHKLAB('WEEK END TRIP LENGTH DISTRIBUTIONS',STDLAB, 
     *      WKND_DIST,SECTION) 
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
     *       IOSTAT=RC,ERR=80) 
        EXTDATA = .TRUE. 
        IF( IWK .EQ. 1) THEN 
          WRITE(IOUREP,130) FILENAME 
  130     FORMAT(/'* Reading non-default WEEKDAY RUNNING LOSS HOURLY ', 
     &          'TRIP LENGTH FRACTIONS'/'* from the following external', 
     &          ' data file: ',A/) 
        END IF 
        IF( IWK .EQ. 2) THEN 
        WRITE(IOUREP,131) FILENAME 
  131     FORMAT(/'* Reading non-default WEEKEND RUNNING LOSS HOURLY ', 
     &          'TRIP LENGTH FRACTIONS'/'* from the following external', 
     &          ' data file: ',A/) 
        END IF 
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
        IF ( (IWK .EQ. 1 .AND. LABTYPE .NE.WKDA_DIST) .OR. 
     *       (IWK .EQ. 2 .AND. LABTYPE .NE.WKND_DIST) ) THEN 
          GOTO 90 
        ENDIF 
C 
      ENDIF 
C 
C 
C     Read in alternate %VMT trip length fractions. 
C 
C     Read the data field from the input record. 
C     And move to the first data line. 
C 
      CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
      IF (RC.NE.0) GOTO 98 
      BACKSPACE(SOURCE) 
C 
      READ(SOURCE,100,END=97,ERR=98) 
     *               ((TLVMTU(ITL,IH,IWK), ITL=1,MAXTL),IH=1,14) 
  100 FORMAT(6(1X,F5.2)) 
C 
C     Verify %VMT trip length fraction values and that they sum to 1.0 
C     for each hourly interval, IH. 
C 
      DO IH = 1,MAXTL_IH 
C 
        TLSUM=0.0 
C 
        DO ITL=1,MAXTL 
C 
          IF( (TLVMTU(ITL,IH,IWK).LT.0.0) .OR. 
     *        (TLVMTU(ITL,IH,IWK).GT.100.0) ) THEN 
            CALL QUITER(TLVMTU(ITL,IH,IWK),0,93,INERR) 
          END IF 
C 
          IF( (TLVMTU(ITL,IH,IWK).LT.0.0) .OR. 
     *        (TLVMTU(ITL,IH,IWK).GT.100.0) ) THEN 
            TLVMTU(ITL,IH,IWK)=0.0 
          END IF 
          TLVMTU(ITL,IH,IWK)=TLVMTU(ITL,IH,IWK)/100.0 
          TLSUM=TLSUM+TLVMTU(ITL,IH,IWK) 
C 
        END DO 
C       end loop over travel lengths 
C       Check that TLSUM = 1.0 
C 
        IF( ((TLSUM+1.E-3).LT.1.) .OR. ((TLSUM-1.E-3).GT.1.) ) 
     *             CALL QUITER(TLSUM,0,94,INERR) 
        IF(TLSUM.LE.0.0) GO TO 99 
C 
      END DO 
C     end loop over hourly intervals 
C     Normalize TLVMTU and set up the global data array, if the data was 
C     entered in the global data section. 
C 
      DO IH = 1,MAXTL_IH 
        DO ITL=1,MAXTL 
          TLVMTU(ITL,IH,IWK)=TLVMTU(ITL,IH,IWK)/TLSUM 
          IF(NSCEN.EQ.0) GLVMTU(ITL,IH,IWK) = TLVMTU(ITL,IH,IWK) 
        END DO 
      END DO 
C 
      IF (EXTDATA) CLOSE(SOURCE) 
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
C     Error: alternate conversion factors have been 
C            read in. 
C     
C   96 WRITE(IOUERR,296) 
C  296 FORMAT(/'*** ERROR: Alternate conversion factors ', 
C     &        'have already been read in this run.', 
C     &       /'*** Only one set of conversion factors ', 
C     &        'can be entered in one run.'/) 
CCC  &      /'*** The factors entered here will replace the ', 
CCC  &       'rates that were entered previously.'/) 
C      RC = 3 
C      INERR = INERR + 1 
C      GOTO 99 
C 
C 
   97 RC = 1 
      GOTO 99 
   98 RC = -1 
C 
      CALL QUITER(0.,IWK,97,INERR) 
      IF (RC.LT.0) THEN 
        WRITE(IOUERR,298) 
  298   FORMAT('*** (Unexpected end of file.)') 
      ELSE 
        BACKSPACE(IOUGEN) 
        READ(INPREC,199,ERR=99,END=99) LASTREC 
  199   FORMAT(A) 
        WRITE(IOUERR,299) LASTREC 
  299   FORMAT('*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
C 
   99 RETURN 
      END 
