      SUBROUTINE GET94E(INPREC,INERR,RC) 
 
C  GET94E reads in, validates and stores user-supplied emission factor 
C  records for 1994+ LDG vehicles. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, UCCOMP, NXTREC, QUITER. 
C 
C  Changes: (Last change first) 
C 
C  11 Oct 00 AIR Task 02: Bug Fix 275: Expand free format reads 
C  03 May 00 @EPA-BG  Changed QUITER message numbers to correspond to messages in 
C            QUITER 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  14 Sep 99 @DynTel-ddj 2-000 Changed check on "deterioration rate" 
C            and "zero mile level" from <= 0.0 to .LT. 0.0.  Added 
C            ISTD index to USEREFS array. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C 
C  Input on call: 
C 
C    argument list:  INPREC 
C    common blocks: 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type Description 
C  ------   ---  ------------------------------------------------------- 
C  DATAFIELD C   One field from the user input record. 
C  DR1       R   The deterioration rate entered by the user. 
C  EF94      C   Numeric code for the 94+ EMISSION FACTORS input label. 
C  FILENAME  C   Name of external data file. 
C  INH       I   Code for the high/normal flag entered by the user. 
C  IOS       I   Return code from Fortran I/O operation. 
C  IP        I   Code for the pollutant entered by the user. 
C  ISR       I   Code for the start/running flag entered by the user. 
C  ISTD      I   The deterioration rate entered by the user. 
C  IV        I   Code for the vehicle type entered by the user. 
C  LABEL     C   Label field from a labeled input record 
C  LABTYPE   I   Numeric token identifying the label type 
C  LASTREC   C   The previous record read from the input file 
C  M6REC     C   Record from the external data file 
C  MYF       I   The 'first model year' entered by the user. 
C  MYL       I   The 'last model year' entered by the user. 
C  PTR1      I   Used to parse the user input record. 
C  PTR2      I   Used to parse the user input record. 
C  SECTION   I   Dummy argument returned by CHKLAB. 
C  STDLAB    C   The standard (uppercase, left justified) input label. 
C  TESTFIELD C   Left justified, uppercase version of DATAFIELD. 
C  ZML       R   The zero mile level enterd by the user. 
C  X         R   Temporary real to integer holder. 
C 
C  Notes: 
C 
C 
      USE EFUSER94 
 
      IMPLICIT NONE 
 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'YEARS4.I' 
  
      INTEGER, EXTERNAL :: JUNIT
  
      CHARACTER(LEN=20)    ::  DATAFIELD 
      CHARACTER(LEN=*)     ::  INPREC 
      INTEGER              ::  INERR 
      INTEGER              ::  RC 
 
      REAL                 ::  DR1 
      INTEGER              ::  EF94 
      CHARACTER(LEN=80)    ::  FILENAME 
      INTEGER              ::  INH 
      INTEGER              ::  IOS 
      INTEGER              ::  IP 
      INTEGER              ::  ISR 
      INTEGER              ::  ISTD 
      INTEGER              ::  IV 
      CHARACTER(LEN=19)    ::  LABEL 
      INTEGER              ::  LABTYPE 
      CHARACTER(LEN=150)   ::  LASTREC 
      CHARACTER(LEN=150)   ::  M6REC 
      INTEGER              ::  MYF 
      INTEGER              ::  MYL 
      INTEGER              ::  PTR1 
      INTEGER              ::  PTR2 
      INTEGER              ::  SECTION 
      CHARACTER(LEN=40)    ::  STDLAB 
      CHARACTER(LEN=20)    ::  TESTFIELD 
      REAL                 ::  ZML 
      REAL                 ::  X 
 
      TYPE (EMISSREC), POINTER  ::  TEMP 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  
      ! Find the token number of the record label that must appear 
      ! on a record somewhere in the header of the external data file. 
  
      CALL CHKLAB('94+ LDG EMISS FACT',STDLAB,EF94,SECTION) 
      ! The next two assignment statements are present solely to suppress 
      ! compiler generated warnings. They have no functional importance. 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
  
      ! Open the external data file. 
  
ccs      CALL UCCOMP(INPREC(21:), FILENAME) 
      CALL LJCOMP(INPREC(21:), FILENAME)

ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME

      IOUALT = JUNIT()

      OPEN(IOUALT,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     &      IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading 94+ LDG EMISSION FACTORS from the ', 
     &        'following external',/'* data file: ',A) 
  
      ! Read the required header line from the file. 
  
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC /= 0) THEN 
        GOTO 85 
      ENDIF 
  
      ! Check for a valid labeled input record.  
 
      CALL UCCOMP(M6REC(1:19),LABEL) 
      CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
      IF (LABTYPE /= EF94) THEN 
        GOTO 90 
      ENDIF 
 
      ! Allocate the array to hold the emission factors, if necessary. 
      IF (.NOT. ALLOCATED(USEREFS)) THEN 
 
        ALLOCATE(USEREFS(5,3,2,2,10),STAT=RC) 
 
        IF (RC /= 0) THEN 
          CALL QUITER(0.0,0,536,INERR) 
          GOTO 98 
        ENDIF 
 
        DO INH = 1,2 
          DO ISR = 1,2 
            DO IP = 1,3 
              DO IV= 1,5 
                DO ISTD = 1, 10 
                  NULLIFY(USEREFS(IV,IP,INH,ISR,ISTD)%ROOT) 
                END DO 
              END DO 
            END DO 
          END DO 
        END DO 
 
      END IF 
   
      ! Read the rest of the data file, storing emisison factors 
      ! in the USEREFS array. The only exit from the loop is a READ  
      ! error or EOF on the input file. 
  
      DO 
  
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC /= 0) GOTO 95 
  
        ! Look for the vehicle type, the first data field on the record. 
        ! Initialize by setting PTR2 to 0 (just before the first character 
        ! in the record). PTR1 is set by the NXTTOK subroutine. 
  
        PTR2 = 0 
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        IF (TESTFIELD == 'LDGV') THEN 
          IV = 1 
        ELSEIF (TESTFIELD == 'LDGT1') THEN 
          IV = 2 
        ELSEIF (TESTFIELD == 'LDGT2') THEN 
          IV = 3 
        ELSEIF (TESTFIELD == 'LDGT3') THEN 
          IV = 4 
        ELSEIF (TESTFIELD == 'LDGT4') THEN 
          IV = 5 
        ELSE 
          CALL QUITER(0.0,0,537,INERR) 
          RC = 4 
        ENDIF  
 
        ! Next, look for the pollutant. 
           
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        IF (TESTFIELD == 'HC') THEN 
          IP = 1 
        ELSEIF (TESTFIELD == 'CO') THEN 
          IP = 2 
        ELSEIF (TESTFIELD == 'NOX') THEN 
         IP = 3 
        ELSE 
          CALL QUITER(0.0,0,538,INERR) 
          RC = 6 
        ENDIF 
 
        ! Next, the start/running flag. 
 
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        IF (TESTFIELD == 'RUN') THEN 
          ISR = 1 
        ELSEIF (TESTFIELD == 'START') THEN 
          ISR = 2 
        ELSE 
          CALL QUITER(0.0,0,539,INERR) 
          RC = 7 
        ENDIF  
 
        ! Next, the high/normal flag. 
          
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        IF (TESTFIELD == 'HIGH') THEN 
          INH = 2 
        ELSEIF (TESTFIELD == 'NORMAL') THEN 
          INH = 1 
        ELSE 
          CALL QUITER(0.0,0,540,INERR) 
          RC = 8 
        ENDIF 
 
        ! Next, the standards category. 
 
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        READ(TESTFIELD,*,ERR=92,IOSTAT=IOS) X 
        ISTD=INT(X) 
        IF (ISTD < 1 .OR. ISTD > 10) THEN 
          CALL QUITER(0.0,0,541,INERR) 
          RC = 9 
        ENDIF 
 
        ! Next, the first model year. 
          
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        READ(TESTFIELD,*,ERR=92,IOSTAT=IOS) X 
        MYF=INT(X) 
        IF (MYF < 1994 .OR. MYF > IYEND) THEN 
          CALL QUITER(0.0,0,29,INERR) 
          RC = 10 
        ENDIF 
  
        ! Next, the last model year. 
  
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        READ(TESTFIELD,*,ERR=92,IOSTAT=IOS) X 
        MYL=INT(X) 
        IF (MYL < 1994 .OR. MYL > IYEND) THEN 
          CALL QUITER(0.0,0,29,INERR) 
          RC = 11 
        ENDIF 
        IF (MYL < MYF) THEN 
          CALL QUITER(0.0,0,29,INERR) 
          RC = 12 
        ENDIF 
 
        ! Next, zero mile level.  
 
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        READ(TESTFIELD,*,ERR=92,IOSTAT=IOS) ZML 
        IF (ZML .LT. 0.0) THEN 
          CALL QUITER(0.0,0,542,INERR) 
          RC = 13 
        ENDIF 
 
        ! Next, the first deterioration rate. 
 
        CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        CALL UCCOMP(DATAFIELD,TESTFIELD) 
        READ(TESTFIELD,*,ERR=92,IOSTAT=IOS) DR1 
        IF (DR1 .LT. 0.0) THEN 
          CALL QUITER(0.0,0,543,INERR) 
         RC = 14 
        ENDIF 
 
        ! Finally, the second deterioration rate. 
 
        ! CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
        ! CALL UCCOMP(DATAFIELD,TESTFIELD) 
        ! READ(TESTFIELD,'(F9.4)',ERR=92,IOSTAT=IOS) DR2 
        ! IF (DR2 <= 0.0) THEN 
        !   CALL QUITER(0.0,0,544,INERR) 
        !   RC = 15 
        ! ENDIF 
 
        ! If all values were read correctly, check for duplicates. 
 
        IF (RC == 0) THEN 
          TEMP => USEREFS(IV,IP,INH,ISR,ISTD)%ROOT 
        ELSE 
          NULLIFY(TEMP) 
        ENDIF 
 
        DO WHILE(ASSOCIATED(TEMP) .AND. RC == 0) 
          IF (MYF <= TEMP%MYLAST .AND. MYL >= TEMP%MYFIRST .AND. 
     &        ISTD == TEMP%STDCAT) THEN 
            CALL QUITER(0.0, 0, 29, INERR) 
            RC = 16 
          ENDIF 
          TEMP => TEMP%NEXT 
        END DO  ! WHILE (ASSOCIATED .AND. RC == 0) 
 
        IF (RC /= 0) THEN 
          WRITE(IOUERR,140) RC, M6REC 
 140      FORMAT(' *** ERROR',I3,' reading the following record', 
     &           ' from the external data file:',/' *** ',A) 
        ENDIF 
 
        ! If everything is OK, store the record. 
 
        IF (INERR /= 0) CYCLE 
 
        ALLOCATE (TEMP) 
        TEMP = EMISSREC (MYF, MYL, ISTD, ZML, DR1, 
     &                   USEREFS(IV,IP,INH,ISR,ISTD)%ROOT) 
        USEREFS(IV,IP,INH,ISR,ISTD)%ROOT => TEMP 
 
      END DO   ! FOREVER 
 
      ! ERROR handling code follows. 
  
      ! Error opening the external data file. 
  
   80 WRITE(IOUERR,180) FILENAME 
  180 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      RC = 1 
      GOTO 99 
  
      ! Error reading the header of the external data file. 
  
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
          RC = 2 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
        CASE DEFAULT 
      END SELECT 
      CLOSE(IOUALT) 
      INERR = INERR + 1 
      GOTO 98 
  
      ! Invalid label on the header record. 
  
   90 WRITE(IOUERR,190) M6REC 
  190 FORMAT(/'*** ERROR: the following record is not allowed', 
     &        ' in the header of the external data file:', 
     &       /'*** ',A) 
      CLOSE(IOUALT) 
      RC = 3   
      INERR = INERR + 1 
      GOTO 98 
 
      ! Fortran numeric I/O error reading internal buffer DATAFIELD. 
  
   92 WRITE(IOUERR,192) IOS, M6REC 
  192 FORMAT(/'*** ERROR reading a data field from the', 
     &        ' following labeled input record:', 
     &       /'*** ',A) 
      RC = IOS 
      INERR = INERR + 1 
      GOTO 99 
 
      ! Error reading the data section of the file. 
 
   95 IF (RC < 0) THEN 
        RC = 0 
      ELSE 
        CALL QUITER(0.,0,128,INERR) 
        BACKSPACE(IOUALT) 
        READ(IOUALT,195,ERR=98,END=98) LASTREC 
  195   FORMAT(A) 
        WRITE(IOUERR,196) LASTREC 
  196   FORMAT('*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
 
   98 CLOSE(IOUALT) 
   99 RETURN 
      END 
