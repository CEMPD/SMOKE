       SUBROUTINE PROCHDR(INPFILE,CMDOUTPUT,CMDACTION,INERR,RC) 
C 
C  PROCHDR reads the header section of a Mobile6 input file. 
C 
C  Called by DRIVER. 
C 
C  Calls DBOPTS,DBSELE,DBSELF,DBSELP,DBSELV,LJCOMP,NXTREC,NXTTOK,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 16 Jan 02 AIR Task 27: Improved error handling via INERR updates and 
C           inclusion in parameter list. 
C 14 Jan 02 EPA-bag bug376: Format statement in FIELDNAMES output for aggregated output 
C                        and remove space, 1X, in FORMAT statement for all FIELDNAMES. 
C 13 Dec 01 AIR Task 27: Changed CO2 so that it is not selected by default. 
C 27 Nov 01 AIR Task 28: Changed SPREADSHEET to apply to all pollutants. 
C           Added dynamic array call for AT. 
C 28 Sep 01 AIR Task 25: Changed TOXICS SPREADSHEET to SPREADSHEET to 
C           accomodate both AT and PM outputs. Also changed spreadsheet 
C           naming methodology. 
C 30 Aug 01 AIR Task 22: Added TOXICS SPREADSHEET 
C 24 Aug 01 AIR Task 15-21: Added warning that the hourly database  
C           output for each AIR TOXICS scenario is 100+mb.  
C 16 Aug 01 AIR Task 15-21: Changed suffix from .AT to .TOX for Air Toxics. 
C 30 Jan 01 AIR Task 11: Added PARTICULATES and AIR TOXICS headers and files. 
C  29 Nov 00 @EPA-BG  renamed database age to database ages 
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection 
C 19 Oct 00 @EPA-djb Task W2, Added AGGREGATED OUTPUT command. 
C 20 Oct 00 AIR Task 03: Added MYR to daily and hourly output 
C 12 Sep 00 @EPA-djb Bug 245, Initialized DBSELFAC, DBSELPOL, DBSELVEH and 
C           DBSELEFT arrays. 
C 14 Sep 00 @EPA-elg bug238  Output of database tables TB2 thru TB7 has been disabled. 
C 28 Jun 00 AIR Task 04: Added EVAP_FLAG processing and warning. Added warning 
C                        about NO DESCRIPTIVE OUTPUT. 
C 15 Jun 00 @EPA-ddj,  Bug120 Removed unused SSACTION variable 
C  8 Jun 00 @EPA-ddj,  Bug120 Removed code associated with spreadsheet 
C                      output. 
C 18 May 00 @EPA-bag   bug#93 removed carriage control from FORMAT statement #'s 280 and 285 
C                      placed carriage return in 'WRITE(IOUOUT,'(" "/A)') LINE2' 
C                      statements following the 280 and 285 FORMAT statements 
C 16 May 00 @EPA-bag   fixed label problem with NO DESC OUTPUT label; 
C                      CHKLAB needs a match for each set of characters 
C                      between the shortened and long version of a label; 
C                      the shortened version is a lower character limit 
C                      which cannot have a blank as the lower character 
C                      limit for any of the character sets in a label 
C                      (character sets are separated by a blank or space) 
C 05 May 00 AIR Task 02: Added NO DESC OUTPUT to input file processing. 
C 05 MAY 00 @EPA TASK E22: ADD DATABASE DAILY AND HEADER RECORD OPTIONS 
C 11 Apr 00 AIR Task 02: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C           so that Lahey Standard for keyboard/monitor I/O is followed. 
C 10 Apr 00 AIR Task 02 Changed DATABASE POLLUTANTS to POLLUTANTS 
C 10 Mar 00 @DynTel-MLA 1-011  Database output 
C 23 Aug 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPFILE, CMDOUTPUT, CMDACTION, INERR, RC 
C    common blocks: 
C    /IOUCOM/ IOUOUT, IOUGEN 
C 
C  Output on return: 
C 
C    argument list:  INERR, RC 
C    common blocks: 
C    /DATABASE/  DBACTION,DBFNAMES,DBSELEFT,DBSELFAC,DBSELPOL,DBSELVEH 
C                DBHEAD, DBDAILY, DBAGGR 
C    /IOUCOM/    WMOUT 
C    /OPCNTL/    NODESC 
C    /FLAGS5/    EVAP_FLAG 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type                      Description 
C  ------     ----  --------------------------------------------------------------- 
C  BACK         L   Tells intrinsic IDX to search backw. from end of string. 
C  DATAFIELD    C   One blank-delimited data item from the input record. 
C  DBTAB        I   Database output table number. 
C  EXTENSION    C   The file name extension of the input file name. 
C  FIELD1       C   One data field from the input record. 
C  FIELD2       C   One data field from the input record. 
C  IDB          I   Temporary flag/do loop 
C  IDX          I   Pointer to position in a character string. 
C  INERR        I   Count of errors returned. 
C  INPLINE      C   A line from the input file. 
C  IOS          I   Return code from FORTRAN I/O operation. 
C  LABEL        C   Label field of the input record. 
C  LABTYPE      I   Numeric code that identifies the label. 
C  LINE1        C   Buffer used to compose error messages. 
C  LINE2        C   Buffer used to compose error messages. 
C  PTR1         I   Pointer used in NXTTOK to parse the input record. 
C  PTR2         I   Pointer used in NXTTOK to parse the input record. 
C  ROOT         C   The root (name minus extension) of the input file name. 
C  RPTFILE      C   The report (descriptive output) file name. 
C  PMRPTFILE    C   The PM report (descriptive output) file name. 
C  RPTACTION    C   Open action for the report (descriptive output) file. 
C  SECTION      I   Section of input file where label is allowed (header=0) 
C  STDLAB       C   Full unabbreviated uppercase version of input label 
C  TNUM         I   Database output table number. 
C 
C 
C  Notes: 
C 
C 
C 
      USE DATABASE 
      IMPLICIT NONE 
      CHARACTER (LEN=*),INTENT(IN)   ::  CMDACTION 
      CHARACTER (LEN=*),INTENT(IN)   ::  CMDOUTPUT 
      CHARACTER (LEN=*),INTENT(IN)   ::  INPFILE 
      INTEGER,INTENT(OUT)            ::  RC 
      INTEGER,INTENT(INOUT)          ::  INERR 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'PRINTC.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'FLAGS5.I' 
C 
C     Declare external function. 
      INTEGER, EXTERNAL :: JUNIT 
C 
      CHARACTER (LEN=19)             ::  LABEL 
      CHARACTER (LEN=20)             ::  EXTENSION 
      CHARACTER (LEN=20)             ::  RPTACTION 
      CHARACTER (LEN=40)             ::  STDLAB 
      CHARACTER (LEN=80)             ::  FIELD1 
      CHARACTER (LEN=80)             ::  FIELD2 
      CHARACTER (LEN=80)             ::  ROOT 
      CHARACTER (LEN=140)            ::  LINE1 
      CHARACTER (LEN=140)            ::  LINE2 
      CHARACTER (LEN=150)            ::  ATRPTFILE 
      CHARACTER (LEN=150)            ::  SSFILE 
      CHARACTER (LEN=150)            ::  DATAFIELD 
      CHARACTER (LEN=150)            ::  INPLINE 
      CHARACTER (LEN=150)            ::  PMRPTFILE 
      CHARACTER (LEN=150)            ::  RPTFILE 
      INTEGER                        ::  DBTAB 
      INTEGER                        ::  DBTABMAX 
      INTEGER                        ::  I 
      INTEGER                        ::  IDB 
      INTEGER                        ::  IDX 
      INTEGER                        ::  IOS 
      INTEGER                        ::  LABTYPE 
      INTEGER                        ::  PTR1 
      INTEGER                        ::  PTR2 
      INTEGER                        ::  SECTION 
      INTEGER                        ::  TNUM 
      INTEGER, SAVE                  ::  INO = 0 
      INTEGER, SAVE                  ::  IYES = 1 
      INTEGER, SAVE                  ::  LOWER_LMT = 1 
      INTEGER, SAVE                  ::  UPPER_LMT = 2 
      LOGICAL                        ::  PREFIX 
      LOGICAL, PARAMETER             ::  BACK = .TRUE. 
C 
      RC = 0 
      DBTABMAX = 1   !Only Table 1 (out of 7) will be printed in this MOBILE6 version. 
C 
C     Look for the batch/run mode header on the first line of the file. 
C 
      CALL NXTREC(IOUGEN, IOUOUT, INPLINE, INERR, RC) 
      IF (RC /= 0) GOTO 99 
C 
      CALL UCCOMP(INPLINE(1:19),LABEL) 
      IF (LABEL == 'MOBILE6 BATCH FILE') THEN 
         RC = 1 
         WRITE(IOUOUT,110) 
         WRITE(IOUERR,110) 
  110    FORMAT(/' *** ERROR: The first record in the input file', 
     &           ' contains a BATCH FILE label.', 
     &          /' Batch files can not be listed inside of other', 
     &           ' batch files in Mobile6.') 
         INERR=INERR+1 
         GOTO 99 
      ENDIF 
      IF (LABEL /= 'MOBILE6 INPUT FILE') THEN 
         RC = 2 
         WRITE(IOUOUT,112) 
         WRITE(IOUERR,112) 
  112    FORMAT(/' *** ERROR: The first record in the input file', 
     &           ' does not contain the required', 
     &          /' *** INPUT FILE label. The characters MOBILE6', 
     &           ' INPUT FILE must appear, in', 
     &          /' *** upper or lower case, in the label field of', 
     &           ' the first data record in', 
     &          /' *** every Mobile6 input file.') 
         INERR=INERR+1 
         GOTO 99 
      ENDIF 
C 
C     Find the input file name. It will be used to build default names. 
C 
      IDX = INDEX(INPFILE,'.',BACK) 
      IF (IDX == 0) THEN 
         ROOT = INPFILE 
         EXTENSION = ' ' 
      ELSE 
         ROOT = INPFILE(1:IDX-1) 
         EXTENSION = INPFILE(IDX:) 
         IF (EXTENSION == '.TXT' .OR. EXTENSION == '.PRN' .OR. 
     &       EXTENSION == '.CSV' .OR. EXTENSION(1:3) == '.TB')  THEN 
            WRITE(IOUOUT,115) 
            WRITE(IOUERR,115) 
  115       FORMAT(/' *** ERROR: An input file was found that had', 
     &              ' an extension of "TXT" or "TB?"', 
     &             /' *** or "PRN" or "CSV." These extensions are', 
     &              ' not allowed for input', 
     &             /' *** files in Mobile6. Please rename the file', 
     &              ' and repeat the run.') 
            RC = 11 
            INERR=INERR+1 
            GOTO 99 
         ENDIF 
      ENDIF 
C 
C     Initialize the database file name arrays and initialize the names 
C     for the descriptive and spreadsheet output files. Initialize the 
C     database output flag to 1 (no database output requested.) 
C 
      EXTENSION = ' ' 
      DO DBTAB=1,DBTABMAX 
         WRITE(EXTENSION,"('.TB',I1)") DBTAB 
         DBFNAMES(DBTAB) = ROOT(1:LEN_TRIM(ROOT)) // EXTENSION(1:4) 
         DBACTION(DBTAB) = 'REPLACE' 
      END DO 
C 
      RPTFILE = ROOT 
      RPTACTION = ' ' 
      WMOUT = ROOT 
      DBFLAG = 1 
      DBHEAD = .FALSE. 
      DBDAILY= .FALSE. 
      DBAGGR = .FALSE. 
C 
C     Initialize the selection arrays for emission type, 
C     facility, pollutant and vehicle class. 
C 
      DBSELEFT = 2 
      DBSELFAC = 2 
      DBSELPOL = 2 
      DBSELPOL(4) = 1 ! CO2 is not selected by default 
      DBSELVEH = 2 
      DBSELAGE(LOWER_LMT) = 0 
      DBSELAGE(UPPER_LMT) = 24 
      DB_MDLYR = INO 
C 
C     Initialize "No Descriptive Output" to false and include evap to "on". 
C     Also set the received database emission flag to off. 
C 
      NODESC=.FALSE. 
      EVAP_FLAG=1 
      IDB=0 
C 
C     Initialize Particulate and Air Toxic flags. 
C 
      PARTFLG=.FALSE. 
      TOXFLG=.FALSE. 
      SSFLAG=.FALSE. 
C 
C     Read and interpret the remaining lines in the header. 
C 
      DO  ! Forever 
C 
         CALL NXTREC(IOUGEN, IOUOUT, INPLINE, INERR, RC) 
         IF (RC < 0) GOTO 90 
         IF (RC /= 0) GOTO 99 
C 
         CALL UCCOMP(INPLINE(1:19), LABEL) 
ccs         CALL UCCOMP(INPLINE(21:), DATAFIELD) 
         CALL LJCOMP(INPLINE(21:), DATAFIELD) 
C 
         CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
         IF (RC /= 0) THEN 
           INERR=INERR+1 
           GOTO 99 
         END IF 
C 
         PTR2 = 0 
         CALL NXTTOK(DATAFIELD,FIELD1,PTR1,PTR2) 
         CALL NXTTOK(DATAFIELD,FIELD2,PTR1,PTR2) 
C 
         DBTAB = 0 
         SELECT CASE (STDLAB) 
            CASE ('REPORT FILE', 
     &            'RPT FILE') 
               RPTFILE = FIELD1 
               RPTACTION = FIELD2 
            CASE ('DATABASE OUTPUT', 
     &            'DBASE OUTPUT') 
               DBFLAG = 2 
            CASE ('NO DESC OUTPUT', 
     &            'NO DESC OUT') 
               NODESC=.TRUE. 
               WRITE(IOUOUT,142) 
  142          FORMAT( 
     &         /' WARNING: User has disabled the DESCRIPTIVE OUTPUT', 
     &         /'          for this entire run file.',/) 
            CASE ('NO EVAPORATIVE') 
               EVAP_FLAG = 0 
               WRITE(IOUOUT,144) 
  144          FORMAT( 
     &         /' WARNING: User has disabled the EVAPORATIVE emissions', 
     &         /'          calculations for this entire run file.',/) 
            CASE ('DATABASE OPTIONS', 
     &            'DBASE OPTIONS') 
               CALL DBOPTS(INPLINE,INERR,RC) 
            CASE ('DATABASE EMISSIONS', 
     &            'DBASE EMISSIONS') 
               CALL DBSELE(INPLINE,INERR,RC) 
               IDB=1 
            CASE ('DATABASE FACILITIES', 
     &            'DBASE FACILITIES') 
               CALL DBSELF(INPLINE,INERR,RC) 
            CASE ('POLLUTANTS', 
     &            'POL') 
               CALL DBSELP(INPLINE,INERR,RC) 
            CASE ('DATABASE VEHICLES', 
     &            'DBASE VEHICLES') 
               CALL DBSELV(INPLINE,INERR,RC) 
            CASE ('WITH FIELDNAMES') 
               DBHEAD = .TRUE. 
            CASE ('DAILY OUTPUT') 
               DBDAILY = .TRUE. 
               DBAGGR  = .FALSE. 
            CASE ('AGGREGATED OUTPUT') 
               DBDAILY = .FALSE. 
               DBAGGR  = .TRUE. 
            CASE ('DATABASE YEARS') 
               DB_MDLYR=IYES 
               CALL DBSELA(INPLINE,INERR,RC) 
            CASE ('DATABASE AGES') 
               DB_MDLYR=INO 
               CALL DBSELA(INPLINE,INERR,RC) 
            CASE ('DATABASE HOURS') 
               CALL DBSELH(INPLINE,INERR,RC) 
            CASE ('EMISSIONS TABLE') 
               DBTAB = 1 
            CASE ('INPUTS TABLE') 
               DBTAB = 2 
            CASE ('DATAFILES TABLE') 
               DBTAB = 3 
            CASE ('MESSAGES TABLE') 
               DBTAB = 4 
            CASE ('RUNDATA TABLE') 
               DBTAB = 5 
            CASE ('SCENARIOS TABLE') 
               DBTAB = 6 
            CASE ('ETYPES TABLE') 
               DBTAB = 7 
            CASE ('PARTICULATES') 
               CALL DBSELPRT(INPLINE,INERR,RC) 
            CASE ('AIR TOXICS') 
               CALL DBSELTX(INPLINE,INERR,RC) 
               CALL ATALLOC 
            CASE ('SPREADSHEET') 
               SSFLAG = .TRUE. 
            CASE ('RUN DATA')   ! End of the header section. 
               EXIT 
            CASE DEFAULT 
               IF (SECTION == 0) THEN 
                  WRITE(IOUOUT,120) INPLINE 
                  WRITE(IOUERR,120) INPLINE 
  120             FORMAT (/' *** The label on the following input', 
     &                     ' record is not recognized:',/' *** ',A) 
               ELSE 
                  WRITE(IOUOUT,121) INPLINE 
                  WRITE(IOUERR,121) INPLINE 
  121             FORMAT (/' *** The following record is not', 
     &                     ' allowed in the header of an input', 
     &                     ' file:', 
     &                    /' *** ',A) 
               ENDIF 
               RC = 3 
               INERR=INERR+1 
               GOTO 99 
         END SELECT 
C 
         IF (DBTAB /= 0) THEN 
            IF (INDEX(FIELD1,'.',BACK) /= 0) THEN 
               DBFNAMES(DBTAB) = FIELD1 
            ELSE 
               WRITE(EXTENSION,"('.TB',I1)") DBTAB 
               DBFNAMES(DBTAB) = 
     &                FIELD1(1:LEN_TRIM(FIELD1)) // EXTENSION(1:4) 
            ENDIF 
            DBACTION(DBTAB) = FIELD2 
         ENDIF 
C 
         IF (RC /= 0) THEN 
            WRITE(IOUOUT,140) INPLINE(1:LEN_TRIM(INPLINE)) 
            WRITE(IOUERR,140) INPLINE(1:LEN_TRIM(INPLINE)) 
  140       FORMAT(' *** The following input record caused the', 
     &             ' error.',/' *** ',A) 
            INERR=INERR+1 
            GOTO 99 
         ENDIF 
C 
      END DO 
C 
C     Write runtime warning about NO EVAP disabling the user-selected 
C     or default database output, then disable it. 
C 
      IF((DBSELEFT(3) .EQ. 2  .OR. 
     &    DBSELEFT(4) .EQ. 2  .OR. 
     &    DBSELEFT(5) .EQ. 2  .OR. 
     &    DBSELEFT(6) .EQ. 2  .OR. 
     &    DBSELEFT(7) .EQ. 2  .OR. 
     &    DBSELEFT(8) .EQ. 2) .AND.  
     &    EVAP_FLAG   .EQ. 0) THEN 
        IF(IDB.EQ.1.AND.DBFLAG.EQ.2) WRITE(IOUOUT,145) 
  145   FORMAT(/ 
     &    ' WARNING: The NO EVAPORATIVE command has disabled the'/ 
     &    ' output of evaporative emissions in the database file.'/) 
        DO IDB=3,8 
          DBSELEFT(IDB)=1 
        END DO 
      ENDIF   
C 
C     Use the 'command line' inputs, if they were given, for all of the 
C     output file names - in preference to any any names specified in the 
C     input file header or to any default name generated from the input 
C     file name. 
C 
ccs      IF (CMDOUTPUT /= ' ') THEN 
C 
ccs         IDX = INDEX(CMDOUTPUT,'.',BACK) 
ccs         IF (IDX == 0) THEN 
ccs            ROOT = CMDOUTPUT 
ccs         ELSE 
ccs            ROOT = CMDOUTPUT(1:IDX-1) 
ccs         ENDIF 
C 
ccs         CALL UCCOMP(ROOT,RPTFILE) 
ccs         CALL UCCOMP(CMDACTION,RPTACTION) 
C 
ccs         WMOUT = ROOT 
C 
ccs         EXTENSION = ' ' 
ccs         DO DBTAB=1,DBTABMAX 
ccs            WRITE(EXTENSION,"('.TB',I1)") DBTAB 
ccs            DBFNAMES(DBTAB) = ROOT(1:LEN_TRIM(ROOT))//EXTENSION(1:4) 
ccs            DBACTION(DBTAB) = RPTACTION 
ccs         END DO 
C 
ccs      ENDIF 
C 
C     Check the 'open' actions for the report file (descriptive output 
C     file), the spreadsheet output file, and the database tables. 
C 
      IF (.NOT. PREFIX(RPTACTION,'APPEND ')  .AND. 
     &    .NOT. PREFIX(RPTACTION,'REPLACE ') .AND. 
     &    RPTACTION /= '') THEN 
         WRITE(IOUOUT,150) RPTACTION 
         WRITE(IOUERR,150) RPTACTION 
  150    FORMAT(/' *** ERROR: "',A7,'" is not a valid open action', 
     &           ' for the descriptive', 
     &          /' *** output (report) file. The action must be', 
     &           ' APPEND or REPLACE.') 
         INERR=INERR+1 
         RC = 5 
      ENDIF 
C 
      DO DBTAB = 1,DBTABMAX 
         IF (.NOT. PREFIX(DBACTION(DBTAB),'APPEND ')  .AND. 
     &       .NOT. PREFIX(DBACTION(DBTAB),'REPLACE ') .AND. 
     &       DBACTION(DBTAB) /= '') THEN 
            WRITE(IOUOUT,160) DBACTION(DBTAB),DBTABLES(DBTAB) 
            WRITE(IOUERR,160) DBACTION(DBTAB),DBTABLES(DBTAB) 
  160       FORMAT(/' *** ERROR: "',A7,'" is not a valid open', 
     &              ' action for the ',A10, 
     &             /' *** database table. The action must be', 
     &              ' APPEND or REPLACE.') 
            INERR=INERR+1 
            RC = 5 
         ENDIF 
      END DO 
C 
      IF (RC /= 0) GOTO 99 
C 
C     If a descriptive output file name was specified on the command line 
C     (as part of the response to the prompt for an input file name) and 
C     we are in run mode (and not in batch mode), use that name for the 
C     ouptut file. If no output file name was specified on the command 
C     line, use the output file name specified in the header of the input 
C     file (assigned in the loop above). Otherwise, generate a default 
C     output file name using the root of the input file name (the file 
C     name with any extension removed - also assigned above). If the 
C     descriptive output file name contains a period (it must be a user- 
C     supplied name), use it exactly as given. If it does not (it may be 
C     a user-supplied name or a default name), add the TXT extension to 
C     the name. 
C 
      IF (INDEX(RPTFILE,'.',BACK) == 0) THEN 
         RPTFILE = RPTFILE(1:LEN_TRIM(RPTFILE)) // '.TXT' 
      ENDIF 
C 
      IOUREP = JUNIT() 
       
      IF (PREFIX(RPTACTION,'APPEND ')) THEN 
         OPEN(IOUREP,FILE=RPTFILE,POSITION='APPEND',ERR=80,IOSTAT=IOS) 
      ELSE 
         OPEN(IOUREP,FILE=RPTFILE,STATUS='REPLACE',ERR=80,IOSTAT=IOS) 
      ENDIF 
C 
      WRITE(IOUOUT,165) RPTFILE(1:LEN_TRIM(RPTFILE)) 
  165 FORMAT(' * Report file: ',A) 
C 
C     Repeat file actions for PM output, if necessary 
C 
      IF(PARTFLG) THEN 
        IDX = INDEX(RPTFILE,'.',BACK) 
        IF (IDX == 0) THEN 
           PMRPTFILE = RPTFILE(1:LEN_TRIM(RPTFILE)) // '.PM' 
        ELSE 
           PMRPTFILE = RPTFILE(1:LEN_TRIM(RPTFILE(1:IDX-1))) // '.PM' 
        ENDIF 
C 
        IOUPM = JUNIT() 
        IF (PREFIX(RPTACTION,'APPEND ')) THEN 
          OPEN(IOUPM,FILE=PMRPTFILE,POSITION='APPEND',ERR=80,IOSTAT=IOS) 
        ELSE 
          OPEN(IOUPM,FILE=PMRPTFILE,STATUS='REPLACE',ERR=80,IOSTAT=IOS) 
        ENDIF 
C 
        WRITE(IOUOUT,166) PMRPTFILE(1:LEN_TRIM(PMRPTFILE)) 
  166   FORMAT(' * PM Report file: ',A) 
      ENDIF 
C 
C     Repeat file actions for Air Toxic (AT) output, if necessary   
C 
      IF(TOXFLG) THEN 
        IDX = INDEX(RPTFILE,'.',BACK) 
        IF (IDX == 0) THEN 
           ATRPTFILE = RPTFILE(1:LEN_TRIM(RPTFILE)) // '.TOX' 
        ELSE 
           ATRPTFILE = RPTFILE(1:LEN_TRIM(RPTFILE(1:IDX-1))) // '.TOX' 
        ENDIF 
C 
        IOUAT = JUNIT() 
        IF (PREFIX(RPTACTION,'APPEND ')) THEN 
          OPEN(IOUAT,FILE=ATRPTFILE,POSITION='APPEND',ERR=80,IOSTAT=IOS) 
        ELSE 
          OPEN(IOUAT,FILE=ATRPTFILE,STATUS='REPLACE',ERR=80,IOSTAT=IOS) 
        ENDIF 
C 
        WRITE(IOUOUT,167) ATRPTFILE(1:LEN_TRIM(ATRPTFILE)) 
  167   FORMAT(' * AT Report file: ',A) 
      ENDIF 
C 
C     Repeat file actions for spreadsheet output, if necessary 
C 
      IF(SSFLAG) THEN 
        IDX = INDEX(RPTFILE,'.',BACK) 
        IF (IDX == 0) THEN 
           SSFILE = RPTFILE(1:LEN_TRIM(RPTFILE)) // '.TAB' 
        ELSE 
           SSFILE = RPTFILE(1:LEN_TRIM(RPTFILE(1:IDX-1))) // '.TAB' 
        ENDIF 
C 
        IOUSS = JUNIT() 
        IF (PREFIX(RPTACTION,'APPEND ')) THEN 
          OPEN(IOUSS,FILE=SSFILE,POSITION='APPEND',ERR=80, 
     *         IOSTAT=IOS) 
        ELSE 
          OPEN(IOUSS,FILE=SSFILE,STATUS='REPLACE',ERR=80, 
     *         IOSTAT=IOS) 
        ENDIF 
C 
        WRITE(IOUOUT,168) SSFILE(1:LEN_TRIM(SSFILE)) 
  168   FORMAT(' * Spreadsheet file: ',A) 
      ENDIF 
C 
C     Write warning about NO EVAP in the report file. 
C 
      IF(EVAP_FLAG.EQ.0) WRITE(IOUREP,144) 
C 
C     If database output is requested, open the database files, 
C     write record containing field names if requested, and 
C     write the ETYPES table. (This table is a constant - it does 
C     not depend on any user input.) 
C 
      IF (DBFLAG == 2) THEN 
         DO DBTAB = 1,DBTABMAX 
            TNUM = DBTAB 
            DBUNIT(DBTAB) = JUNIT() 
            IF (PREFIX(DBACTION(DBTAB),'APPEND ')) THEN 
               OPEN(DBUNIT(DBTAB),FILE=DBFNAMES(DBTAB), 
     &           POSITION='APPEND',ERR=85,IOSTAT=IOS) 
            ELSE IF (PREFIX(DBACTION(DBTAB),'REPLACE ') .OR. 
     &               DBACTION(DBTAB) == '') THEN 
               OPEN(DBUNIT(DBTAB),FILE=DBFNAMES(DBTAB), 
     &           STATUS='REPLACE',ERR=85,IOSTAT=IOS) 
            ENDIF 
         END DO 
C 
C     WRITE FIELDNAMES RECORD IF REQUESTED 
C 
      IF (DBHEAD) THEN 
        IF (DBDAILY) THEN 
                WRITE(DBUNIT(1),180) 
     &          (DBDAYCOLS(I), TABCHAR, I=1,17), DBDAYCOLS(18) 
180             FORMAT(17(A8,A1),A8) 
        ELSEIF(DBAGGR) THEN 
                WRITE(DBUNIT(1),185) 
     &          (DBAGGRCOLS(I), TABCHAR, I=1,12), DBAGGRCOLS(13) 
185             FORMAT(12(A8,A1),A8) 
        ELSE 
                WRITE(DBUNIT(1),190) 
     &          (DBHOURCOLS(I), TABCHAR, I=1,21), DBHOURCOLS(22) 
190             FORMAT(21(A8,A1),A8) 
        END IF 
 
      END IF 
C 
C        Produce the EFTYPES table. This table is a constant that does 
C        not depend on user input.  This option has been disabled for the 
C        time being.  Database output tables 2 through 7 do not contain 
C        meaningful or fully tested results. 
C 
C        WRITE(DBUNIT(7),170) (IDX,TABCHAR,DBEFTNAMES(IDX),IDX=1,8) 
C 170    FORMAT(I2,A1,A16) 
C        CLOSE(DBUNIT(7)) 
C 
      END IF 
C 
      GOTO 99 
C 
   80 WRITE(LINE1,280) IOS,RPTFILE 
  280 FORMAT(' *** ERROR ',I9,' opening the descriptive output', 
     &       ' file ',A80) 
      CALL LJCOMP(LINE1,LINE2) 
      WRITE(IOUOUT,'(" "/A)') LINE2 
      WRITE(IOUERR,'(" "/A)') LINE2 
      RC = 7 
      INERR=INERR+1 
      GOTO 99 
C 
   85 WRITE(LINE1,285) IOS,DBFNAMES(TNUM) 
  285 FORMAT(' *** ERROR ',I9,' opening database output file ',A80) 
      CALL LJCOMP(LINE1,LINE2) 
      WRITE(IOUOUT,'(" "/A)') LINE2 
      WRITE(IOUERR,'(" "/A)') LINE2 
      INERR=INERR+1 
      RC = 8 
      GOTO 99 
C 
   90 WRITE(IOUOUT,290) INPFILE 
      WRITE(IOUERR,290) INPFILE 
  290 FORMAT(/' *** ERROR: Unexpected end of file before header', 
     &        ' processing is complete', 
     &       /'*** for file ',A80) 
C 
   99 RETURN 
      END SUBROUTINE PROCHDR 
