      SUBROUTINE GETTAM(INPREC,INERR,RC) 
C 
C  GETTAM reads in, validates and stores the replacement tampering  
C  rates that were read from an external data file. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, UCCOMP, NXTREC, QUITER. 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  11 Apr 00 AIR Task 02: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C   2 Mar 99 @DynTel-MLA 2-659  Cleaned up comments and error messages. 
C  29 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS2/ IMFLAG 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS1/ TAMFLG 
C    /TAMEQ1/ TAMZML,TAMDR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  EXTDATA  L    Is external data file open? 
C  FILENAME C    Name of external data file 
C  ID       I    Tampering type  
C  IM       I    Are rates for I/M or non-I/M area? 
C  IGD      I    Model year group index for tampering rates 
C  IVTAM    I    Vehicle type 
C  LABEL    C    Label field from a labeled input record 
C  LABTYPE  I    Numeric code for an input label 
C  LASTREC  C    The previous record read from the input file 
C  M6REC    C    Record from the external data file 
C  NTAM     I    Number of sets of tampering rates to read 
C  SECTION  I    Dummy argument for CHKLAB. 
C  SOURCE   I    Logical unit number for the input file 
C  STDLAB   C    Dummy argument for CHKLAB. 
C  TAMPRATE I    Numeric code for the TAMPERING RATES input label 
C  ZEROML   R    user entered tampering zml to be checked: ZEROML <= 1.0 
C 
C  Notes: 
C 
C  11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C            input formats. Additional records are required for new 
C            vehicle types. 
C 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C 
C  27 Mar 98 @DynTel-MLA 2-659  Updated for Mobile6 (the transitional 
C            labeled inputs version). This subroutine now sets TAMFLG. 
C            The call to GETTAM has been moved to DSPTCH. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'TAMEQ1.I' 
C 
      INTEGER, EXTERNAL :: JUNIT   
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     IM, IGD, IVTAM, ID 
      REAL        ZEROML 
      LOGICAL     EXTDATA 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
      INTEGER     NTAM 
      INTEGER     SECTION 
      INTEGER     SOURCE 
      CHARACTER   STDLAB*40 
      INTEGER     TAMPRATE 
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=50) M6DIR
      CHARACTER(LEN=80) M6FILE
C 
      RC = 0 
      IF (TAMFLG.EQ.1) THEN 
        TAMFLG = 2 
      ELSE 
        GOTO 96 
      ENDIF 
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
        CALL CHKLAB('TAMPER RATES',STDLAB,TAMPRATE,SECTION) 
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
     &       IOSTAT=RC,ERR=80) 
        WRITE(IOUREP,120) FILENAME 
  120   FORMAT(/'* Reading non-default TAMPERING RATES from the ', 
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
        IF (LABTYPE.NE.TAMPRATE) THEN 
          GOTO 90 
        ENDIF 
C 
      ENDIF 
C 
C     Read the data field from the input record. 
C 
      READ(M6REC,110,ERR=92,END=92,IOSTAT=RC) NTAM 
  110 FORMAT(T22,BN,I3) 
      IF (NTAM.NE.1.AND.NTAM.NE.2) GOTO 95 
C 
C     Find the first data line. 
C 
      CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
      IF (RC.NE.0) GOTO 98 
      BACKSPACE(SOURCE) 
C 
C     Read the tampering rates. 
C 
      IF (NTAM.EQ.2) THEN 
        GOTO 10 
      ENDIF 
C 
      READ(SOURCE,100,ERR=98,END=98,IOSTAT=RC) 
     *   (((TAMZML(ID,IVTAM,IGD,1),ID=1,8),IVTAM=1,14),IGD=1,3) 
  100 FORMAT(8F8.4) 
      READ(SOURCE,100,ERR=98,END=98,IOSTAT=RC) 
     *   (((TAMDR(ID,IVTAM,IGD,1),ID=1,8),IVTAM=1,14),IGD=1,3) 
      GOTO 20 
C 
C  Read in all but cell 9 of TAMZML & TAMDR (if picked up cell 9 too, could 
C  just put in READ(---) TAMZML / TAMDR, since user inputs in FTN array order. 
C 
   10 CONTINUE 
      READ(SOURCE,100,ERR=98,END=98,IOSTAT=RC) 
     *  ((((TAMZML(ID,IVTAM,IGD,IM),ID=1,8),IVTAM=1,14),IGD=1,3),IM=1,2) 
      READ(SOURCE,100,ERR=98,END=98,IOSTAT=RC) 
     *  ((((TAMDR(ID,IVTAM,IGD,IM),ID=1,8),IVTAM=1,14),IGD=1,3),IM=1,2) 
C 
      IF (EXTDATA) CLOSE(SOURCE) 
C 
C  Move overall misfueling ZML/DR to ID=9.  Compute other misfueling ZML/DR 
C  for ID=4 as (overall - filler neck).  For every disablement case, 
C  insure tampering zero mile levels (intercepts) do not exceed unity. 
C  There is no range check for deterioration rates. 
C 
   20 DO IM=1,NTAM 
   25   DO IGD=1,3 
          DO IVTAM=1,14 
            TAMZML(9,IVTAM,IGD,IM)=TAMZML(4,IVTAM,IGD,IM) 
            TAMDR(9,IVTAM,IGD,IM)=TAMDR(4,IVTAM,IGD,IM) 
            TAMZML(4,IVTAM,IGD,IM)= 
     *        TAMZML(9,IVTAM,IGD,IM)-TAMZML(3,IVTAM,IGD,IM) 
            TAMDR(4,IVTAM,IGD,IM)=TAMDR(9,IVTAM,IGD,IM)- 
     *        TAMDR(3,IVTAM,IGD,IM) 
            DO ID=1,9 
              ZEROML=TAMZML(ID,IVTAM,IGD,IM) 
              IF(ZEROML.GT.1.0) CALL QUITER(ZEROML,0,50,INERR) 
            END DO 
          END DO 
        END DO 
      END DO 
C 
   40 GOTO 99 
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
     &            'the following external data file:',/'***',A) 
          RC = 1 
          INERR = INERR + 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
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
      RC = 2   
      INERR = INERR + 1 
      GOTO 99 
C 
   92 WRITE(IOUERR,192) M6REC 
  192 FORMAT(/'*** ERROR reading the data field of the following ', 
     &        'labeled input record:', 
     &       /'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
   95 WRITE(IOUERR,195) M6REC 
  195 FORMAT('*** ERROR: The data field of the following record ', 
     &       'must contain either 1 or 2:',/'*** ',A) 
      RC = 3 
      INERR = INERR + 1 
      GOTO 99 
C 
   96 WRITE(IOUERR,196) 
  196 FORMAT(/'*** ERROR: Alternate tampering rates ', 
     &        'have already been read in this run.', 
     &       /'*** Only one set of tampering rates ', 
     &        'can be entered in each run.'/) 
      RC = 3   
      INERR = INERR + 1 
      GOTO 99 
C 
   98 CALL QUITER(0.,0,126,INERR) 
      IF (RC.LT.0) THEN 
        WRITE(IOUERR,198) 
  198   FORMAT('*** (Unexpected end of file.)') 
      ELSE 
        BACKSPACE(SOURCE) 
        READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199   FORMAT(A) 
        WRITE(IOUERR,200) LASTREC 
  200   FORMAT('*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
C 
   99 RETURN 
      END 
