      SUBROUTINE GETT2CRT(INPREC,INERR,RC) 
C 
C 
C  GET_T2_CERT reads user supplied Tier2 50K certification standards. 
C  There are 12 bins for each vehicle type, pollutant, model year 
C  combination.  Values read in will replace the default values 
C  located in BD_TIER2. 
C 
C 
C  Called by: DSPTCH. 
C 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  23 Apr 01 AIR Task 11: Redirected error messages to IOUERR. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  15 Mar 00 @DynTel-ddj 1-032  Implemented echo print of TIER2_50K 
C            input file. 
C   1 Feb 00 @DynTel-ddj 1-021, Echo Print user supplied TIER2_50K 
C  26 Jan 00 @DynTel-ddj 1-021, Echo Print user supplied TIER2_50K 
C  30 Dec 99 @DynTel-ddj 1-021, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C    /IOUCOM/    IOUALT, IOUOUT, IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C    /TIER2/     TIER2_50K 
C 
C    common block array subscripts: 
C    TIER2_50K     TIER2_50K(IV, IBIN, IP) 
C 
C 
C    common variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  ------------------------------------------------------- 
C TIER2_50K    R    Array containing Tier2 50K certification standards. 
C 
C 
C 
C   Local array subscripts: 
C    DUMMY         DUMMY(IV, IBIN, IP) 
C 
C 
C   Local variable / array dictionary: 
C 
C   Name    Type              Description 
C  ------   ----  ------------------------------------------------------ 
C   DUMMY    R    Array for temporary storage of user supplied values 
C                 for Tier2 50K Certification Standards prior to 
C                 validation of user supplied values. 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'TIER2.I' 
C  
      INTEGER, EXTERNAL :: JUNIT    
C 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     SOURCE 
      INTEGER     I 
      INTEGER     IBIN 
      INTEGER     IP 
      INTEGER     IV 
      INTEGER     LABTYPE 
      INTEGER     LDGIMP 
      INTEGER     NCHAR 
      INTEGER     SECTION 
C 
      REAL DUMMY(5,12,3) 
C 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*40 
      CHARACTER   STDLAB*40 
      CHARACTER   INPREC*(*) 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=50) M6DIR 
      CHARACTER(LEN=80) M6FILE 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('T2 CERT',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Find out whether we are reading from an external data file. 
C 
      READ(INPREC(21:),*,IOSTAT=RC,END = 10,ERR=97) FILENAME 
C 
C       Open the external data file. 
C 
  10  IF ( RC.EQ.0 ) THEN 
         SOURCE = JUNIT() 
 
ccs      All M6 input files live in SMK_M6PATH 
         CALL GETENV( "SMK_M6PATH", M6DIR ) 
         M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME 
         
         OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *        IOSTAT=RC,ERR=80) 
C 
C     Read the required header line from the file. 
C 
         CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
         IF (RC.EQ.0) THEN 
C 
C     Check for a valid labeled input record. 
C 
           CALL UCCOMP(M6REC(1:19),LABEL) 
           CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
           IF (LABTYPE.EQ.LDGIMP) THEN 
C 
             READ(SOURCE,*,ERR=98,END=98,IOSTAT=RC) DUMMY 
             CLOSE(SOURCE) 
             WRITE(IOUREP,300) 
  300        FORMAT(/,'  Reading User Supplied Tier2 50K certification', 
     &          ' standards',/) 
C 
             DO IV = 1, 5 
               DO IP = 1, 3 
                 DO IBIN = 1, 12 
                   IF (DUMMY(IV,IBIN,IP).LT.0.0) THEN 
                     CALL QUITER(DUMMY(IV,IBIN,IP),IP,530,INERR) 
 
C      Unacceptable input for Tier2 50K certification standards 
C      (array TIER2_50K).  Execution returned to DSPTCH without 
C      updating TIER2_50K. 
C 
                     RETURN 
                   END IF 
                END DO 
              END DO 
            END DO 
C 
C      Input values all .GE. zero.  Assign values to TIER2_50K. 
C 
            TIER2_50K = DUMMY 
C 
C      Echo Print FILENAME for user supplied TIER2_50K 
C 
            WRITE(IOUREP,301) FILENAME 
  301       FORMAT('     Data read from file: ',A80,/ ) 
C 
          ELSE   !  LABTYPE.NE.LDGIMP 
             WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
             RC = 2 
           END IF 
         ELSE 
           SELECT CASE (RC) 
           CASE (:-1) 
             WRITE(IOUERR,185) TRIM(FILENAME) 
  185        FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
             RC = 1 
           CASE (1:) 
             WRITE(IOUERR,186) TRIM(FILENAME) 
  186        FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
           END SELECT 
         END IF 
C 
      ELSE ! RC not equal to zero 
C 
C     Error opening the external data file. 
C 
   80    WRITE(IOUERR,180) FILENAME 
  180    FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
   97    RC = MOD(RC,256) 
         IF (RC.EQ.71) THEN 
           WRITE(IOUERR,181) 
  181      FORMAT('*** (The file name or path could not be found.)') 
           DO I = LEN(FILENAME),1,-1 
             NCHAR = I 
             IF (FILENAME(I:I).NE.' ') EXIT 
           END DO 
           IF (NCHAR.GT.50) THEN 
             WRITE(IOUERR,182) 
  182        FORMAT('*** (The file name must not be longer than 50 ', 
     &        'characters.)') 
             INERR = INERR + 1 
           END IF 
         ELSE 
   98      IF (RC.NE.0) THEN 
              CALL QUITER(0.,0,177,INERR) 
              IF (RC.LT.0) THEN 
                 WRITE(IOUERR,198) 
  198            FORMAT('*** (Unexpected end of file.)') 
              ELSE 
                 BACKSPACE(SOURCE) 
                 READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199            FORMAT(A) 
                 WRITE(IOUERR,200) LASTREC 
  200            FORMAT('*** The following record caused the error:', 
     &              /'*** ',A) 
              END IF 
           END IF 
         ENDIF 
      END IF 
C 
  99  RETURN 
      END 
