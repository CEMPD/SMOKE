      SUBROUTINE GETNH3SDR(INERR,RC) 
C 
C  GETNH3SDR reads the Ammonia (NH3) PM Sulfur Deterioration Rates 
C  (g/mi/ppm) from the file PMNH3SDR.D. 
C 
C  Changes: (Last change first) 
C 
C  09 Apr 01 AIR Task 11: New Subroutone 
C 
C  Called by: PMEMFAC. 
C 
C  Calls CHKLAB 
C 
C  Input on call: 
C 
C    argument list: INERR 
C 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C 
C  Common block array Subscripts: 
C 
C     REAL NH3SDR(ICAT,IROAD,IV) - 
C 
C  ICAT is for catalyst, as follows: 
C    1 = no catalyst or diesel 
C    2 = oxidation catalyst 
C    3 = 3-way catalyst 
C 
C  IROAD and R are the roadway type, as follows: 
C    1 = freeways 
C    2 = arterial 
C    3 = local 
C    4 = ramp 
C 
C  IV is the MOBILE6 vehicle class 
C 
C  NH3SDR values are given in g/mi/ppm. 
C 
C  Local variable Subscripts: 
C 
C    DUMMY       DUMMY( MY ) 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store 3 catalyst type (NoCat/Dsl,  
C                   OxCat, 3-way) NH3 Sulfur Deterioration Rates (g/mi/ppm). 
C   I          I    Do Loop index 
C   RDFLG      L    Saved flag to prevent multiple reads of PMNH3SDR.D 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART6.I' 
C  
      INTEGER, EXTERNAL :: JUNIT
C
      INTEGER       INERR 
      INTEGER       RC 
C 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       ICAT 
      INTEGER       IROAD 
      INTEGER       IV 
      INTEGER       LABTYPE 
      INTEGER       LDGIMP 
      INTEGER       SECTION 
      LOGICAL,SAVE::RDFLG=.FALSE. 
      REAL          DUMMY(3) 
       
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     If data has already been read in, return 
C 
      IF(RDFLG) RETURN 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('PMNH3SDR',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file PMNH3SDR.D. 
C 
      IOUALT = JUNIT()

      OPEN(IOUALT,FILE='PMNH3SDR.D',STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120)  
  120 FORMAT(/'* Reading Ammonia (NH3) Sulfur Deterioration Rates '/ 
     &        '* from the external data file PMNH3SDR.D') 
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
C     Read the file containing the ammonia sulfur deterioration rates (g/mi/ppm) 
C     Read data a line a time, skipping any blank lines or comments. 
C 
C       Find the first data line in the next data block. 
C 
      DO I = 1,32767 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C       Read and validate the basic emission rates. 
C 
        READ (IOUALT,*,IOSTAT=RC) IROAD,IV,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
        IF (IV.GT.MAXVEH.OR.IV.LT.1 .OR. 
     *      IROAD.GT.4.OR.IROAD.LT.1) GOTO 92 
C 
        DO ICAT=1,3 
          NH3SDR(ICAT,IROAD,IV) = DUMMY(ICAT) 
        END DO 
C 
      END DO 
C 
      RC = 0 
      RDFLG=.TRUE. 
C 
      GOTO 99 
C 
C     Close the datafile and write the appropriate warning message. 
C 
  80  WRITE(IOUOUT,180,ERR=98)  
      WRITE(IOUERR,180,ERR=98)  
 180  FORMAT('*** Error opening external data file PMNH3SDR.D') 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185)  
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            'the external data file:',/'*** PMNH3SDR.D') 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186)  
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            'external data file:',/'*** PMNH3SDR.D') 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in'/ 
     &        'the header of the external data file PMNH3SDR.D:'/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   92 WRITE(IOUERR,192) TRIM(M6REC) 
  192 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        'indexes or data in the external data file PMNH3SDR.D:'/ 
     &        1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
