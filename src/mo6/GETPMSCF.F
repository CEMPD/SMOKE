      SUBROUTINE GETPMSCF(INERR,RC) 
C 
C  GETPMSCF reads the EMFAC2000v2.02-based PM Speed Correction 
C  Factor Coefficients from the file PMSCFC.D. 
C 
C  Changes: (Last change first) 
C 
C  27 Mar 01 AIR Task 11: New Subroutone 
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
 
C  PMSCFC(46,4,3,2,2,28) - PMSCFC (MY,ICOEF,IFDS,ICAT,MODE,IV) 
C 
C  MY is Model Year Group, 1=1975-, 46=2020+ 
C 
C  ICOEF are for SCF equation coefficient, as follows 
C    A, B, C, D 
C 
C  IFDS is for fuel delivery system, as follows 
C    1 = Carb/Diesel 
C    2 = TBI 
C    3 = MPFI 
C 
C  ICAT is for catalyst, as follows: 
C    1 = no catalyst 
C    2 = catalyst 
C 
C  MODE is for cyle mode, as follows: 
C    1 = FTP 
C    2 = Unified Cycle 
C 
C  IV is the MOBILE6 Vehicle Class 
C 
C  Local variable Subscripts: 
C 
C    DUMMY       DUMMY( MY ) 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C   DUMMY      R    Temporary array to store 46 model years (1975- through 
C                   2020+ of the EMFAC2000v2.02-based Speed Correction Factor 
C                   Equation Coefficients. 
C   I          I    Do Loop index 
C   RDFLG      L    Saved flag to prevent multiple reads of PMSCFC.D 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART5.I' 
C 
      INTEGER       INERR 
      INTEGER       RC 
C 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       ICAT 
      INTEGER       ICOEF 
      INTEGER       IFDS 
      INTEGER       IV 
      INTEGER       LABTYPE 
      INTEGER       LDGIMP 
      INTEGER       MODE 
      INTEGER       MY 
      INTEGER       SECTION 
      LOGICAL,SAVE::RDFLG=.FALSE. 
      REAL          DUMMY(46) 
       
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
      CALL CHKLAB('PMSCFC',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
C 
C     Open the external data file PMSCFC.D. 
C 
      OPEN(IOUALT,FILE='PMSCFC.D',STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120)  
  120 FORMAT(/'* Reading PM Speed Correction Factor Equation '/ 
     &        '* Coefficients from the external data file PMSCFC.D') 
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
C     Read the file containing the EMFAC2000v2.02-based Speed Correction 
C     Factor Equation Coefficients.. Read data a line a time, skipping 
C     any blank lines or comments. 
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
C       Read and validate the speed correction factor equation coefficients. 
C 
        READ (IOUALT,*,IOSTAT=RC) ICOEF,IFDS,ICAT,MODE,IV,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
        IF (IV.GT.MAXVEH.OR.IV.LT.1 .OR. 
     *      IFDS.GT.3.OR.IFDS.LT.1  .OR. 
     *      MODE.GT.2.OR.MODE.LT.1  .OR. 
     *      ICAT.GT.2.OR.ICAT.LT.1  .OR. 
     *      ICOEF.GT.4.OR.ICOEF.LT.1) GOTO 92 
         
C 
        DO MY=1,46 
          PMSCFC(MY,ICOEF,IFDS,ICAT,MODE,IV) = DUMMY(MY) 
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
 180  FORMAT('*** Error opening external data file PMSCFC.D') 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185)  
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            'the external data file:',/'*** PMSCFC.D') 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186)  
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            'external data file:',/'*** PMSCFC.D') 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in'/ 
     &        'the header of the external data file PMSCFC.D:'/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   92 WRITE(IOUERR,192) TRIM(M6REC) 
  192 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        'indexes or data in the external data file PMSCFC.D:'/ 
     &        1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
