      SUBROUTINE GETPMGZML(FILENAME,INERR,RC) 
C 
C  GETPMGZML reads the gasoline vehicle carbon ZML levels, g/mi, from the 
C  user-supplied file. 
C 
C  Changes: (Last change first) 
C 
C  15 Oct 01 AIR Task 25: New Subroutone 
C 
C  Called by: GETPMEF. 
C 
C  Calls CHKLAB, UCCOMP 
C 
C  Input on call: 
C 
C    argument list: INERR,FILENAME 
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
C  PMGZML(71,4,28) - PMGZML (MY,IFC,IV) 
C 
C  MY is Model Year Group, 1=1950-, 71=2020+ 
C 
C  IFC is the PART5 Fuel/Catalyst Index, as follows: 
C    1 = Leaded Fuel (regardless of Cat or Air) 
C    2 = Unleaded Fuel, Catalyst, No Air Pump 
C    3 = Unleaded Fuel, Catalyst, Air Pump 
C    4 = Unleaded Fuel, No Catalyst 
C 
C  CAT is for catalyst equipment: Yes, No 
C 
C  FUEL is for fuel type: Leaded, Unleaded 
C 
C  AIR is for air pump: Yes, No 
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
C   DUMMY      R    Temporary array to store 71 model years (1950- through 
C                   2020+ of the PART5-based carbon levels (g/mi). 
C   I          I    Do Loop index 
C   L          I    Lines read in 
C   IFC        I    Fuel/Catalyst index 
C   IDUMMY     I    Temporary input for PMDRAGE 
C   MAPIV      I    Maps the gasoline only input into the PMDRAGE array 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EFGAS.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'PART4.I' 
C 
      INTEGER, EXTERNAL :: JUNIT 
C 
      CHARACTER*(*) FILENAME 
      INTEGER       INERR 
      INTEGER       RC 
C 
      CHARACTER     AIR*3 
      CHARACTER     CAT*3 
      CHARACTER     FUEL*8 
      CHARACTER     INAIR*3 
      CHARACTER     INCAT*3 
      CHARACTER     INFUEL*8 
      CHARACTER     LABEL*40 
      CHARACTER     M6REC*150 
      CHARACTER     STDLAB*40 
      INTEGER       I 
      INTEGER       IDUMMY(15) 
      INTEGER       IFC 
      INTEGER       IV 
      INTEGER       L 
      INTEGER       LABTYPE 
      INTEGER       LDGIMP 
      INTEGER       MY 
      INTEGER       SECTION 
      INTEGER, SAVE, DIMENSION(15) ::  
     *              MAPIV=(/1,2,3,4,5,6,7,8,9,10,11,12,13,24,25/) 
      REAL          DUMMY(71) 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Open the external data file for PM carbon. 
C       
      IOUALT = JUNIT() 
       
      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=80) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading PM Gas Carbon ZML Levels '/ 
     &        '* from the external data file ',A) 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('PMGZML',STDLAB,LDGIMP,SECTION) 
C 
C     The next two assignment statements are present solely to suppress 
C     compiler generated warnings. They have no functional importance. 
C 
      LABEL = STDLAB(1:1) 
      LABTYPE = SECTION 
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
C     Find, read and verify the age when PMDDR1 ends and PMDDR2 begins. 
C     Input is for gasoline vehicles only. 
C 
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC.LT.0) GOTO 91 
      BACKSPACE(IOUALT) 
C 
      READ (IOUALT,*,IOSTAT=RC) IDUMMY 
C 
      DO I=1,15 
      IV=MAPIV(I) 
      PMDRAGE(IV)=IDUMMY(I) 
      IF(PMDRAGE(IV).LT. 0 .OR. 
     *   PMDRAGE(IV).GT.25) GOTO 96 
      END DO 
C       
C     Read the file containing the PM Gas Carbon Levels (g/mi) 
C     Read data a line a time, skipping any blank lines or comments. 
C 
      L = 0  ! Initialize lines read 
C 
      DO I = 1,32767 
C 
        CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
        IF (RC.LT.0) THEN 
          EXIT 
        ENDIF 
        BACKSPACE(IOUALT) 
C 
C  Read and validate the carbon levels. 
C 
        READ (IOUALT,*,IOSTAT=RC) IV,INFUEL,INCAT,INAIR,DUMMY 
C 
        IF(RC.LT.0) EXIT 
C 
C  Verify the indexing input 
C 
        IF (IV.GT.MAXVEH.OR.IV.LT.1) GOTO 95 
C 
        IF(VVGASMC(IV).EQ.0) GOTO 94  ! Data must be for Gasoline Vehicles Only 
C 
        CALL UCCOMP(INFUEL,FUEL) 
        IF(FUEL.NE.'LEADED' .AND. FUEL.NE.'UNLEADED') GOTO 91  
C 
        CALL UCCOMP(INCAT,CAT) 
        IF(CAT.NE.'YES' .AND. CAT.NE.'NO' .AND. CAT.NE.'NA') GOTO 92 
C 
        CALL UCCOMP(INAIR,AIR) 
        IF(AIR.NE.'YES' .AND. AIR.NE.'NO' .AND. AIR.NE.'NA') GOTO 93 
C 
C  Set the PART5 fuel/catalyst index 
C 
        IF(FUEL.EQ.'LEADED') THEN 
          IFC=1  ! Leaded 
        ELSE 
          IF(CAT.EQ.'NO') THEN 
            IFC=4  ! Unleaded, No Cat 
          ELSE 
            IF(AIR.EQ.'NO') THEN 
              IFC=2  ! Unleaded, Cat, No Air 
            ELSE 
              IFC=3  ! Unleaded, Cat, Air 
            END IF 
          END IF 
        END IF 
C 
C  Store the carbon values 
C 
        DO MY=1,71 
          PMGZML(MY,IFC,IV) = DUMMY(MY) 
        END DO 
C 
        L = L + 1  ! Update lines read in 
C 
      END DO 
C 
C     Make sure all carbon lines are read in since EXCEL does not terminate 
C     the file properly. 
C 
      IF(L.EQ.60) THEN 
        RC = 0 
        GOTO 99 
      ELSE 
        WRITE(IOUOUT,170,ERR=98) TRIM(FILENAME), L 
        WRITE(IOUERR,170,ERR=98) TRIM(FILENAME), L 
 170    FORMAT('*** ERROR:  invalid number of data lines in the'/ 
     &         '    external PM carbon ZML data file ',A,'. The number'/ 
     &         '    of lines read in was ',I3,', out of 60. Please'/ 
     &         '    correct the file, especially checking for the'/ 
     &         '    proper termination of the last line.') 
        RC = 2 
        GOTO 98 
      END IF 
C 
C     Close the datafile and write the appropriate warning message. 
C 
  80  WRITE(IOUOUT,180,ERR=98) FILENAME 
      WRITE(IOUERR,180,ERR=98) FILENAME 
 180  FORMAT('*** Error opening external PM carbon ZML data file ',A) 
      RC = 2 
      GOTO 98 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in'/ 
     &            '    the external PM carbon ZML data file: ',A) 
          RC = 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the '/ 
     &            '    external PM carbon ZML data file ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 98 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) FILENAME, TRIM(M6REC) 
  190 FORMAT(/'*** ERROR: the following record is not allowed in the'/ 
     &        '    header of the external PM carbon ZML data file ',A/ 
     &        1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Missing or invalid fuel data. 
C 
   91 WRITE(IOUERR,191) FILENAME, TRIM(M6REC) 
  191 FORMAT(/'*** ERROR: the following record contains invalid Fuel'/ 
     &        '    data in the external PM carbon ZML data file ',A/ 
     &       1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid catalyst data. 
C 
   92 WRITE(IOUERR,192) FILENAME, TRIM(M6REC) 
  192 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    Catalyst data in the external PM carbon ZML'/ 
     &        '    data file ',A 
     &       /1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid air pump data. 
C 
   93 WRITE(IOUERR,193) FILENAME, TRIM(M6REC) 
  193 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    Air Pump data in the external PM carbon ZML'/ 
     &        '    data file ',A 
     &       /1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid air pump data. 
C 
   94 WRITE(IOUERR,194) FILENAME, TRIM(M6REC) 
  194 FORMAT(/'*** ERROR: the following record contains gasoline'/ 
     &        '    carbon data for a diesel vehicle class in the '/ 
     &        '    external PM carbon ZML data file ',A 
     &       /1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   95 WRITE(IOUERR,195) FILENAME, TRIM(M6REC) 
  195 FORMAT(/'*** ERROR: the following record contains invalid'/ 
     &        '    indexes or data in the external PM carbon ZML'/ 
     &        '    data file ',A 
     &       /1X,A) 
      RC = 2 
      GOTO 98 
C 
C     Invalid indexes or data. 
C 
   96 WRITE(IOUERR,196) FILENAME, TRIM(M6REC) 
  196 FORMAT(/'*** ERROR: the following record contains invalid PM'/ 
     &        '    DR age indexes in the external PM carbon ZML'/ 
     &        '    data file ',A 
     &       /1X,A) 
      RC = 2 
C 
  98  INERR = INERR + 1 
C 
  99  CLOSE(IOUALT) 
      RETURN 
      END 
