      SUBROUTINE GETATHAPS(INPREC,INERR,RC) 
C 
C  GETATHAPS obtains the user-supplied Air Toxic Exhaust emission factors. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,UCCOMP,GETPMZML,GETPMDR1,GETPMDR2 
C 
C  Changes: (Last change first) 
C 
C  31 Aug 01: AIR Task 18: New Subroutine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  FILENAME   C   Input filename for PM ZML or DR values.  
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'STRING.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER DATAFIELD*80 
      CHARACTER FILENAME*80 
      CHARACTER LABEL*40 
      CHARACTER M6REC*150 
      CHARACTER STDLAB*40 
      INTEGER   I 
      INTEGER   IATH 
      INTEGER   IB 
      INTEGER   IE 
      INTEGER   IOS 
      INTEGER   L 
      INTEGER   LABTYPE 
      INTEGER   LDGIMP 
      INTEGER   NCOM 
      INTEGER   PTR1 
      INTEGER   PTR2 
      INTEGER   SECTION 
      INTEGER   TLOC(10) 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C    Initially set ATHAPS to abort in routine ATCALX. 
C 
      ATHAPS = .FALSE. 
C 
C    Look for the first data field (Add HAPS filename) following the delimiter 
C    in column 20. Initialize by setting PTR2 to 21 (the position of the blank 
C    following the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 

ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)

C 
C     Open the external data file for ADDITIONAL HAPS. 
C 
      OPEN(IOUALT,FILE=FILENAME,STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC,ERR=999) 
      WRITE(IOUREP,120) FILENAME 
  120 FORMAT(/'* Reading the Additional HAPS Rates '/ 
     &        '* from the external data file ',A) 
C 
C     Find the token number of the record label that must appear 
C     on a record somewhere in the header of the external data file. 
C 
      CALL CHKLAB('ATHAPS',STDLAB,LDGIMP,SECTION) 
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
        GOTO 999 
      ENDIF 
C 
C     Check for a valid labeled input record. 
C 
      CALL UCCOMP(M6REC(1:19),LABEL) 
      CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
      IF (LABTYPE.NE.LDGIMP) THEN 
        GOTO 999 
      ENDIF 
C 
C    Initialize the number of HAPS read in 
C 
      NATH=0 
C 
C    Loop until end of file. 
C 
      DO IATH = 1,50 
C 
C    Look for next valid record 
C 
      CALL NXTREC(IOUALT,IOUREP,M6REC,INERR,RC) 
      IF (RC.EQ.-1) EXIT  ! EOF 
C 
C    Find the locations of the comma characters in the input line. 
C 
      L = LEN_TRIM(M6REC) 
C 
      NCOM = 1 
      TLOC(NCOM)=1 
C 
      DO I=1,L 
      IF(M6REC(I:I).EQ.COMMA) THEN 
      NCOM=NCOM+1 
      TLOC(NCOM)=I 
      END IF 
      END DO 
C 
C  If not enough commas are found, file is not constructed correctly. 
C 
      IF(NCOM.LE.7) THEN 
      WRITE(IOUERR,77) NCOM,IATH,TRIM(M6REC) 
   77 FORMAT(' Error, the data in the Additional HAPS file must be'/ 
     *       ' comma-delimited and have at least 8 fields per line.'/ 
     *       ' NCOM=',I2,' IATH=',I2,' M6REC=',A) 
      GOTO 999 
      END IF 
C 
C  Force a comma as end of line character. Need to do this since EXCEL 
C  does not always terminate the last line properly. 
C 
      NCOM=NCOM+1 
      TLOC(NCOM)=L+1 
C 
C  PARSE THE DATA 
C 
C  Parse and check the pollutant number. It must be between 50 and 99. 
C 
      IB=TLOC(1) 
      IE=TLOC(2)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHPOL(IATH) 
      ELSE 
      WRITE(IOUERR,101) 
  101 FORMAT(' Error: The Additional HAPS Pollutant Number Field is', 
     *       ' null') 
      GOTO 999 
      END IF 
C 
      IF(ATHPOL(IATH).LT.50 .OR. 
     *   ATHPOL(IATH).GT.99) THEN 
      WRITE(IOUERR,11) ATHPOL(IATH) 
   11 FORMAT(' Error, the Additional HAPS Pollutant Number must be'/ 
     *       ' between 50 and 99, inclusive. ATHPOL(IATH)=',I3) 
      GOTO 999 
      END IF 
C 
C Parse the pollutant name. 
C 
      IB=TLOC(2)+1 
      IE=TLOC(3)-1 
      IF(IE-IB.GE.0) THEN 
      ATHPNAME(IATH)=TRIM(M6REC(IB:IE)) 
      ELSE 
      WRITE(IOUERR,102) 
  102 FORMAT(' Error: The Additional HAPS Pollutant Name Field is null') 
      GOTO 999 
      END IF 
C 
C Parse and verify the EKIND. Must be EXH or EVAP 
C 
      IB=TLOC(3)+1 
      IE=TLOC(4)-1 
      IF(IE-IB.GE.0) THEN 
      CALL UCCOMP(TRIM(M6REC(IB:IE)),ATHEKIND(IATH)) 
      ELSE 
      WRITE(IOUERR,103) 
  103 FORMAT(' Error: The Additional HAPS EKIND Field is null') 
      GOTO 999 
      END IF 
C 
      IF(ATHEKIND(IATH).NE.'EXH' .AND. 
     *   ATHEKIND(IATH).NE.'EVAP') THEN 
      WRITE(IOUERR,12) TRIM(ATHEKIND(IATH)) 
   12 FORMAT(' Error, the Additional HAPS EKIND must be either'/ 
     *       ' EXH or EVAP. ATHEKIND=',A) 
      GOTO 999 
      END IF 
C 
C Parse and check the vehicle number. It must be between 1 and 28. 
C 
      IB=TLOC(4)+1 
      IE=TLOC(5)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHVTYPE(IATH) 
      ELSE 
      WRITE(IOUERR,104) 
  104 FORMAT(' Error: The Additional HAPS VType Field is null') 
      GOTO 999 
      END IF 
C 
      IF(ATHVTYPE(IATH).LT.1 .OR. 
     *   ATHVTYPE(IATH).GT.28) THEN 
      WRITE(IOUERR,13) ATHVTYPE(IATH) 
   13 FORMAT(' Error, the Additional HAPS vehicle number must be'/ 
     *       ' between 1 and 28, inclusive. ATHVTYE=',A) 
      GOTO 999 
      END IF 
C 
C Parse and check the beginning model year. If below 1950- or 2051+, reset to the limits. 
C 
      IB=TLOC(5)+1 
      IE=TLOC(6)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHBEGMY(IATH) 
      ELSE 
      WRITE(IOUERR,105) 
  105 FORMAT(' Error: The Additional HAPS BegMY Field is null') 
      GOTO 999 
      END IF 
C 
      IF(ATHBEGMY(IATH).LT.1951) ATHBEGMY(IATH)=1951 
      IF(ATHBEGMY(IATH).GT.2050) ATHBEGMY(IATH)=2050 
C 
C Read and check the ending model year. If below 1950- or 2051+, reset to the limits. 
C 
C 
      IB=TLOC(6)+1 
      IE=TLOC(7)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHENDMY(IATH) 
      ELSE 
      WRITE(IOUERR,106) 
  106 FORMAT(' Error: The Additional HAPS EndMY Field is null') 
      GOTO 999 
      END IF 
C 
      IF(ATHENDMY(IATH).LT.1951) ATHENDMY(IATH)=1951 
      IF(ATHENDMY(IATH).GT.2050) ATHENDMY(IATH)=2050 
C 
C  Check the model year order 
C 
      IF(ATHENDMY(IATH).LT.ATHBEGMY(IATH)) THEN 
      WRITE(IOUERR,110) ATHBEGMY(IATH),ATHENDMY(IATH) 
  110 FORMAT(' Error: The model years are out of order'/ 
     *       ' BegMy=',I4,' EndMy=',I4) 
      GOTO 999 
      END IF 
C 
C Parse and verify the RatioBEF. Must be RATIOVOC, RATIOTOG, RATIOPM or BEF. 
C 
      IB=TLOC(7)+1 
      IE=TLOC(8)-1 
      IF(IE-IB.GE.0) THEN 
      CALL UCCOMP(TRIM(M6REC(IB:IE)),ATHRATIO(IATH)) 
      ELSE 
      WRITE(IOUERR,107) 
  107 FORMAT(' Error: The Additional HAPS Ratio/BEF Field is null') 
      GOTO 999 
      END IF 
C 
      IF(ATHRATIO(IATH).NE.'RATIOVOC' .AND. 
     *   ATHRATIO(IATH).NE.'RATIOTOG' .AND. 
     *   ATHRATIO(IATH).NE.'RATIOPM'  .AND. 
     *   ATHRATIO(IATH).NE.'BEF') THEN 
      WRITE(IOUERR,14) TRIM(ATHRATIO(IATH)) 
   14 FORMAT(' Error, the Additional HAPS Ratio/BEF must be either'/ 
     *       ' RATIOVOC, RATIOTOT, RATIOPM or BEF. ATHRATIO=',A) 
      GOTO 999 
      END IF 
C 
C  If the RATIOBEF is RATIOPM, then verify that the PARTICULATES command was also entered. 
C 
      IF(ATHRATIO(IATH).EQ.'RATIOPM') THEN 
      IF(.NOT.PARTFLG) THEN 
      WRITE(IOUERR,15)  
   15 FORMAT(' Error, the Additional HAPS entry specifies RATIOPM but'/ 
     *       ' no PARTICULATES command was entered.') 
      GOTO 999 
      END IF 
      END IF 
C 
C  If the RATIOBEF is RATIOVOC or RATIOTOG, then verify that the AIR TOXICS command was also entered. 
C 
      IF(ATHRATIO(IATH).EQ.'RATIOVOC' .OR. 
     *   ATHRATIO(IATH).EQ.'RATIOTOG') THEN 
      IF(.NOT.TOXFLG) THEN 
      WRITE(IOUERR,16) ATHRATIO(IATH),TOXFLG 
   16 FORMAT(' Error, the Additional HAPS entry specifies RATIOVOC or'/ 
     *       ' RATIOTOG but no AIR TOXICS command was entered.'/ 
     *       ' ATHRATIO=',A,' TOXFLG=',L7) 
      GOTO 999 
      END IF 
      END IF 
C 
C  Catch where EVAP is specified for a PM run 
C 
      IF(ATHRATIO(IATH).EQ.'RATIOPM' .AND. 
     *   ATHEKIND(IATH).EQ.'EVAP') THEN 
      WRITE(IOUERR,111) 
  111 FORMAT(' Error: the Additional HAPS entry RATIOPM cannot be'/ 
     *       'used with the EVAP parameter') 
      GOTO 999 
      END IF 
C 
C  Now obtain the Ratio or the ZML. The ZML value is in ug/mi. 
C 
      IB=TLOC(8)+1 
      IE=TLOC(9)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHRATZML(IATH) 
      ELSE 
      WRITE(IOUERR,108) 
  108 FORMAT(' Error: The Additional HAPS Ratio/ZML Field is null') 
      GOTO 999 
      END IF 
C 
C  If the RATIOBEF is BER, then obtain the deterioration rate, in ug/10K Miles. 
C 
      IF(ATHRATIO(IATH).EQ.'BER') THEN 
      IB=TLOC(9)+1 
      IE=TLOC(10)-1 
      IF(IE-IB.GE.0) THEN 
      READ(M6REC(IB:IE),*,ERR=999,END=999,IOSTAT=IOS) ATHDET(IATH) 
      ELSE 
      WRITE(IOUERR,109) 
  109 FORMAT(' Error: The Additional HAPS DET Field is null') 
      GOTO 999 
      END IF 
      ELSE 
      ATHDET(IATH)=0. 
      END IF 
C 
C  Update the number of HAPS successfully processed 
C 
      NATH=NATH+1 
C 
      END DO 
C 
      IF(NATH.GT.0) THEN 
      ATHAPS = .TRUE. 
      WRITE(IOUREP,125) NATH,FILENAME 
  125 FORMAT(/'* Successfully read in',I3,' Additional HAPS lines'/ 
     &        '* from the external data file ',A) 
      RC = 0 
      GOTO 99 
      END IF 
C 
 999  INERR=INERR+1 
      RC=2 
C 
   99 CLOSE (IOUALT) 
C 
      RETURN 
      END 
        
