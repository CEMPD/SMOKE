      SUBROUTINE GETPRS(INPREC,INERR,RC) 
C 
C  GETPRS gets and verifies pressure check program parameters and effectiveness 
C  rates. 
C 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK, QUITER and YRTEST. 
C 
C  Changes: (Last change first) 
C 
C  13 Oct 00 AIR Task 02: Bug Fix 279: Changed routine to free format read 
C                         with 2 or 4 digit years. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  22 Oct 1998 @DynTel-ZK 2-000 Explicit Typing  
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IOUCOM/ IOUOUT 
C    /YEARS4/ IY1960,IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR,INPREC,RC 
C    common blocks: 
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS 
C    /PRSCH2/ PRSFLG,PRSPGM,PRSFQT,CRPRS 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  LVT      I    flag indicating whether or not a vehicle type is 
C                covered by the Pressure Check program: value must be 
C                either 1 (no) or 2 (yes). 
C  DATAFIELD  C   One blank-delimited data item from the input record. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  X          R   Temporary real to integer variable 
C 
C 
C  Notes: 
C 
C  it was added in MOBILE5.3 so that GETPRS is called when PRSFLG=2 
C 
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C       Change PRSPGM range checking in the IF statement 
C 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. Vehicle type vectors were expanded. 
C 
C Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
      IMPLICIT NONE       
      INCLUDE 'PRSCH1.I' 
      INCLUDE 'PRSCH2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'YEARS4.I' 
C       
      CHARACTER*10   DATAFIELD 
      CHARACTER INPREC*(*) 
      INTEGER   IVTAM, LVT 
      INTEGER   INERR 
      INTEGER   RC 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
C 
      RC=0 
      PRSFLG = 2 
C 
C  Look for the first data field following the delimiter in column 20. 
C  Initialize by setting PTR2 to 21 (the position of the blank following 
C  the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
C  Acquire and verify the start year, first year covered, and last year covered. 
C  Years can be 2 or 4 digits long. 
C 
      PTR2 = 21 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPRSY=INT(X) 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPR1ST=INT(X) 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPRLST=INT(X) 
C 
C  Verify parameters: 
C 
      CALL YRTEST(LPRSY,75,IY1960,IYEND,INERR) 
      CALL YRTEST(LPR1ST,76,IY1960,IYEND,INERR) 
      CALL YRTEST(LPRLST,77,IY1960,IYEND,INERR) 
C 
C  read and verify the vehicle types covered. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.5) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVPRS(IVTAM),IVTAM=1,5) 
  110 FORMAT(8I1) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVPRS(IVTAM),IVTAM=6,13) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.1) GOTO 96 
      READ (DATAFIELD,110,ERR=96) LVPRS(14) 
C 
      DO 10 IVTAM=1,14 
      LVT=LVPRS(IVTAM) 
      IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR) 
   10 CONTINUE 
C 
C  Obtain and verify the inspection program type, inspection frequency, 
C  and compliance rate. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.2) GOTO 96 
      READ (DATAFIELD,110,ERR=96) PRSPGM,PRSFQT 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) CRPRS 
C 
      IF(PRSPGM.LT.1.OR.PRSPGM.GT.3) CALL QUITER(0.,PRSPGM,103,INERR) 
      IF(PRSFQT.LT.1.OR.PRSFQT.GT.2) CALL QUITER(0.,PRSFQT,104,INERR) 
      IF(CRPRS.LT.0..OR.CRPRS.GT.100.) CALL QUITER(CRPRS,1,105,INERR) 
C 
      RETURN 
C 
C  "Fatal" error in trying to read ATP program description 
C 
   96 RC = 1 
C 
   99 CALL QUITER(0.,0,107,INERR) 
      RETURN 
      END 
