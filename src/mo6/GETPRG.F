      SUBROUTINE GETPRG(INPREC,INERR,RC) 
C 
C  GETPRG gets and verifies purge check program parameters and effectiveness 
C  rates. 
C  This checks if a Purge Program is covered by an IM240. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK, QUITER, YRTEST 
C 
C  Changes: (Last change first) 
C 
C  13 Oct 00 AIR Task 02: Bug Fix 279: Changed routine to free format read 
C                         with 2 or 4 digit years. 
C  21 OCT 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS2/ IMFLAG 
C    /IMPAR1/ MODYR 
C    /IMPAR2/ ILDT, ITEST 
C    /YEARS4/ IY1960,IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR,INPREC,RC 
C    common blocks: 
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG 
C    /PRGCH2/ PRGFLG,PRGPGM,PRGFQT,CRPRG 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  LVT      I    flag indicating whether or not a vehicle type is 
C                covered by the Pressure Check program value: must be 
C                either (1) no or (2) yes. 
C  MY       I    see top of program. 
C  DATAFIELD  C   One blank-delimited data item from the input record. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  X          R   Temporary real to integer variable 
C 
C 
C  Notes: 
C 
C  GETPRG was added in MOBILE5.3 so that it is called when PRGFLG=2 
C  If the purge program is not covered by an IM240 the user gets 
C  a warning about dynomometer. 
C 
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C       Change PRGPGM range checking in the IF statement 
C 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. Vehicle type vectors were expanded. 
C 
C Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
      IMPLICIT NONE  
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'PRGCH1.I' 
      INCLUDE 'PRGCH2.I' 
      INCLUDE 'YEARS4.I' 
C       
      CHARACTER INPREC*(*) 
      INTEGER   INERR 
      INTEGER   RC 
C       
      INTEGER   IVTAM, LVT, ICALL, MY  
      CHARACTER*10   DATAFIELD 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
C 
      RC = 0 
      PRGFLG = 2 
C 
C  Look for the first data field following the delimiter in column 20. 
C  Initialize by setting PTR2 to 21 (the position of the blank following 
C  the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
C  Read and verify the start year, first year covered, and last year covered. 
C  Years can be 2 or 4 digits long. 
C 
      PTR2 = 21 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPGSY=INT(X) 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPG1ST=INT(X) 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LPGLST=INT(X) 
C 
C  Verify parameters: 
C 
      CALL YRTEST(LPGSY,75,IY1960,IYEND,INERR) 
      CALL YRTEST(LPG1ST,76,IY1960,IYEND,INERR) 
      CALL YRTEST(LPGLST,77,IY1960,IYEND,INERR) 
C 
C  Read and verify the vehicle types covered. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.5) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVPRG(IVTAM),IVTAM=1,5) 
  110 FORMAT(8I1) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVPRG(IVTAM),IVTAM=6,13) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.1) GOTO 96 
      READ (DATAFIELD,110,ERR=96) LVPRG(14) 
C 
      DO 10 IVTAM=1,14 
        LVT=LVPRG(IVTAM) 
        IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR) 
   10 CONTINUE 
C 
C  Read and verify the inspection program type, inspection frequency, 
C  and compliance rate. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.2) GOTO 96 
      READ (DATAFIELD,110,ERR=96) PRGPGM,PRGFQT 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) CRPRG 
      IF(PRGPGM.LT.1.OR.PRGPGM.GT.3) CALL QUITER(0.,PRGPGM,103,INERR) 
      IF(PRGFQT.LT.1.OR.PRGFQT.GT.2) CALL QUITER(0.,PRGFQT,104,INERR) 
      IF(CRPRG.LT.0..OR.CRPRG.GT.100.) CALL QUITER(CRPRG,1,105,INERR) 
C 
      ICALL=0 
      DO 60 IVTAM=1,14 
      IF(LVPRG(IVTAM).EQ.1) GOTO 60 
      IF(IMFLAG.EQ.1) THEN 
      ICALL=1 
      GOTO 60 
      ENDIF 
C 
      DO 50 MY=LPG1ST, LPGLST 
C 
      IF(IMFLAG.EQ.3) GOTO 30 
C 
C  Assume the one I/M program case 
C 
      IF(ILDT(IVTAM,1).EQ.1)  ICALL=1 
C 
      IF(ITEST(1).NE.4)  ICALL=1 
C 
      IF(MY.LT.MODYR(1,1))  ICALL=1 
C 
      IF(MY.GT.MODYR(2,1))  ICALL=1 
      GOTO 50 
C 
C  Assume the two I/M program case 
C 
   30 IF(ILDT(IVTAM,1).EQ.1) GOTO 40 
C 
      IF(ITEST(1).NE.4) GOTO 40 
C 
      IF(MY.LT.MODYR(1,1)) GOTO 40 
C 
C  Assume 1st I/M program is IM240 
C 
      IF(MY.LE.MODYR(2,1)) GOTO 50 
   40 IF(ILDT(IVTAM,2).EQ.1)  ICALL=1 
C 
      IF(ITEST(2).NE.4)  ICALL=1 
C 
      IF(MY.LT.MODYR(1,2))  ICALL=1 
C 
      IF(MY.GT.MODYR(2,2))  ICALL=1 
C 
   50 CONTINUE 
   60 CONTINUE 
C 
C  Warn users about dynomometer use with Purge Program. 
C 
      IF(ICALL.EQ.1) CALL QUITER(0.,0,114,INERR) 
      RETURN 
C 
C  "Fatal" error in trying to read ATP program description 
C 
   96 RC = 1 
C 
   99 CALL QUITER(0.,0,107,INERR) 
      RETURN 
      END 
