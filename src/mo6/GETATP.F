      SUBROUTINE GETATP(INPREC,INERR,RC) 
C 
C  GETATP gets and verifies anti-tampering program parameters and effectiveness 
C  rates. 
C 
C  Called by DSPTCH. 
C 
C  Calls ATPAER, CALAER, NXTTOK, ,NXTREC, QUITER, and YRTEST. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; 
C            added dummy variables for model year, MYGTP, and vehicle yearly 
C            index, IDXGTP, for call to ATPAER 
C  26 Oct 00 AIR Task 02 Bug Fix 279: Allow for 1 or 2 line input. 
C  13 Oct 00 AIR Task 02: Bug Fix 279: Changed routine to free format read 
C                         with 2 or 4 digit years. 
C  07 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IOUCOM/ IOUGEN, IOUREP 
C    /YEARS4/ IY1960,IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG 
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP 
C    /EGSCAL/ AER 
C    /EVPDAT/ EVPAER 
C    /FLAGS3/ ATPFLG 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  COLSUM   R    column sum values from AER: each column's sum must be 
C                within the closed interval from 0.0 to 1.0. 
C  EVP      R    evaporative A-T effectiveness rate: each rate must be 
C                within the closed inetrval from 0.0 to 1.0. 
C  LVT      I    flag indicating whether or not a vehicle type is 
C                covered by the A-T program: value must be either 1 (no) 
C                or 2 (yes). 
C  RATE     R    A-T effectiveness rate against exhaust emissions: each 
C                rate must be within the closed interval from 0.0 to 1.0 
C  DATAFIELD  C   One blank-delimited data item from the input record. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  X          R   Temporary real to integer variable 
C 
C  Notes: 
C 
C  GETATP was modified in MOBILE5. 
C  ATP vs I/M error checking is not done it GETATP.  This is now done in 
C  three places, ONESEC, QUITER and CHCKIM. (In increasing order of 
C  responsibility and importance). 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. 
C  OCt-27-97 DynTel@gjr 2-645 Removed code for interactive input. 
C 
      IMPLICIT NONE  
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'EVPDAT.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'YEARS4.I' 
C 
      CHARACTER INPREC*150 
      INTEGER   INERR 
      INTEGER   RC 
      REAL EVP, COLSUM, RATE 
      INTEGER IVTAM, LVT, I, IAY, IHG, IPG, ICOL, IROW 
      CHARACTER*10   DATAFIELD 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
      INTEGER        MYGTP,IDXGTP 
C 
C  Identify program parameters: start year, first year covered, last 
C  year covered, vehicle types covered, inspection program type, 
C  inspection frequency, compliance rate, and disablement inspections 
C  conducted. 
C  M6 transition: users enter 5 values for LDGV/LDGT and 8 for HDGV. 
C 
C  Initializations: 
      RC = 0 
      ATPFLG = 2 
C 
C  The ATP data can be on the same line as the command, or on the next line. 
C  Therefore, first look for the first data field following the delimiter 
C  in column 20.  Initialize by setting PTR2 to 21 (the position of the  
C  blank following the delimiter). PTR1 is set by the NXTTOK subroutine. 
C  If no data is found (TRIM(LEN(DATAFIELD))=0), then skip to next line 
C  and read data starting at column 1 (PTR2=0).  
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
C 
      IF(LEN(TRIM(DATAFIELD)).EQ.0) THEN 
        INPREC = ' ' 
        CALL NXTREC(IOUGEN, IOUREP, INPREC, INERR, RC) 
        SELECT CASE (RC) 
          CASE (:-1) 
            GOTO 95 
          CASE (1:) 
            GOTO 96 
          CASE DEFAULT 
        END SELECT 
        PTR2 = 0 
      ELSE 
        PTR2 = 21 
      ENDIF 
C 
C  Acquire and verify the start year, first year covered, and last year covered. 
C  Years can be 2 or 4 digits long. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LAPSY=INT(X) 
      MYGTP=LAPSY 
      IDXGTP=1 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LAP1ST=INT(X) 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) X 
      LAPLST=INT(X) 
C 
C  Verify parameters: 
C 
      CALL YRTEST(LAPSY,75,IY1960,IYEND,INERR) 
      CALL YRTEST(LAP1ST,76,IY1960,IYEND,INERR) 
      CALL YRTEST(LAPLST,77,IY1960,IYEND,INERR) 
C 
C  Obtain and verify the vehicle types covered. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.5) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVTFLG(IVTAM),IVTAM=1,5) 
  110 FORMAT(8I1) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (LVTFLG(IVTAM),IVTAM=6,13) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.1) GOTO 96 
      READ (DATAFIELD,110,ERR=96) LVTFLG(14) 
C 
      DO 10 IVTAM=1,14 
      LVT=LVTFLG(IVTAM) 
      IF(LVT.LT.1.OR.LVT.GT.2) CALL QUITER(0.,LVT,78,INERR) 
   10 CONTINUE 
C 
C  Obtain and verify the inspection program type, inspection frequency, 
C  and compliance rate. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.2) GOTO 96 
      READ (DATAFIELD,110,ERR=96) ATPPGM,ATPFQT 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96) CRATP 
C 
      IF(ATPPGM.LT.1.OR.ATPPGM.GT.2) CALL QUITER(0.,ATPPGM,103,INERR) 
      IF(ATPFQT.LT.1.OR.ATPFQT.GT.2) CALL QUITER(0.,ATPFQT,104,INERR) 
      IF(CRATP.LT.0..OR.CRATP.GT.100.) CALL QUITER(CRATP,1,105,INERR) 
C 
C  Obtain and verify the disablement types 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 96 
      READ (DATAFIELD,110,ERR=96) (DISTYP(IVTAM),IVTAM=1,8) 
C 
      DO 20 I=1,8 
      IF(DISTYP(I).LT.1.OR.DISTYP(I).GT.2) 
     *   CALL QUITER(0.,DISTYP(I),106,INERR) 
   20 CONTINUE 
      IF(DISTYP(6).EQ.2.AND.DISTYP(8).NE.2) 
     *   CALL QUITER(0.,0,121,INERR) 
C 
C  Assign base ATP effectiveness rates. 
C 
      CALL ATPAER(MYGTP,IDXGTP,1,ATPPGM,ATPFQT,CRATP,DISTYP) 
C 
C  Calculate ATP effectiveness rates. 
C 
      CALL CALAER 
C 
C  EVPAER and AER checks. 
C 
      DO 40 IAY=1,2 
      DO 30 IHG=1,3 
      EVP=EVPAER(IHG,IAY) 
      IF(EVP.LT.0.0.OR.EVP.GT.1.0) CALL QUITER(EVP,0,79,INERR) 
   30 CONTINUE 
   40 CONTINUE 
C 
C  Note: ICOL / IROW here refers to the input file layout of the arrays. 
C  In subroutine ATPAER the labels are switched, to indicate the actual 
C  orientation of the matrices in the effectiveness computation process. 
C 
      DO 80 IAY=1,2 
      DO 70 IPG=1,2 
      DO 60 ICOL=1,11 
      COLSUM=0.0 
      DO 50 IROW=1,11 
      RATE=AER(ICOL,IROW,IPG,IAY) 
      IF(RATE.LT.0.0.OR.RATE.GT.1.0) CALL QUITER(RATE,0,80,INERR) 
      COLSUM=COLSUM+RATE 
   50 CONTINUE 
      IF(COLSUM.LT.0.0.OR.COLSUM.GT.1.0) CALL QUITER(COLSUM,0,81,INERR) 
   60 CONTINUE 
   70 CONTINUE 
   80 CONTINUE 
C 
      RETURN 
C 
C  "Fatal" error in trying to read ATP program description 
C 
   95 RC = -1 
      GOTO 99 
C 
   96 RC = 1 
C 
   99 CALL QUITER(0.,0,107,INERR) 
      RETURN 
      END 
