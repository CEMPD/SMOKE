      SUBROUTINE GETNAT(INPREC,INERR,RC) 
C 
C  GETNAT reads in and processes refueling loss related HC parameters 
C  that are to apply nationwide - that is, to all scenarios of the run. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER, YRTEST, NXTREC, NXTTOK 
C 
C  Changes: (Last change first) 
C 
C  26 Oct 00 AIR Task 02 Bug Fix 279: Allow for 1 or 2 line input. 
C  12 Oct 00 AIR Task 02 Bug Fix 279: Input changed to free format 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  20 Oct 1998 @DynTel 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS3/ RLFLAG 
C    /IOUCOM/ IOUGEN 
C    /YEARS4/ IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /RLCOM2/ IS2SY,NPHASE,S2EFF 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  -------------------------------------------------------------- 
C  LOWVYR    I    lowest allowed S2 or ob vrs installation start year 
C  PCTEL     I    percent efficiency of S2 vrs in reducing LDGV/T refueling loss 
C  PCTEH     I    percent efficiency of S2 vrs in reducing HDGV   refueling loss 
C  PTR1      I    NXTTOK pointer 
C  PRT2      I    NXTTOK pointer 
C  DATAFIELD C    temporary character to number holder 
C  X         R    temporary real to integer holder 
C 
C  Notes: 
C 
C    5-May-1994 @ CSC-tmm request 2-411 
C       Add code to inform the user that onboard vrs parameters are ignored. 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. Additional records are required for new vehicle 
C           types. 
C 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C  Jun-24-98 @ DynTel-bag 2-000 Changed QUITER argument RVALUE (9 --> 10) 
C             to reflect changes made to the NAMFLG array. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'RLCOM2.I' 
      INCLUDE 'YEARS4.I' 
C 
      INTEGER, INTENT(INOUT) :: INERR 
      INTEGER, INTENT(OUT) :: RC 
      INTEGER     I 
      INTEGER, PARAMETER :: LOWVYR=1989 
      CHARACTER INPREC*(*) 
      CHARACTER*10   DATAFIELD 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           PCTEL 
      REAL           PCTEH 
      REAL           X 
C 
      RC = 0 
      RLFLAG = 2 
C 
C  The Stage II data can be on the same line as the command, or on the next line. 
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
C  Vapor Recovery Systems (vrs) input parameters section: 
C 
C  Stage II vrs input: start year, no. of phase-in years & percent efficiency 
C  for LDG & HDG are entered.  Years can be 2 or 4 digits long. 
C   
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96,END=95) X 
      IS2SY =INT(X) 
      CALL YRTEST(IS2SY,87,LOWVYR,IYEND,INERR) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96,END=95) X 
      NPHASE=INT(X) 
      IF(NPHASE.EQ.0) CALL QUITER(0.,NPHASE,86,INERR) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96,END=95) PCTEL 
      IF(PCTEL.LT.0.0.OR.PCTEL.GT.100.0) CALL QUITER(PCTEL,1,99,INERR) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=96,END=95) PCTEH 
      IF(PCTEH.LT.0.0.OR.PCTEH.GT.100.0) CALL QUITER(PCTEH,2,99,INERR) 
C 
      S2EFF(1)=PCTEL/100. 
      S2EFF(2)=S2EFF(1) 
      S2EFF(3)=S2EFF(2) 
      S2EFF(4)=S2EFF(3) 
      S2EFF(5)=S2EFF(4) 
      S2EFF(6)=PCTEH/100. 
      DO 5 I = 7,14 
        S2EFF(I)=S2EFF(6) 
    5 CONTINUE 
C 
C  The OBD code is obsolete in Mobile6 and has been removed. 
C 
   30 RETURN 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 1 
      GOTO 99 
C 
   96 RC = -1 
C 
   99 CALL QUITER(0.,0,130,INERR) 
      RETURN 
      END 
