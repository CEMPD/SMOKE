C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETSC1(ICY,INERR,*)
C
C  GETSC1 reads in and validates the 1st required scenario parameter record.
C  Also read in the TECH1&2 I/M credits if the altitude is "new" and there is
C  an I/M program.
C
C  Called by PARSEC.
C
C  Calls FLOAT, IFIX, QUITER and YRTEST.  FLOAT & IFIX are FORTRAN library
C  functions.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT,SPDFLG,VMFLAG
C    /FLAGS2/ IMFLAG
C    /IMPAR6/ IFILET
C    /IMPAR7/ TNAME,IREAD,SVALTC
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IY1960,IYEND
C    /MAXIMA/ MAXVEH
C    /REGION/ FEET
C    /RSDINF/ RSDFLG
C    /TTCINF/ IHYB
C
C  Output on return:
C
C    parameter list: ICY,INERR
C    common blocks:
C    /CITCIN/ PFRATE,IPFCU
C    /EVAL/   MEVAL
C    /IMPAR4/ SMRTL1 
C    /IMPAR7/ IREAD,SVALTC
C    /IM12HC/ CR12HC
C    /IM12CO/ CR12CO
C    /LEVBLK/ LEVFLG,IMFLG,LEVYRS,APPLEV
C    /REGION/ IREJN,ALT
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /TEMPS/  AMBT
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BUFFER    C    To keep "Smart line" record 
C  CSMART(1) C   Tells if Hybrid program  
C  CSMART(2) C   Tells if RSD program
C  CSMART(3) C   Reserv 
C  DTNAME   C    Holds the description of TECH1&2
C  ICTECH   I    Used to check the TECH1&2 validity
C  IMSKIP   I    Counts how many lines to skip in the I/M file
C  ITALD    I    Used to compare the altitude
C  ITMPLV   I    Temporary local for LEVYRS
C  IV, ITECH, ISTR, IAGE, IBY
C           I    See top of program
C  LEVYRS   I    LEV Years as read in when Region=4
C  NSPDS    I    number of speeds entered by user for current scenario:
C                  SPDFLG = 1 => NSPDS = 1;  SPDFLG = 2 => NSPDS = 8
C  OPCHK1   R    operating modes check 1: 0.0 <= PCHC + PCCC <= 100.
C  OPCHK2   R    operating modes check 2: 0.0 <= PCHC + PCCC - PCCN <= 100.
C  TAPPLV   I    Temporary local for APPLEV flag
C
C
C  Notes:
C
C  AMBT is rounded to the nearest 0.1 degrees F.
C  GETSC1 was modified for MOBILE5 version 2 to allow the TECH1&2 credits to
C  be read after comparing with the altitude flag.  The logical place was to
C  make that check here.
C  OCTOBER-3-1994 @CSC-pme. Apply LEV stds flag (APPLEV) is now being read from
C  input file and its value is stored on LEVBLK common. Its default value is
C  '1' that means 'do not apply LEV stds to LDGT2'
C
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added checking "smart line" if Hybrid or RSD program.
C
C  30 May 96 @DynTel-MLA. 2-620  Added variable SVALTC to show whether
C    the previous request to read I/M credits was for standard credits
C    (from TECH12.D) or alternate credits (from a user-designated file).
C
      CHARACTER*120  BUFFER
      CHARACTER*1    CSMART(3)
      CHARACTER*80   DTNAME
      INTEGER        TAPPLV
C
      INCLUDE 'CITCIN.I'
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IMPAR7.I'
      INCLUDE 'IM12HC.I'
      INCLUDE 'IM12CO.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RSDINF.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'TTCINF.I'
      INCLUDE 'YEARS4.I'

cmh  LMOS Modification from Mobile5a
cmh  Add definition of ndi_di, TDNAME, ENVL
c
      CHARACTER*1 ndi_di
      CHARACTER*132 ENVL,TDNAME
c
cmh  LMOS Modification from Mobile5a
cmh  Add ndi_di in common block SCENE3
c
      COMMON /SCENE3/ ndi_di
c
cmh  LMOS Modification from Mobile5a
cmh  Added FLAGS3 common block so the the temperature flag could be modified
c
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
c
cmh  Add declarations and common block for temperature overide
      REAL             TMPRTR, TMIN, TMAX
      COMMON /TOVRIDE/ TMPRTR, TMIN, TMAX

cmh  Add declarations and common block for flag for writing messages
      LOGICAL WFLAG
      COMMON /WARNING/ WFLAG
C
C  initilized values for the LEVs
C
C  start model year for LEV vehicles
      LEVYRS=1994
      IF (IREAD .EQ. 0) SMRTL1=' '
C
C
C  flag to control writes showing LEV I/M credits
C
C  tempory local variable
      ITMPLV=LEVYRS
C
      IMFLG=2
C
C  defaults to ONE (1) .FALSE.
      LEVFLG=1
C
C  Default APPLEV to 1 (no LEV standards for LDGT2)
      APPLEV=1
C
C  GETSC1 is called once for each scenario of the run.
C
C  Get the first required scenario record.
C
      IF(SPDFLG.EQ.2) GOTO 10
      NSPDS=1
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter region, CY,SPD(1),AMBT,PCCN,PCHC,PCCC,Month: ')
c
cmh  LMOS Modification from Mobile5a
cmh  Add the non-diurnal/diurnal flag to the end of the scenario record
c
c     READ(IOUGEN,100,ERR=99,END=99)IREJN,ICY,SPD(1),AMBT,PCCN,PCHC,PCCC
c    *,MEVAL,(PFRATE(KEYEQ,2),KEYEQ=1,3)
c 100 FORMAT(I1,1X,I2,5(1X,F4.1),I3,3F5.0)
      READ(IOUGEN,100,ERR=99,END=99)IREJN,ICY,SPD(1),AMBT,PCCN,PCHC,PCCC
     *,MEVAL,(PFRATE(KEYEQ,2),KEYEQ=1,3),ndi_di
  100 FORMAT(I1,1X,I2,5(1X,F4.1),I3,3F5.0,A1)
c
cmh add temperature override
      AMBT = TMPRTR
c
C
C     check for irejn = 3 or 4 and the Processing of the LEVs
      IF(IREJN.EQ.3.OR.IREJN.EQ.4) THEN
	IF(IREJN.EQ.3) GOTO 15
	  IF(IREJN.EQ.4) THEN
	    IF(PROMPT.EQ.2) WRITE(IOUASK,220)
	    READ(IOUGEN,110,ERR=97,END=99) ITMPLV, IMFLG, TAPPLV
  110       FORMAT(I2,1X,I1,1X,I1)
	  ENDIF
      ENDIF
      GOTO 15
C
   10 NSPDS=MAXVEH
      IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ','Enter region, CY,SPD(8),AMBT,PCCN,PCHC,PCCC,Month: ')
c
cmh  LMOS Modification from Mobile5a
cmh  Add the non-diurnal/diurnal flag to the end of the scenario record
c
c     READ(IOUGEN,120,ERR=97,END=99) IREJN,ICY,SPD,AMBT,PCCN,PCHC,PCCC
c    *,MEVAL,(PFRATE(KEYEQ,2),KEYEQ=1,3)
c 120 FORMAT(I1,1X,I2,12(1X,F4.1),I3,3F5.0)
      READ(IOUGEN,120,ERR=97,END=99) IREJN,ICY,SPD,AMBT,PCCN,PCHC,PCCC
     *,MEVAL,(PFRATE(KEYEQ,2),KEYEQ=1,3),ndi_di
  120 FORMAT(I1,1X,I2,12(1X,F4.1),I3,3F5.0,A1)
c
C
C     check for irejn = 4 and the Processing of the LEVs
      IF(IREJN.EQ.4) THEN
	 IF(PROMPT.EQ.2) WRITE(IOUASK,220)
  220    FORMAT(' ','Enter Starting Model Year set the ',
     *   'I/M Flag and choose to apply LEV stds'/' ','for LDGT2',
     *   /,' ','I2 1X I1 1X I1 ','format ')
	 READ(IOUGEN,130,ERR=97,END=99) ITMPLV, IMFLG, TAPPLV
  130    FORMAT(I2,1X,I1, 1X, I1)
      ENDIF
C  Check for "end of run" record, a row of 60 zeroes.  This record acts as
C  visible break between 2 input streams, providing dummy values for all first
C  scenario record fields and leaving the file pointer at the first record of
C  the next input stream.  Multiple input streams in a single file or a set of
C  concatenated files could occur in using MOBILE4 as a subroutine, for example.
C  A normal end of run also occurs, of course, when the first scenario record
C  read encounters an end of file marker.
C
c
cmh  LMOS Modification from Mobile5a
cmh  Add assignment of temperature flag based on value of ndi_di
c
c  15 IF(IREJN.EQ.0) GOTO 99
   15 IF(ndi_di .EQ. 'N') THEN TEMFLG=2
        ELSE TEMFLG = 1
      IF(IREJN.EQ.0) GOTO 99
c
C
C  Otherwise, validate inputs.
C
C     if IREJN equals 3 or 4 Process the Low Emitting Vehicles
C     set the LEVFLG to 2 (TRUE) and set the IMFLG to the default of 2
C
C     if IREJN = 1: Low  Altitude, No LEVs just normal Federal
C
C     if IREJN = 2: High Altitude  No LEVs just normal Federal
C
C     if IREJN = 3: Low  Altitude  With LEVs
C                   Set LEVFLG => 2 (TRUE) intergrate LEVs into the Federal
C                   Standard Emission Factors starting in MY 1993
C                   with the Default Special I/M Credit.
C
C     if IREJN = 4: Low  Altitude With LEVs
C                   Set LEVFLG => 2 (TRUE) read in new LEV Record.
C                   The Lev Record is user defined. This LEV Record
C                   immediately follows the Scenairo Record in a
C                   I2,1X,I1 format, reading two Variables. 1) LEVYRS
C                   the Start Year of LEV integration into the Model
C                   and 2) the IMFLG this flag does or does not add in
C                   the Special I/M Credit.
C                   if IMFLG = 1: No Special I/M for the LEV *
C                   * if IREJN=4 and IMFLG=1 then the LEVs get the
C                      the 1992 Federal LDGV I/M Credit, but adjusted
C                      difference in the zero mile levels.
C                    if IMFLG = 2: Special I/M Credit for LEV
C
C
      IF (IREJN.EQ.3.OR.IREJN.EQ.4) THEN
C
C     process the LEVs starting with setting
C     the LEVFLG TO A TWO (2) .TRUE.
	 LEVFLG=2
C
	 IF(IREJN.EQ.4) THEN
C
C           imflg is initilized to a two (2)
C           if IMFLG = 1 => No Special I/M Credit ( just the 92 LDGV's )
C           if IMFLG = 2 => Special I/M test
C
C           call Quiter "FOR ERROR EXIT"
C
	    IF(IMFLG.NE.1.AND.IMFLG.NE.2)
     *      CALL QUITER(0.,IMFLG,10,INERR)
C
C           test the levyrs (in ITMPLV) from 1994 to 2051
	    CALL YRTEST(ITMPLV,32,IY1960,IYEND,INERR)
C
C           call Quiter "FOR ERROR EXIT"
C
	    IF(ITMPLV.LT.1994) CALL QUITER(0.,ITMPLV,11,INERR)
C
C           load the tested LEV Year into LEVYRS
	    LEVYRS=ITMPLV
C
C           If user chose to apply LEV stds to LDGT2
			IF(TAPPLV.EQ.2) APPLEV = 2 
C
	 ENDIF
C
C        Region is always Low Altitude if processing LEV's LEVFLG=2
C
	 IREJN=1
C
      ENDIF
C
C     Region error exit
      IF(IREJN.LT.1.OR.IREJN.GT.2) CALL QUITER(0.,IREJN,24,INERR)
C
C  Do not read any I/M credits unless IMFLG is greater than 1. Read new
C  credits whenever alternate credits are requested. Also read new credits 
C  whenever region (high or low altitude) changes, and when standard 
C  credits are requested for the current scenario but alternate credits 
C  were used for the last scenario. The only time we will skip reading new 
C  credits is when standard credits were used for the last scenario and 
C  standard credits (for the same region) are requested for the current 
C  scenario.
C
      IF(IMFLAG.EQ.1) GOTO 65
      IF (IFILET(1).EQ.1.AND.SVALTC.EQ.1.AND.IREAD.EQ.IREJN) GOTO 65
C
      IREAD = IREJN
      SVALTC = IFILET(1)
C
C     If the user supplied an alternate credit file name, it will be in
C     the variable TNAME. Otherwise, we can use TNAME to hold the standard
C     credit file name 'TECH12.D.'
C
cmh  LMOS Modification from Mobile5a
cmh
cmh  Direct user-supplied filename to be supplied from the EMS_LOC directory.
cmh  If the default file, TECH12.D, is supplied, get it from EMS_LOC
c
         IF(IFILET(1).EQ.2) THEN
            CALL GETENV("EMS_LOC",ENVL)
            TDNAME=ENVL(1:LNBLNK(ENVL))//'/'//TNAME
            OPEN(3,FILE=TDNAME,STATUS='OLD')
c            OPEN(3,FILE=TNAME,STATUS='OLD')
            GOTO 20
         ENDIF
         CALL GETENV("EMS_LOC",ENVL)
         TDNAME=ENVL(1:LNBLNK(ENVL))//'/TECH12.D'
         OPEN(3,FILE=TDNAME,STATUS='OLD')
c
c     IF(IFILET(1).EQ.1) TNAME='TECH12.D'
c     OPEN(3,FILE=TNAME,STATUS='OLD')
c
cmh end
C
   20 READ(3,141,END=95,ERR=96) BUFFER
      READ(UNIT=BUFFER,FMT=140,END=95,ERR=96) ITALD,ICTECH,
     *  (CSMART(ISMART), ISMART=1,3), DTNAME
C     
C Check MOBILE5b or MOBILE5a format
C
	 IF ((CSMART(1) .EQ. '1') .OR. (CSMART(1) .EQ. '2')) THEN 
C
C        Check the tech type and altitude and if they are not right then
C        jump back and get the next one after skipping 400 lines.
C
	    IF( (ITALD.NE.IREJN) .OR.(ICTECH.NE.1) .OR.
C                 Check if Hybrid program
     *         ((IHYB(1) .EQ. 2) .AND. (CSMART(1) .NE. '2')) .OR.
C                 Check if RSD program     
     *         ((RSDFLG .EQ.2) .AND. (CSMART(2) .NE. '2'))  ) THEN
	       GOTO 21
	    ELSE   
	       GOTO 30
	    ENDIF
	  ELSE     
C MOBILE5a format
	    READ(UNIT=BUFFER,FMT=142,END=95,ERR=96) ITALD,ICTECH,DTNAME
	    IF (ITALD.EQ.IREJN .AND. ICTECH.EQ.1 .AND. IHYB(1).NE.2
     *          .AND. RSDFLG.NE.2) GOTO 30
	  ENDIF  
C     
   21       DO 25 IMSKIP=1,400
	      READ(3,150)
   25       CONTINUE
  140    FORMAT(I1,1X,I1,T21,3A1,1X,A80)
  141    FORMAT(120A)
  142    FORMAT(I1,1X,I1,17X,A80) 
  150    FORMAT(1X)
	 GOTO 20
C
C        Reading the data in this manner accomplishes the same task as
C        Block Data Subprograms 06, and 07 once did.
C
C        Loop for the two (2) ITECHs to be stored in the record.
C
cmh add if statement
   30    IF( WFLAG ) THEN
             WRITE(IOUASK,240) DTNAME
         ENDIF
  240    FORMAT(7X,' Reading TECH1&2 I/M Credit Data File...'/ 7X,A80)

	 SMRTL1=DTNAME
	 DO 45 ITECH=1,2
C
C          Loop for the five (5) stringency values in the record.
C
	   DO 40 ISTR=1,5
C
C            Loop for the values in the record
C
	     DO 35 IAGE1=1,20
	       READ(3,160,END=95,ERR=96)
     *             (CR12HC(IBY,IAGE1,ISTR,ITECH),IBY=1,19)
   35        CONTINUE
C
   40      CONTINUE
   45    CONTINUE
C
C        Loop for the two (2) ITECHs to be stored in the record.
C
	 DO 60 ITECH=1,2
C
C          Loop for the five (5) stringency values in the record.
C
	   DO 55 ISTR=1,5
C
C            Loop for the values in the record
C
	     DO 50 IAGE1=1,20
	       READ(3,160,END=95,ERR=96)
     *             (CR12CO(IBY,IAGE1,ISTR,ITECH),IBY=1,19)
  160    FORMAT(19F4.3)
   50        CONTINUE
C
   55      CONTINUE
   60    CONTINUE
C
	 CLOSE(3)
C
C  Converts ICY to a four digit integer year
C
   65 CALL YRTEST(ICY,38,IY1960,IYEND,INERR)
C
C  Check month. Allow only Jan (MEVAL=0 or 1) and July (MEVAL=7)
C
      IF(MEVAL.LE.0) MEVAL=1
      IF(MEVAL.NE.1.AND.MEVAL.NE.7) CALL QUITER(0.,MEVAL,145,INERR)
      DO 70 IV=1,NSPDS
      IF(SPD(IV).LE.0.0) CALL QUITER(SPD(IV),0,39,INERR)
      IF(SPD(IV).GT.0.0.AND.SPD(IV).LT.2.5)
     *    CALL QUITER(SPD(IV),0,52,INERR)
      IF(SPD(IV).GT.0.0.AND.SPD(IV).LT.2.5) SPD(IV)=2.5
      IF(SPD(IV).GT.65.0) CALL QUITER(SPD(IV),0,96,INERR)
      IF(SPD(IV).GT.65.0) SPD(IV)=65.0
C
   70 CONTINUE
C
      IF(SPDFLG.NE.2) THEN
      DO 75 IV=2,MAXVEH
C
      SPD(IV)=SPD(1)
   75 CONTINUE
      ENDIF
C
      AMBT=FLOAT(IFIX(AMBT*10.+.5))/10.
      IF(AMBT.LT.0.0.OR.AMBT.GT.110.0) CALL QUITER(AMBT,0,40,INERR)
      IF(PCCN.LT.0.0.OR.PCCN.GT.100.0) CALL QUITER(PCCN,0,41,INERR)
      IF(PCHC.LT.0.0.OR.PCHC.GT.100.0) CALL QUITER(PCHC,0,42,INERR)
      IF(PCCC.LT.0.0.OR.PCCC.GT.100.0) CALL QUITER(PCCC,0,43,INERR)
      OPCHK1=PCHC+PCCC
      IF(OPCHK1.LT.0.0.OR.OPCHK1.GT.100.0)
     *    CALL QUITER(OPCHK1,0,63,INERR)
      OPCHK2=OPCHK1-PCCN
      IF(OPCHK2.LT.0.0.OR.OPCHK2.GT.100.0)
     *    CALL QUITER(OPCHK2,0,64,INERR)
C
C  Check user supplied Pass, Purge Fail and Pressure Fail rates
C  Since this is optional input, any value > 0 triggers user-rate usage.
C
      IPFCU=1
      PFSUM=0.
      DO 80 KEYEQ=1,3
      IF(PFRATE(KEYEQ,2).LT.0..OR.PFRATE(KEYEQ,2).GT.100.)
     *    CALL QUITER(PFRATE(KEYEQ,2),0,150,INERR)
      IF(PFRATE(KEYEQ,2).GT.0.) IPFCU=2
      PFSUM=PFSUM+PFRATE(KEYEQ,2)
      PFRATE(KEYEQ,2)=PFRATE(KEYEQ,2)/100.
   80 CONTINUE
      IF(IPFCU.EQ.2.AND.PFSUM.NE.100.) CALL QUITER(PFSUM,0,151,INERR)
C
C Assign altitude in feet corresponding to selected (low,high) region.
C
      ALT=FEET(IREJN)
C
      GOTO 98
C
   95 CALL QUITER(0.,0,54,INERR)
      CLOSE(3)
      GOTO 98
   96 CALL QUITER(0.,0,55,INERR)
      CLOSE(3)
      GOTO 98
   97 CALL QUITER(0.,0,132,INERR)
      GOTO 99
   98 RETURN
C
C  Everything below here is fatal
C
   99 RETURN 1
      END
