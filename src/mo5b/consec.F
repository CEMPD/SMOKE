C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CONSEC(*,*)
C
C  CONSEC reads in and validates the Control Section flags.  It also sets 
C  the input mode for the run, resets the output units (if requested) and 
C  reads in and stores the run title.
C
C  Called by MOBILE.
C
C  Calls GETHOR, GETIOU, GETVER, MOD (FORTRAN library function), OUTHDR,
C  and QUITER.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUGEN,IOUREP,IOUASK
C
C  Output on return:
C
C    common blocks:
C    /FLAGS1/ PROMPT, PPTFLG
C    /LEVINF/ LEVNAM
C    /PROJEC/ PROJID
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  FEXIST   L    flag that indicates wheather the filename provided by the 
C                user exists.
C  FLGERR   I    number of flag values found to be out of range
C  INFORM   I    dummy variable: QUITER call purpose is to write out run
C                'inform'ation, not an error or a warning message
C  IPRERR   I    PROMPT flag error code: 1 = in range (1-6)
C                                        2 = out of range
C  VHMODE   I    decides prompt & read format for Control Section flags:
C                  1 = vertical format     2 = horizontal format
C
C  Notes:
C
C  CONSEC was modified in Version 01 to provide a second alternative 
C  RETURN for a GETIOU failure.  21-May-1993 @ ARC-bsg Subtask 222, 
C  preventing the project id from being written to the CSV file.
C
C  August-11-1994 @ CSC-pme request 2-455. CONSEC was modified to maximize
C  the flexibility of Mobile5b  by letting the user to provide Zero
C  Mile Levels (ZML) and Deterioration rates (DET) for the LEV
C  standards. A new level of the Prompt flag was added.
C
C  August-23-1994 @CSC-pme request 2-455. Major changes were done to 
C  subroutine CONSEC to ease the flow of data.
C
      INCLUDE 'BASEQ9.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVINF.I'
      INCLUDE 'PROJEC.I'
C
      INTEGER FLGERR, IOLEV, VHMODE
      CHARACTER*40 BLANK
      LOGICAL FEXIST
C
C     Initialize the LEV phase-in exist file flag, the prompt error flag,
C     the permanent prompt flag, and the blank constant; and set the I/O
C     LEV file unit number.
C
      FEXIST = .FALSE.
      IOLEV  = 18
      IPRERR = 1
      PPTFLG = 0
      BLANK = ' '
C
C     First pick up the I/O parameters (prompting or not, vertical or 
C     horizontal control flag input, and alternate output units) and 
C     identify & assign (if any) any allowable entered output units.
C
      CALL GETIOU(IOUSUM,*999)
C
C     Store prompt value into permanent prompt value flag. The value
C     is stored for future reference even if it is out of range. 
C
      PPTFLG = PROMPT
C
C     Check if entered value is in range.  Skip parsing if it is not.
C 
      IF(PROMPT .GE. 1 .AND. PROMPT .LE. 6) THEN
C
C       Lev Testing Option Reading in the LEVIMP
C
        IF (PROMPT.EQ. 5 .OR. PROMPT .EQ. 6) THEN
C
    5     WRITE (IOUASK,200)
  200     FORMAT (/'&', 'Please enter the LEVIMP input filename:')
          READ (IOUASK,100) LEVNAM
  100     FORMAT (A40)
C
C         Exit  if <return> is hit for filename prompts.
C
          IF (LEVNAM .EQ. BLANK) GOTO 99
C
C         Verify that the program exists...
C
          INQUIRE(FILE=LEVNAM,EXIST=FEXIST)
C
C         If the file does not exist ask the filename again
C
          IF (.NOT. FEXIST) GOTO 5
C
C         Open the existing file...
C
          OPEN (IOLEV,FILE=LEVNAM,STATUS='OLD')
C
C         Reading the file that contains the LEV phase-in schedule.
C
          DO 10 I = 1, 5
            DO 20 J = 1, 12
              READ (IOLEV,110,ERR=80,END=80) (LEVIMP(K,J,I),K=1,8)
  110         FORMAT (F5.3, 7(1X,F5.3))
   20       CONTINUE
   10     CONTINUE
C
C         Identify if it is necessary to read the second block of
C         the data file (PROMPT = 6)
C
          IF(PROMPT.EQ.6) THEN
C
C           Reading ZML, DET1, and DET2 from the given 
C           input data file.
C
            DO 30 IIP = 1, 3
              DO 40 IIV = 1, 5
                DO 50 ILEV = 1, 6
                READ (IOLEV,120,ERR=90,END=90) 
     *                                (LEVSTD(IX,ILEV,IIV,IIP),IX=1,3)
 120              FORMAT (F6.4,1X,F6.4,1X,F6.4)
  50            CONTINUE
  40          CONTINUE
  30        CONTINUE
          ENDIF
C
C         Reset flag to 1 a normal LEV run and close the data file.
C
          PROMPT = 1
          CLOSE(IOLEV)
        ENDIF
C
C       Parse prompt/mode flag.
C
C       Default mode to vertical.
C
  60    VHMODE=1
C
C       If prompt value greater than 2 set mode to horizontal.
C
        IF(PROMPT.GE.3)VHMODE=2
C
C       Re-define variable PROMPT to 1: do not prompt for input  or 
C                                    2: prompt for input.
C
        PROMPT=MOD(PROMPT-1,2)+1
C
C     Otherwise the PROMPT value is out of range... 
C
      ELSE
C
C       Set the PROMPT error flag to two.
C
        IPRERR=2
      ENDIF
C
C     Read title, even if PROMPT flag is out of range.
C
   70 IF(PROMPT.EQ.2.OR.IPRERR.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(///'&','Enter project id:')
C
      READ(IOUGEN,130,END=99) PROJID
  130 FORMAT(20A4)
C
C     Write run title to each i/o unit (up to 3, depending on 
C     GETIOU's IOUCHG).
C
      CALL OUTHDR(IOUERR)
C
CC     IF(IOUERR.NE.IOUREP) CALL OUTHDR(IOUERR)
CC     IF(PROMPT.EQ.2.AND.IOUASK.NE.IOUREP.AND.IOUASK.NE.IOUERR)
CC    *   CALL OUTHDR(IOUASK)
C
C         Display unit number values if changed.
C
      IF(IOUSUM.GT.0) CALL QUITER(0.,0,71,INFORM)
C
C     If PROMPT is out of range, tell user and stop run.
C
      FLGERR=0
      IF(IPRERR.EQ.2) CALL QUITER(0.,PROMPT,89,FLGERR)
C
C     Ready to read in flags in either vertical or horizontal mode.
C
      IF(VHMODE.EQ.1) CALL GETVER(FLGERR,*99)
      IF(VHMODE.EQ.2) CALL GETHOR(FLGERR,*99)
C
C     One or more flag values out of range => warn user & then stop run.
C
      IF(FLGERR.GT.0) CALL QUITER(0.,0,53,FLGERR)
      RETURN
C
C     EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
C
C     Close the datafile and write the appropriate warning message.
C
  80  CLOSE(IOLEV)
      WRITE(IOUASK,220) IOLEV, LEVNAM
      WRITE(IOUREP,220) IOLEV, LEVNAM
      WRITE(IOUREP,230)K,J,I,LEVIMP(K,J,I)
C
 220  FORMAT(' ','Warning ',/,
     *        'Error reading the LEV phase-in IO =',I2,' FILE =',A40)
 230  FORMAT((3(2X,I2),1X,F10.4)) 
C
      GOTO 999
C
C     Close the datafile and write the appropriate warning message.
C
  90  IF (PROMPT .EQ. 6) THEN
C
        CLOSE(IOLEV)
        WRITE(IOUREP,240) IOLEV, LEVNAM
        WRITE(IOUASK,240) IOLEV, LEVNAM
        WRITE(IOUREP,250)(LEVSTD(IX,ILEV,IIV,IIP),IX=1,3)
C
 240    FORMAT(' ','Warning ',/,
     *        'Error reading the LEV emission factors IO =',
     *                                     I2,' FILE =',A40)
      ENDIF
 250  FORMAT (F6.4,1X,F6.4,1X,F6.4)
C
 999  RETURN2
      END
