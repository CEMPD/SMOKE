C Version @(#)$Id$ $Source$ $Date$ 

C     SUBROUTINE CONSEC(*,*)
C
C  CONSEC reads in and validates the Control Section flags.  It also sets the
C  input mode for the run, resets the output units (if requested) and reads in
C  and stores the run title.
C
C  Called by MOBILE.
C
C  Calls GETHOR, GETIOU, GETVER, MOD (FORTRAN library function), OUTHDR and
C  QUITER.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUGEN,IOUREP,IOUASK
C
C  Output on return:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /PROJEC/ PROJID
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  FLGERR   I    number of flag values found to be out of range
C  INFORM   I    dummy variable: QUITER call purpose is to write out run
C                'inform'ation, not an error or a warning message
C  IPRERR   I    PROMPT flag error code: 1 = in range (1-4)  2 = out of range
C  VHMODE   I    decides prompt & read format for Control Section flags:
C                  1 = vertical format     2 = horizontal format
C
C  Notes:
C
C  CONSEC was modified in Version 01 to provide a second alternative RETURN,
C  for a GETIOU failure.
C
C
C     CHARACTER*4 PROJID
C     INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
C     COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C     COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C     COMMON /PROJEC/ PROJID(20)
C
C     INTEGER FLGERR,VHMODE
C
C  First pick up the I/O parameters (prompting or not, vertical or horizontal
C  control flag input, and alternate output units) and identify & assign (if
C  any) any allowable entered output units.
C
C     CALL GETIOU(IOUSUM,*999)
C
C  Now parse prompt / input mode flag.
C
C  Check if entered value is out of range.  Skip parsing if it is.
C
C     IPRERR=1
C     IF(PROMPT.LT.1.OR.PROMPT.GT.4) IPRERR=2
C     IF(IPRERR.EQ.2) GOTO 10
C
C  Parse prompt/mode flag.
C
C     VHMODE=1
C     IF(PROMPT.GE.3)VHMODE=2
C     PROMPT=MOD(PROMPT-1,2)+1
C
C  Read title, even if PROMPT flag is out of range.
C
C  10 IF(PROMPT.EQ.2.OR.IPRERR.EQ.2) WRITE(IOUASK,200)
C 200 FORMAT(///'&','Enter project id:')
C     READ(IOUGEN,100,END=99) PROJID
C 100 FORMAT(20A4)
C
C  Write run title to each i/o unit (up to 3, depending on GETIOU's IOUCHG).
C
C     CALL OUTHDR(IOUREP)
C     IF(IOUERR.NE.IOUREP) CALL OUTHDR(IOUERR)
C     IF(PROMPT.EQ.2.AND.IOUASK.NE.IOUREP.AND.IOUASK.NE.IOUERR)
C    *   CALL OUTHDR(IOUASK)
C
C  Display unit number values if changed.
C
C     IF(IOUSUM.GT.0) CALL QUITER(0.,0,71,INFORM)
C
C  If PROMPT is out of range, tell user and stop run.
C
C     FLGERR=0
C     IF(IPRERR.EQ.2) CALL QUITER(0.,PROMPT,89,FLGERR)
C
C  Ready to read in flags in either vertical or horizontal mode.
C
C     IF(VHMODE.EQ.1) CALL GETVER(FLGERR,*99)
C     IF(VHMODE.EQ.2) CALL GETHOR(FLGERR,*99)
C
C  1 or more flag values out of range => warn user & then stop run.
C
C     IF(FLGERR.GT.0) CALL QUITER(0.,0,53,FLGERR)
C     RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
C  99 RETURN1
C 999 RETURN2
C     END
      SUBROUTINE CONSEC(*,*)
C
C  CONSEC reads in and validates the Control Section flags.  It also sets the
C  input mode for the run, resets the output units (if requested) and reads in
C  and stores the run title.
C
C  Called by MOBILE.
C
C  Calls GETHOR, GETIOU, GETVER, MOD (FORTRAN library function), OUTHDR and
C  QUITER.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUGEN,IOUREP,IOUASK
C
C  Output on return:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /PROJEC/ PROJID
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  FLGERR   I    number of flag values found to be out of range
C  INFORM   I    dummy variable: QUITER call purpose is to write out run
C                'inform'ation, not an error or a warning message
C  IPRERR   I    PROMPT flag error code: 1 = in range (1-4)  2 = out of range
C  VHMODE   I    decides prompt & read format for Control Section flags:
C                  1 = vertical format     2 = horizontal format
C
C  Notes:
C
C  CONSEC was modified in Version 01 to provide a second alternative RETURN,
C  for a GETIOU failure.
C
C
      CHARACTER*4 PROJID
c
c  LMOS Modification 04/08/93
c  Add definition of ENVL
c
      CHARACTER*132 ENVL
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      REAL LEVSTD, LEVIMP, LDTAB
C
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /PROJEC/ PROJID(20)
C
      INTEGER FLGERR,VHMODE, IOLEV
c
c  LMOS Modification 04/08/93
c  Change length of LEVNAM to allow full directory pathname
c
c      CHARACTER*40 LEVNAM
      CHARACTER*132 LEVNAM
C
      REAL LEVGET
      DIMENSION LEVGET(8,12,5)
C
      EQUIVALENCE ( LEVIMP(1,1,1), LEVGET(1,1,1) )
C
CCCC  DATA LEVNAM /'INLEV.D'/
      IOLEV = 18
C
C
C  First pick up the I/O parameters (prompting or not, vertical or horizontal
C  control flag input, and alternate output units) and identify & assign (if
C  any) any allowable entered output units.
C
      CALL GETIOU(IOUSUM,*999)
C
C  Now parse prompt / input mode flag.
C
C  Check if entered value is out of range.  Skip parsing if it is.
C
      IPRERR=1
CC  Lev Testing Option Reaading in the LEVIMP
CC
      IF(PROMPT.LT.1.OR.PROMPT.GT.5) IPRERR=2
      IF(IPRERR.EQ.2) GOTO 10
C
      IF(PROMPT.NE.5) GOTO 225
C
      IF (PROMPT.EQ.5) THEN
CC
c  LMOS Modification 04/08/93
c  Comment out input of alternate LEV input file.  Assume file will be named
c  LEVIMP and be located in EMS_LOC
c
c        WRITE (IOUASK,60)
c   60   FORMAT (/'&', 'Please enter the LEVIMP input filename:')
c        READ (IOUASK,70) LEVNAM
c   70   FORMAT (A40)
       CALL GETENV("EMS_LOC",ENVL)
       LEVNAM=ENVL(1:LNBLNK(ENVL))//'/LEVIMP'
C
C       IOLEV is set to 18
        OPEN (IOLEV,FILE=LEVNAM,STATUS='OLD')
C
CC      LEVIMP(8,12,5) LEV Implemenitation
CC
        DO 23 I = 1, 5
        DO 25 J = 1, 12
        READ (IOLEV,170,ERR=224,END=225) (LEVGET(K,J,I),K=1,8)
  170   FORMAT (F5.3, 7(1X,F5.3))
C
   25   CONTINUE
   23   CONTINUE
C
        CLOSE (IOLEV)
C
C  Reset flag to 1 a normal LEV ren
        PROMPT = 1
C
        GOTO 225
C
      ENDIF
C
 224  CONTINUE
c
c  LMOS Modification 04/08/93
c
c  We'll leave this write statement in for IOUREP.  If such an error
c  is encountered, the data set will cause the SAS program to bomb, making
c  it certain the error is brought to operator attention (should it not be
c  caught by checking the run status report)
c
      WRITE(IOUREP,190) IOLEV, LEVNAM
      WRITE(IOUASK,190) IOLEV, LEVNAM
 190  FORMAT(' ','Warning ',/,
     *'Error reading the LEVIMP IO =',I2,' FILE =',A10 )
      GOTO 999
C
 225  CONTINUE
C
C  Parse prompt/mode flag.
C
      VHMODE=1
      IF(PROMPT.GE.3)VHMODE=2
      PROMPT=MOD(PROMPT-1,2)+1
C
C  Read title, even if PROMPT flag is out of range.
C
   10 IF(PROMPT.EQ.2.OR.IPRERR.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(///'&','Enter project id:')
      READ(IOUGEN,100,END=99) PROJID
  100 FORMAT(20A4)
C
C  Write run title to each i/o unit (up to 3, depending on GETIOU's IOUCHG).
C
c
c  LMOS Modification 04/08/93
c  Comment out writing of header records
c
c      CALL OUTHDR(IOUREP)
c      IF(IOUERR.NE.IOUREP) CALL OUTHDR(IOUERR)
c
C     IF(PROMPT.EQ.2.AND.IOUASK.NE.IOUREP.AND.IOUASK.NE.IOUERR)
C    *   CALL OUTHDR(IOUASK)
C
C  Display unit number values if changed.
C
      IF(IOUSUM.GT.0) CALL QUITER(0.,0,71,INFORM)
C
C  If PROMPT is out of range, tell user and stop run.
C
      FLGERR=0
      IF(IPRERR.EQ.2) CALL QUITER(0.,PROMPT,89,FLGERR)
C
C  Ready to read in flags in either vertical or horizontal mode.
C
      IF(VHMODE.EQ.1) CALL GETVER(FLGERR,*99)
      IF(VHMODE.EQ.2) CALL GETHOR(FLGERR,*99)
C
C  1 or more flag values out of range => warn user & then stop run.
C
      IF(FLGERR.GT.0) CALL QUITER(0.,0,53,FLGERR)
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
  999 RETURN2
      END
