C Version "@(#)$Id$ $Source$ $Date$ 

      PROGRAM DRIVER
C
C
C  DRIVER calls MOBILE5 as a subroutine.
C
C  Calls MOBILE.
C
C  Local array subscripts
C
C  M5OUT(40)  -  M5OUT ( I )
C
C
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  BATCH    C*1  flag for batch mode ('Y') versus not batch
C                mode ('N')
C  BATIN   C*40  batch file name
C  BLANK   C*40  blank filename detector
C  CEOPEN   L    tells if alternative diagnostic message file has
C                been opened, T=open, F=not open
C  DIAGMF   C*1  variable that holds alternative diagnostic
C                messge file name if CSV is selected
C  FEXIST   L    temporary detection for status of input file,
C                if input file or batch file exist, T=exist
C                F=does not exist.
C  NFILE    I    Number of total files in batch file
C  NUMDID   I    Number of the total files in the batch file that
C                were processed.
C  INTER    C    flag for interactive mode ('Y') versus file
C                mode ('N')
C  ISBLNK   L    flag for indicating status of output filenames,
C                F=does not contain blanks, T=contains blanks.
C  M5CALL   I    loop index for MOBILE5 call loop
C  M5IN    C*40  filename for input to MOBILE5
C  M5OUT    C*1  filename characters for output from MOBILE5.
C  WM5OUT  C*40  whole version of the M5OUT for use in the open
C                statement
C
C
C  Notes:
C
C  DRIVER was modified for MOBILE5 to add the CSV option and to
C  inform the user of errors due to invalid files, including
C  file not found errors.
C  14-June-1993 @ ARC-bsg Subtask 222
C
C  DRIVER was modified for MOBILE5 to comment out the CSV option.
C  16-March-1994 @ CSC-tmm Subtask 392
C  2-May-1994 @ ARC-rwk Request 2-425 Exit if <return> is hit for filename
C  prompts
C
C  12-Dec-1994 @ CSC-rwk Req. 2-000 Include files VDATA and VNAME were added
C  to make it easier to change the name of the version of the model.  This was
C  done for quality control.  COMMTS.I was added so that the comments would be
C  'attached' to the model.
C
C   6 Mar 96 @DynTel-MLA 2-000 'Changes for MAC driver' To prepare this code
C      for the MAC, change all occurrences of "C MAC" (beginning in column
C      one) to "" (the null string). Remember to reassign the I/O units,
C      then compile and execute as usual.
C 
C  29 Apr 96 @DynTel-MLA 2-000 'Changes for Mac driver' Removed the 'click
C      cancel' text from the completion message that is displayed on the
C      Mac. The new message reads 'Processing complete. Hit Return to exit
C      from the program.'
C
C  3-Sep-96 @DynTel-yc Request 0-448  'Code change for removing warning 
C           messages'  Include statement for COMMTS.I was moved to just 
C           before the END statement.
C
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'VDATA.I'
C MACC Begin system dependent code for Macintosh
C MAC      INCLUDE 'fpc.inc'
C MACC End of system dependent code for Macintosh
      CHARACTER*40 M5IN,WM5OUT,BATIN,BLANK
CC    CHARACTER*40 M5IN,WM5OUT,BATIN,DIAGMf,BLANK
      CHARACTER*1 INTER,BATCH,M5OUT
      DIMENSION M5OUT(40)
      LOGICAL FEXIST
CC    LOGICAL CEOPEN,FEXIST
      EQUIVALENCE (WM5OUT(1:40),M5OUT)
      DATA BLANK /'                                        '/
C
C On the MAC, set the floating point control reg to control 
C rounding precision used for floating point calculations.
C MACC Begin system dependent code for Macintosh
C MAC      IMASK = getfpcontrol()
C MAC      IMASK = (IMASK.AND.RPCLEAR)
C MAC      IMASK = (IMASK.OR.RPSINGLE)
C MAC      CALL setfpcontrol(IMASK)
C MACC End of system dependent code for Macintosh
C Change for different computers
C on the MAC IOUGEN=9 and on the IBM IOUGEN=5 and on MTS IOUGEN=1
C on the MAC IOUREP=6 and on the IBM IOUREP=6 and on MTS IOUREP=7
C on the MAC IOUASK=9 and on the IBM IOUASK=0 and on MTS IOUASK=6
C on the MAC IOUERR=6 and on the IBM IOUERR=6 and on MTS IOUERR=7
C
      IOUGEN=5
      IOUREP=6
      IOUASK=0
      IOUERR=6
CC    CEOPEN=.FALSE.
      IOUPFF = 10
      IOUTTC = 11
      IOBAT=19
      FEXIST=.FALSE.
      OUTFMT=3
C
C  Initialize the name of the program.
C
      INCLUDE 'VNAME.I'
C
   10 CONTINUE
      WRITE(IOUASK,200) VERSON
  200 FORMAT (///' ',A80)
C MACC Begin system dependent code for Macintosh
C MAC      WRITE(IOUASK,201)
C MAC  201 FORMAT(/' ','Portions of this program contain',
C MAC     &       ' material copyrighted (c) by Absoft Corp 1988')
C MACC End of system dependent code for Macintosh
      WRITE(IOUASK,202)
  202 FORMAT(///'&',
     &        'Do you want to use the interactive input mode (Y/N)?')
      READ(IOUASK,100) INTER
  100 FORMAT(A1)
C
      IF(INTER.EQ.'y') INTER='Y'
      IF(INTER.EQ.'n') INTER='N'
      IF(INTER.NE.'Y'.AND.INTER.NE.'N') GOTO 10
C
C  IBM-PC Mode
C
      IF(INTER.EQ.'Y'.AND.IOUGEN.EQ.5) IOUGEN=0
C
C  MTS-Mainframe Mode
C
      IF(INTER.EQ.'Y'.AND.IOUGEN.EQ.1) IOUGEN=5
      BATCH='N'
      IF(INTER.EQ.'N') THEN
C
   15   CONTINUE
        WRITE(IOUASK,205)
  205   FORMAT(/'&','Do you want to use the batch file',
     *  ' input mode (Y/N)?')
        READ(IOUASK,100) BATCH
C
        IF(BATCH.EQ.'y') BATCH='Y'
        IF(BATCH.EQ.'n') BATCH='N'
        IF(BATCH.NE.'Y'.AND.BATCH.NE.'N') GOTO 15
      ENDIF
C
      IF(INTER.EQ.'Y') GOTO 25
      IF(BATCH.EQ.'Y') GOTO 35
C
   20 CONTINUE
      WRITE(IOUASK,210)
  210 FORMAT(/'&','Please enter the MOBILE5 input filename:')
      READ(IOUASK,110,END=99) M5IN
  110 FORMAT(A40)
C
C  2-May-1994 @ ARC-rwk Request 2-425 Exit if <return> is hit for filename
C  prompts
C
      IF(M5IN.EQ.BLANK) GOTO 99
C
      INQUIRE(FILE=M5IN,EXIST=FEXIST)
      IF(FEXIST) GOTO 25
      WRITE(IOUASK,215)
  215 FORMAT(' File does not exist, 1) if on a PC, in this directory',
     *' or'/' 2) in this account, if on a mainframe, ',
     *'please try again.'/' ')
      GOTO 20
C
   25 IF(INTER.EQ.'N'.AND.IOUGEN.EQ.9) IOUGEN=5
      IF(INTER.EQ.'N') OPEN(IOUGEN,FILE=M5IN,STATUS='OLD')
C
   30 CONTINUE
      WRITE(IOUASK,220)
  220 FORMAT(/'&','Please enter the MOBILE5 output filename:')
      READ(IOUASK,110,END=99) WM5OUT
      IF(WM5OUT.NE.BLANK) GOTO 40
C
C  2-May-1994 @ ARC-rwk Request 2-425 Exit if <return> is hit for filename
C  prompts
C
      GOTO 99
C
   35 CONTINUE
      WRITE(IOUASK,230)
  230 FORMAT(/'&','Please enter the MOBILE5 batch filename:')
      READ(IOUASK,110,END=99) BATIN
C
C  5-May-1994 @ ARC-rwk Request 2-425 Exit if <return> is hit for filename
C  prompts
C
      IF(BATIN.EQ.BLANK) GOTO 99
C
      INQUIRE(FILE=BATIN,EXIST=FEXIST)
      IF(FEXIST) GOTO 60
      WRITE(IOUASK,215)
      GOTO 35
C
   40 CONTINUE
      WRITE(IOUASK,235)
  235 FORMAT(/)
C
      OPEN(IOUREP,FILE=WM5OUT)
      IOUERR=IOUREP
      OUTFMT=3
C
C  Checking for .CSV format spreadsheet version
C
C  TMM@CSC, 16 March, 1994 Request 2-392, commented out CSV code.
C
CC      DO 55 I=1,37
CC        IF(M5OUT(I).EQ.'.'.AND.
CC     *    (M5OUT(I+1).EQ.'C'.OR.M5OUT(I+1).EQ.'c').AND.
CC     *    (M5OUT(I+2).EQ.'S'.OR.M5OUT(I+2).EQ.'s').AND.
CC     *    (M5OUT(I+3).EQ.'V'.OR.M5OUT(I+3).EQ.'v')) THEN
CC                IOUERR=8
CC                        OUTFMT=7
CC                ENDIF
CC        INQUIRE(UNIT=IOUERR,OPENED=CEOPEN)
CC        IF(IOUERR.EQ.8.AND..NOT.CEOPEN) THEN
CC   45     WRITE(IOUASK,240)
CC  240     FORMAT(/'&','Please enter the diagnostic message filename:')
CC          READ(IOUASK,110) DIAGMF
CC          IF(DIAGMF.NE.BLANK) GOTO 50
CC          WRITE(IOUASK,225)
CC  225     FORMAT(' Invalid filename, please try again.'/' ')
CC          GOTO 45
CC   50     OPEN(IOUERR,FILE=DIAGMF)
CC        ENDIF
C
CC   55 CONTINUE
   60 IF(BATCH.EQ.'Y') OPEN(IOBAT,FILE=BATIN,STATUS='OLD')
C
      IF(INTER.EQ.'Y') WRITE(IOUASK,245)
  245 FORMAT(//' ','Please enter the PROMPT flag...'/'&',
     *       '(Enter 2 for Vertical or 4 for Horizontal Prompting)')
C
      NFILE=0
      NUMDID=0
C
   65 CONTINUE
C
      IF (BATCH.EQ.'Y') THEN
        IF (INTER.EQ.'N'.AND.IOUGEN.EQ.9) IOUGEN=5
        IF (NFILE.GT.0) CLOSE(IOUGEN)
        IF (NFILE.GT.0) CLOSE(IOUREP)
   70   READ(IOBAT,120,END=93) M5IN,WM5OUT
  120   FORMAT(A40/A40)
        INQUIRE(FILE=M5IN,EXIST=FEXIST)
        IF (FEXIST.AND.WM5OUT.NE.BLANK) GOTO 75
        WRITE(IOUASK,250)
        NFILE=NFILE+1
  250   FORMAT(' Error in batch run, batch file #',I4,': was skipped.')
        GOTO 70
   75   OPEN(IOUGEN,FILE=M5IN)
        OPEN(IOUREP,FILE=WM5OUT)
        IOUERR=IOUREP
        OUTFMT=3
C
C  Checking for .CSV format spreadsheet version
C
C  TMM@CSC, 16 March, 1994 Request 2-392, commented out CSV code.
C
CC        DO 90 I=1,37
CC          IF(M5OUT(I).EQ.'.'.AND.
CC     *      (M5OUT(I+1).EQ.'C'.OR.M5OUT(I+1).EQ.'c').AND.
CC     *      (M5OUT(I+2).EQ.'S'.OR.M5OUT(I+2).EQ.'s').AND.
CC     *      (M5OUT(I+3).EQ.'V'.OR.M5OUT(I+3).EQ.'v')) THEN
CC               IOUERR=8
CC               OUTFMT=7
CC          ENDIF
CC          INQUIRE(UNIT=IOUERR,OPENED=CEOPEN)
CC          IF(IOUERR.EQ.8.AND..NOT.CEOPEN) THEN
CC   80       WRITE(IOUASK,240)
CC            READ(IOUASK,110) DIAGMF
CC            IF(DIAGMF.NE.BLANK) GOTO 85
CC            WRITE(IOUASK,225)
CC            GOTO 80
CC   85       OPEN(IOUERR,FILE=DIAGMF,STATUS='NEW')
CC          ENDIF
C
CC   90   CONTINUE
        NFILE=NFILE+1
        NUMDID=NUMDID+1
        WRITE(IOUASK,255) NFILE,M5IN
  255   FORMAT(/' ','Processing batch file #',I4,': ',A40)
      END IF
C
      M5CALL=0
      INERR=0
C         
   91 M5CALL=M5CALL+1
      CALL MOBILE(INERR,*92)
      WRITE(IOUASK,260) M5CALL,INERR
  260 FORMAT(' ','Run #',I3,'  INERR =',I2/)
      GOTO 91
C
   92 CONTINUE
      WRITE(IOUASK,265)
  265 FORMAT(' ','DRIVER calls completed.'/)
C
      IF (BATCH.EQ.'Y') GOTO 65
C
   93 CONTINUE
      IF (BATCH.EQ.'Y') THEN
        WRITE(IOUASK,270) NUMDID
        WRITE(IOUASK,275) NFILE
  270   FORMAT(//' ',I4,' Batch file(s) processed')
  275   FORMAT(//' ',I4,' Batch file(s) read')
      ENDIF
C
      IF (IOUASK.EQ.9) THEN
        WRITE(IOUASK,280)
  280   FORMAT(' Processing complete.',
     &             ' Hit Return to exit from the program.')
        READ(IOUASK,285)
  285   FORMAT(1X)
      ENDIF
C
   99 STOP
C
      INCLUDE 'COMMTS.I'
C
      END

