C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETIOU(IOUSUM,*)
C
C  GETIOU reads in the prompt / mode flag and, if present, user supplied
C  output unit numbers for the run.  GETIOU then identifies the numbers and
C  for each number determines if it can be used for the requested purpose and,
C  if so, assigns the number to the appropriate unit mnemonic.
C
C  Called by CONSEC.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C
C  Output on return:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUREP,IOUERR,IOUASK
C
C  Local array subscripts:
C
C  IOUCHG(9)  -  IOUCHG ( JU )
C  IOUDEF(3)  -  IOUDEF ( JOU )
C  IOUNEW(3)  -  IOUNEW ( JOU )
C  IOUREF(9)  -  IOUREF ( JU )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  IOUCHG   C*1  user entered changes to output unit assignments
C  IOUNEW   I    default iounit values: IOUREP,IOUERR,IOUASK and
C                output unit numbers to be used this run
C  IOUREF   C*1  character representations of i/o unit numbers 1 - 9
C  IOUSUM   I    sum of acceptable changed i/o unit numbers **** MOVE TO (1,999)
C
C
C  Notes:
C
C  For MOBILE4.1
C  GETIOU was modified to reinitialize IOUREP/ERR/ASK.
C  For MOBILE5
C  GETIOU was modified to make changes for microcomputer i/o
C  units.
C  March 1, 1993, GETIOU was modified to prevent over-riding of
C  IOUERR if the IOUCHG option was not used and the CSV option
C  was specified.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
C
      CHARACTER*1 IOUREF(9),IOUCHG(3)
      DIMENSION IOUNEW(3)
C
      DATA IOUREF/'1','2','3','4','5','6','7','8','9'/
C
C  Initialize all logical output units to the defaults.
C
      IOUNEW(1)=IOUREP
      IOUNEW(2)=IOUERR
      IOUNEW(3)=IOUASK
C
C  Read in prompt/mode flag and output unit changes (if any).
C
      READ(IOUGEN,100,END=99) PROMPT,IOUCHG
  100 FORMAT(I1,3A1)
C
C  Output unit numbers parsed here.
C
C  Change alphanumeric to numeric values.
C
   90 IOUSUM=0
      DO 30 JOU=1,3
C
C  Test for allowed values (no input units).
C
        DO 10 JU=1,9
          IF(JU.EQ.IOUIMD.OR.JU.EQ.IOUGEN) GOTO 10
          IF(IOUCHG(JOU).EQ.IOUREF(JU)) GOTO 20
   10   CONTINUE
        GOTO 30
C
C  Save allowed new value.
C
   20   IOUNEW(JOU)=JU
        IOUSUM=IOUSUM+JU
   30 CONTINUE
C
C  Assign output units.
C
      IOUREP=IOUNEW(1)
      IOUERR=IOUNEW(2)
      IOUASK=IOUNEW(3)
C
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
      END
