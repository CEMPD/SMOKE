C Version @(#)$Id$ $Source$ $Date$ 

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
C  ------  ----  -------------------------------------------------------
C  IOUCHG   C*1  user entered changes to output unit assignments
C  IOUDEF   I    default iounit values: IOUREP,IOUERR,IOUASK
C  IOUNEW   I    output unit numbers to be used this run
C  IOUREF   C*1  character representations of i/o unit numbers 1 - 9
CC IOUSUM   I    sum of acceptable changed i/o unit numbers **** MOVE TO (1,999)
C
C  Notes:
C
C  GETIOU was modified in Version 01 to reinitialize IOUREP/ERR/ASK.
C  GETIOU was modified to make changes for microcomputer i/o units.
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C
      CHARACTER*1 IOUREF(9),IOUCHG(3)
      DIMENSION IOUNEW(3),IOUDEF(3)
C
      DATA IOUREF/'1','2','3','4','5','6','7','8','9'/
C
C  Initialize all logical output units to the default '6'.
C
C  Change for different computers
C  on the MAC DATA IOUDEF/6,6,9/ and on the IBM DATA IOUDEF/6,6,0/
C  and on MTS DATA IOUDEF/7,7,6/
C
      DATA IOUDEF/6,6,0/
C
C  Read in prompt/mode flag and output unit changes (if any).
C
      READ(IOUGEN,100,END=99) PROMPT,IOUCHG
  100 FORMAT(I1,3A1)
c
c  LMOS Modification 04/08/93
c  Comment out potential reassignment of unit numbers.
c
cC
cC  Output unit numbers parsed here.
cC
cC  Change alphanumeric to numeric values.
cC
c      IOUSUM=0
c      DO 30 JOU=1,3
cC
cC  First set IOUDEF values as defaults, in case input number is rejected.
cC
c      IOUNEW(JOU)=IOUDEF(JOU)
cC
cC  Test for allowed values (no input units).
cC
c      DO 10 JU=1,9
c      IF(JU.EQ.IOUIMD.OR.JU.EQ.IOUGEN) GOTO 10
c      IF(IOUCHG(JOU).EQ.IOUREF(JU)) GOTO 20
c   10 CONTINUE
c      GOTO 30
cC
cC  Save allowed new value.
cC
c   20 IOUNEW(JOU)=JU
c      IOUSUM=IOUSUM+JU
c   30 CONTINUE
cC
cC  Assign output units.
cC
c      IOUREP=IOUNEW(1)
c      IOUERR=IOUNEW(2)
c      IOUASK=IOUNEW(3)
cC
c  LMOS Modification 04/08/93
c  End of commented out code
c
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
      END
