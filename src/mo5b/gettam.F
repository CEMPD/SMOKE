C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETTAM(INERR,*)
C
C  GETTAM reads in, validates and stores replacement tampering rates.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS2/ IMFLAG
C    /IOUCOM/ IOUGEN,IOUASK
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /TAMEQ1/ TAMZML,TAMDR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ZEROML   R    user entered tampering zml to be checked: ZEROML <= 1.0
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'TAMEQ1.I'
C
      IF(IMFLAG.NE.1) GOTO 10
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT('&','Enter tampering intercepts (zero mile levels)')
      IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ',' for non-I/M cases only:')
      READ(IOUGEN,100,ERR=99,END=99)
     *   (((TAMZML(ID,IVTAM,IGD,1),ID=1,8),IVTAM=1,4),IGD=1,3)
  100 FORMAT(8F8.4)
      IF(PROMPT.EQ.2) WRITE(IOUASK,220)
  220 FORMAT('&','Enter tampering slopes (deterioration rates)')
      IF(PROMPT.EQ.2) WRITE(IOUASK,210)
      READ(IOUGEN,100,ERR=99,END=99)
     *   (((TAMDR(ID,IVTAM,IGD,1),ID=1,8),IVTAM=1,4),IGD=1,3)
      GOTO 20
C
C  Read in all but cell 9 of TAMZML & TAMDR (if picked up cell 9 too, could
C  just put in READ(---) TAMZML / TAMDR, since user inputs in FTN array order.
C
   10 IF(PROMPT.EQ.2) WRITE(IOUASK,200)
      IF(PROMPT.EQ.2) WRITE(IOUASK,230)
  230 FORMAT(' ',' for (non-I/M and/or non-IM240) & (I/M and/or IM240)')
      READ(IOUGEN,100,ERR=99,END=99)
     *  ((((TAMZML(ID,IVTAM,IGD,IM),ID=1,8),IVTAM=1,4),IGD=1,3),IM=1,2)
      IF(PROMPT.EQ.2) WRITE(IOUASK,220)
      IF(PROMPT.EQ.2) WRITE(IOUASK,230)
      READ(IOUGEN,100,ERR=99,END=99)
     *  ((((TAMDR(ID,IVTAM,IGD,IM),ID=1,8),IVTAM=1,4),IGD=1,3),IM=1,2)
C
C  Move overall misfueling ZML/DR to ID=9.  Compute other misfueling ZML/DR
C  for ID=4 as (overall - filler neck).  For every disablement case,
C  insure tampering zero mile levels (intercepts) do not exceed unity.
C  There is no range check for deterioration rates.
C
   20 IM=1
   25 DO 30 IGD=1,3
      DO 30 IVTAM=1,4
      TAMZML(9,IVTAM,IGD,IM)=TAMZML(4,IVTAM,IGD,IM)
      TAMDR(9,IVTAM,IGD,IM)=TAMDR(4,IVTAM,IGD,IM)
      TAMZML(4,IVTAM,IGD,IM)=
     *    TAMZML(9,IVTAM,IGD,IM)-TAMZML(3,IVTAM,IGD,IM)
      TAMDR(4,IVTAM,IGD,IM)=TAMDR(9,IVTAM,IGD,IM)-TAMDR(3,IVTAM,IGD,IM)
      DO 30 ID=1,9
      ZEROML=TAMZML(ID,IVTAM,IGD,IM)
      IF(ZEROML.GT.1.0) CALL QUITER(ZEROML,0,50,INERR)
   30 CONTINUE
      IM=IM+1
C
C When IM=3 is the sentinel value.
C
      IF(IMFLAG.EQ.1.OR.IM.EQ.3) GOTO 35
      GOTO 25
C
   35 RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,126,INERR)
C
      RETURN1
      END
