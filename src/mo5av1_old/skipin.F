C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE SKIPIN(NREC,INERR,*)
C
C  SKIPIN skips NREC input records in the file attached to IOUGEN.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: NREC,INERR
C    common blocks:
C    /IOUCOM/ IOUGEN
C
C  Output on return:
C
C    parameter list: none
C    common block: none
C
C  Notes:
C
C  SKIPIN was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C
      IF(NREC.LT.1.OR.NREC.GT.999)
     *    CALL QUITER(0.,NREC,11,INERR)
      IF(NREC.LT.1.OR.NREC.GT.999) GOTO 99
C
      DO 10 IREC=1,NREC
      READ(IOUGEN,100,ERR=99,END=99)
  100 FORMAT(' ')
   10 CONTINUE
C
      RETURN
C
C  NREC out of range or EOF on read => take alternate return 1 => run aborts.
C
   99 CALL QUITER(0.,0,135,INERR)
C
      RETURN1
      END
