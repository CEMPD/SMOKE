C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE PARSEC(ICY,INERR,*)
C
C  PARSEC reads in the Parameter Section variables and validates them.
C
C  Called by MOBILE.
C
C  Calls GETALH, GETCOL, GETSC1, GETVMT , LOCAL, and WTLRUL.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ VMFLAG
C    /FLAGS4/ COLDFG
C
C  Output on return:
C
C    parameter list: ICY,INERR
C
C  Notes:
C
C  For MOBILE4.1, the COLD CO read in routine GETCOL is skipped,
C  The CAA COLD CO values are set in BD17.
C  Added call to WTLRUL (SPDFLG=3) in MOBILE4.1, D02.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS4.I'
C
C  PARSEC is called once for each scenario of the run.
C
C  Get the first required scenario record.
C
      CALL GETSC1(ICY,INERR,*99)
C
C  Read in and process the second required scenario record.  This record
C  and the parameters calculated from it supply all the local (city specific)
C  data needed for the current scenario.
C
      CALL LOCAL(ICY,INERR,*99)
C
C  Check for user supplied vmt mix to be used for this scenario only.
C
      IF(VMFLAG.EQ.2) CALL GETVMT(INERR,*99)
C
C  Check for user supplied %VMT trip len info for this scenario only.
C
      IF(SPDFLG.EQ.3) CALL WTLRUL(INERR,*99)
C
C  Here read in BIGALH parameters.
C
      CALL GETALH(INERR,*99)
C
C  If cold CO standard(s) requested, read in cold CO parameters record(s).
C  (For MOBILE4.1, these values are set in BD17.)
C
C     IF(COLDFG.GE.2) CALL GETCOL(INERR,*99)
C
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C  End of run record => take alternate return => stop run (normal termination)
C
   99 RETURN1
      END
