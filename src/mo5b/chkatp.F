C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION CHKATP(MYR,IVTAM,ICY)
C
C  CHKATP discovers the status of the Anti-Tampering Program
C  for a given MYR, Vehicle, and CY for functionality.
C
C  Called by FAIL, TAMPER.
C
C  Input on call:
C
C    Parameter list: MYR,IVTAM,ICY
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
C    /FLAGS3/ ATPFLG
C
C  Output on return:
C
C    Function: CHKATP
C    1=none
C    2=ATP applies but EVAP does not apply
C    3=ATP applies and EVAP applies
C
C  Notes:
C
C  This function was developed for the FAIL code re-write ordered as
C  part of a March 31, 1992 assignment to allow multiple
C  functionality checking.  It was added to job-out some of the
C  "thought" in the new MOBILE5.3 FAIL.  This was added in MOBILE5.3.
C  Also there is a built in temporary fix to an ATP evap inconsistency
C
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'FLAGS3.I'
C
      CHKATP = 1
C
C  Here is one possible exclusion if the vehicle is not available.
C
      IF(IVTAM.GT.4) RETURN
C
      IF(ATPFLG.EQ.1) RETURN
C
C  The following is kept as four separate IFs for clarity only
C
      IF(ICY.LT.LAPSY) RETURN
      IF(MYR.LT.LAP1ST) RETURN
      IF(MYR.GT.LAPLST) RETURN
      IF(LVTFLG(IVTAM).EQ.1) RETURN
C
      CHKATP = 2
      IF(DISTYP(6).EQ.2.AND.DISTYP(8).EQ.(2)) CHKATP = 3
C
      RETURN
      END
