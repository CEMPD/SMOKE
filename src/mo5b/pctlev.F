C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION PCTLEV(MYR,IV,ILEV)
C
C  PCTLEV discovers the percent of LEV's implementation
C
C  Called by BEF, CCEVRT, EFCALX, and HCCALX.
C
C  Input on call:
C
C    Parameter list: MYR,IV,ILEV
C    common blocks:
C    /BASEQ9/ LEVIMP,LDTAB
C    /LEVBLK/ LEVFLG,LEVYRS,APPLEV
C
C  Output on return:
C
C    Function: PCTLEV
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  LMY      I    Index into the LEV application years
C
C
C  Notes:
C
C  This function was developed for the BEF but can be applied
C  wherever there is need to know the percent of the LEV implementation
C  functionality checking.  It was added to job-out some of the
C  "thought" in the new MOBILE5v11 BEF.  This was added in MOBILE5v11
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C  Aug-18-1994 @ CSC-bsg Request 2-462 added checks for LEV start year and
C  LEV flag to indicate LEV's exist.
C
C  October-5-1994, @CSC pme request 2-472. Check for APPLEV flag was added 
C  to identify whether or not the LEV stds should be applied to LDGT2.
C  APPLEV=1 (do not apply stds), APPLEV=2 (apply stds)
C
      INCLUDE 'BASEQ9.I'
      INCLUDE 'LEVBLK.I'
C
      PCTLEV=0.0
      IF(LEVFLG.EQ.2.AND.IV.LE.3.AND.MYR.GE.LEVYRS) THEN
C
      LMY=(MYR-1994)+1
      IF(LMY.LT.1) LMY=1
      IF(LMY.GT.12) LMY=12
C
      IF(IV.EQ.1) THEN
         PCTLEV=LEVIMP(ILEV,LMY,1)
      ELSEIF(IV.EQ.2) THEN
         PCTLEV=(LEVIMP(ILEV,LMY,2)*LDTAB(1,1)+
     *           LEVIMP(ILEV,LMY,3)*LDTAB(1,2))
      ELSEIF(IV.EQ.3) THEN
         IF(APPLEV.EQ.2)THEN
           PCTLEV=(LEVIMP(ILEV,LMY,4)*LDTAB(2,1)+
     *             LEVIMP(ILEV,LMY,5)*LDTAB(2,2))
         ENDIF
      ENDIF
C
      ENDIF
      RETURN
      END
