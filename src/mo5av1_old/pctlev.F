C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION PCTLEV(MYR,IV,ILEV)
C
C  PCTLEV discovers the percent of LEV's implementation
C
C  Called by BEFUN, CCEVRT, and HCCALX.
C
C  Input on call:
C
C    Parameter list: MYR,IV,ILEV
C    common blocks:
C    /BASEQ9/ LEVIMP,LDTAB
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
C  This function was developed for the BEFUN but can be applied
C  wherever there is need to know the percent of the LEV implementation
C  functionality checking.  It was added to job-out some of the
C  "thought" in the new MOBILE5v11 BEFUN.  This was added in MOBILE5v11
C
C
      REAL LEVSTD,LEVIMP,LDTAB
C
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
C
      IF(IV.GT.2) RETURN
C
      LMY=(MYR-1994)+1
      IF(LMY.LT.1) LMY=1
      IF(LMY.GT.12) LMY=12
      PCTLEV=0.0
C
      IF(IV.EQ.1) PCTLEV=LEVIMP(ILEV,LMY,1)
      IF(IV.EQ.2) PCTLEV=(LEVIMP(ILEV,LMY,2)*LDTAB(1,1)+
     *                    LEVIMP(ILEV,LMY,3)*LDTAB(1,2))
C
      RETURN
      END
