      INTEGER FUNCTION CHKATP(MYR,IV,ICY) 
C 
C  CHKATP discovers the status of the Anti-Tampering Program 
C  for a given MYR, Vehicle, and CY for functionality. 
C 
C  Called by FAIL, TAMPER. 
C 
C  Changes: (Last change first) 
C 
C  13 Dec 00 @EPA-bag t-x6 added INTENT(IN) declaration of 
C            function parameters 
C  26 january 99 @DynTel-zk 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MYR,IV,ICY 
C    common blocks: 
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(14) 
C    /FLAGS3/ ATPFLG 
C    /VVNAMS/ VVGAS  (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
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
C  23 Apr 97 @DynTel-MLA  Added type declaration to FUNCTION statement. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: MYR 
      INTEGER IVTAM 
C 
      CHKATP = 1 
C 
C  Here is one possible exclusion if the vehicle is not available. 
C 
      IF(VVGAS(IV).EQ.0) RETURN 
C 
      IF(ATPFLG.EQ.1) RETURN 
C 
C  Convert IV to an index for LVTFLG. 
C 
      IVTAM = IV 
      IF (IV.EQ.VTGASBUS) THEN 
        IVTAM = 14 
      ENDIF 
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
