      FUNCTION PCTLEV(MYR,IV,ILEV) 
C 
C  PCTLEV discovers the percent of LEV's implementation 
C 
C  Called by BEF, EVP_HRLY, EFCALX, and HCCALX. 
C 
C  Changes: (Last change first) 
C 
C  22 Jun 99 @DynTel-MLA 2-694  Removed unused LEV variables. This function 
C            always loops over all 5 IV types for all model years (94+). 
C  18 dec 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    Constant: MAXLYR 
C    Argument list: MYR,IV,ILEV 
C    common blocks: 
C    /BASEQ9/ LEVIMP,MAXLYR 
C    /VVNAMS/ VVLDGAS, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C    /VVTYPS/ VTLDGV  (IVTYPE.I) 
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
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
C  10 Oct 97 DynTel@MLA 2-622  Removed LDTAB weighting for LDGT12/34. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'IVTYPE.I' 
C 
      REAL PCTLEV 
      INTEGER  MYR, LMY, IV, ILEV 
C 
      PCTLEV=0.0 
C 
      IF (MYR.GE.1994.AND.VVLDGAS(IV).EQ.1) THEN 
        LMY=(MYR-1994)+1 
        IF(LMY.GT.MAXLYR) LMY=MAXLYR 
        PCTLEV = LEVIMP(ILEV,LMY,IV) 
      ENDIF 
C    
      RETURN 
      END 
