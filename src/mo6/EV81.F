      REAL FUNCTION EV81(IFDS,IVGAS,RIVAL,MY) 
C 
C  EV81 calculates a hot soak or diurnal rate for the 1981+ LDGV and 
C  LDGT emission standards, given the fuel delivery system and the 
C  fuel RVP or Wade Index, respectively. 
C 
C  Called by EVMAIN, EVSTD and EV85. 
C 
C  Changes: (Last change first) 
C 
C  24 Aug 00 @EPA-djb Bug 233, Applied Arcadis equations to fuel 
C            injected trucks as well as carbureted. 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  28 Sep  1999 @EPA-EG          Updated EVUNC function to reflect Arcadis Hot Soak Analysis. 
C  22 Jan  1999 @Dyntel-BG 2-664 Removed diurnal calculations 
C  23 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: IFDS,IVGAS,RIVAL 
C 
C    common blocks: 
C    /EVAHS1/ HSEQ 
C    /EVAPHD/ IVVMAP 
C 
C  Output on return: 
C 
C    function: EV81 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IVV      I    vehicle (IVV=1) or truck (IVV=2) classification switch 
C 
C 
C  Notes: 
C 
C  Modified for MOBILE4.1 to handle pressure/purge fail effects. 
C  Modified for MOBILE5v13 to reflect changes in diurnal. 
C     The changes account for the difference 
C     between a 1 hour measurement for diurnal and an 8 hour 
C     measurement of diurnal emissions. 
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar 
C  year removed IY and ICY from the argument list 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAHS1.I' 
      INCLUDE 'EVAPHD.I' 
C 
      INTEGER, INTENT(IN) :: IVGAS 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: MY 
      REAL, INTENT(IN)    :: RIVAL 
      INTEGER IVV 
C 
C  Map for car or truck. 
C 
      IVV=IVVMAP(IVGAS) 
C 
C  Hot soak 
C 
      IF(RIVAL.LT.9.0) THEN 
C 
C  All 1981 and newer model year vehicles get the new MOBILE6 correction. 
C  All others keep the MOBILE5 correction (HSEQ). 
C 
        IF(MY.GE.1981 .AND. MY.LE.1985) THEN 
          EV81 = HSEQ_81(IFDS,IVV,1) + 
     *           HSEQ_81(IFDS,IVV,2)*RIVAL 
        ELSEIF(MY.GE.1986) THEN 
          EV81 = HSEQ_86(IFDS,IVV,1) + 
     *           HSEQ_86(IFDS,IVV,2)*RIVAL 
        ELSE 
          EV81 = HSEQ(IFDS,IVV,1) + 
     *           HSEQ(IFDS,IVV,2)*RIVAL 
        ENDIF 
C 
      ELSEIF(RIVAL.GE.9.0) THEN 
            EV81 = HSEQ(IFDS,IVV,3) + 
     *             HSEQ(IFDS,IVV,4)*RIVAL + 
     *             HSEQ(IFDS,IVV,5)*RIVAL**2 
      ENDIF 
C 
      EV81=AMAX1(EV81,0.0) 
C 
      RETURN 
      END 
