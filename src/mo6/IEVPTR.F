      FUNCTION IEVPTR(MY,IVGAS,ICASE) 
C 
C  IEVPTR selects pointers in support of evaporative hot soak / diurnal HC 
C  calculations.  Case 1 returns the emission standard subscript value. 
C  Case 2 returns the pointer to the fraction of cars that are fuel 
C  injected, given the model year. 
C 
C  Called by BIGTCF, EVP_HRLY, FUEL. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  09 May 01 @EPA-bag Added logic for post 2004 fuel injection vehicle percentages 
C  04 Dec 00 @EPA-bag Task X6, declare MY with INTENT(IN) 
C  27 Jul 00 @EPA-djb Bug 208, removed MAXSTD. 
C   2 Jun 00 @EPA-ddj Bug 86, Changed conditional from VVLDGAS to 
C            VVGAS (HDGVs included). 
C  18 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IVGAS,ICASE 
C    common block: 
C    /EVAPAR/ IFDS 
C    /EVASTD/ MYGSTD,KSTD81 
C    /INJECT/ MAXFIY,MINFIY 
C    /VVNAMS/ VVLDGAS  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: IEVPTR 
C 
C  Notes: 
C 
C  In either Case, a default assignment is not included because there will 
C  always be a match. 
C 
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed in MOBILE4.1. 
C  In fact, all dependency on KEYEQ was removed and it was deleted from the 
C  argument list for IEVPTR. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'INJECT.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'EVAPAR.I' 
      INCLUDE 'EVASTD.I' 
      INTEGER, INTENT(IN) :: MY 
      INTEGER ICASE, IVGAS, MYQTR, IG, IEVPTR 
C 
      IF(ICASE.EQ.2) GOTO 60 
C 
C  Select emission standard subscript value. 
C 
      IF(MY.GE.1981.AND.VVGAS(IVGAS).EQ.1) GOTO 50 
C 
C  Untampered or uncontrolled subcase, with no fuel injection. 
C 
      DO 20 IG=1,5 
      IF(MY.LE.MYGSTD(IG,IVGAS)) GOTO 80 
   20 CONTINUE 
C 
C  If execution reaches here, IG index is greater than 5.  BLOCK 
C  DATA 09's initialization of MYGSTD needs to be corrected.  Tell user 
C  and stop execution within QUITER. 
C 
      MYQTR=MY 
      CALL QUITER(0.0,1,60,MYQTR) 
C 
   50 IG=KSTD81(IVGAS,IFDS) 
      GOTO 80 
C 
C  Select pointer to fuel injection (TBI & PFI) fractions. 
C 
   60 IF(MY.LE.MINFIY) IG=1 
      IF(MINFIY.LT.MY.AND.MY.LT.MAXFIY) IG=MY-MINFIY+1 
C      IF(MY.GE.MAXFIY) IG=MAXFIY-MINFIY+1 
      IF(MY>=MAXFIY.AND.MY<2004) IG=MAXFIY-MINFIY+1 
      IF(MY>=2004) IG=17 
C 
C  Assign index number. 
C 
   80 IEVPTR=IG 
C 
      RETURN 
      END 
