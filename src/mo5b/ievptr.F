C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION IEVPTR(MY,IVGAS,ICASE)
C
C  IEVPTR selects pointers in support of evaporative hot soak / diurnal HC
C  calculations.  Case 1 returns the emission standard subscript value.
C  Case 2 returns the pointer to the fraction of cars that are fuel
C  injected, given the model year.
C
C  Called by BIGTCF, CCEVRT, FUEL.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IVGAS,ICASE
C    common block:
C    /EVAPAR/ IFDS
C    /EVASTD/ MYGSTD,MAXSTD,KSTD81
C    /INJECT/ MAXFIY,MINFIY
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
C  parameter list for IEVPTR.
C
      INCLUDE 'INJECT.I'
      INCLUDE 'EVAPAR.I'
      INCLUDE 'EVASTD.I'
C
      IF(ICASE.EQ.2) GOTO 60
C
C  Select emission standard subscript value.
C
      IF(MY.GE.1981.AND.IVGAS.LE.3) GOTO 50
C
C  Untampered or uncontrolled subcase, with no fuel injection.
C
      DO 20 IG=1,MAXSTD
      IF(MY.LE.MYGSTD(IG,IVGAS)) GOTO 80
   20 CONTINUE
C
C  If execution reaches here, MY higher than MYGSTD(MAXSTD,IVGAS).  BLOCK
C  DATA 09's initialization of MYGSTD needs to be corrected.  Tell user
C  and stop execution within QUITER.
C
      CALL QUITER(0.0,1,60,MY)
C
   50 IG=KSTD81(IVGAS,IFDS)
      GOTO 80
C
C  Select pointer to fuel injection (TBI & PFI) fractions.
C
   60 IF(MY.LE.MINFIY) IG=1
      IF(MINFIY.LT.MY.AND.MY.LT.MAXFIY) IG=MY-MINFIY+1
      IF(MY.GE.MAXFIY) IG=MAXFIY-MINFIY+1
C
C  Assign index number.
C
   80 IEVPTR=IG
C
      RETURN
      END
