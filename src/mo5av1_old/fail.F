C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION FAIL(MYR,KDX,IV,IK)
C
C  FAIL handles the effect of Purge and Pressure failures as well
C  as evap system tampering, returning a fraction for
C  purge/pressure.  Using old code for evaluation years less than
C  1992 and new code to fix overlaps for evaluation years greater
C  than 1991.
C
C  Called by HCCALX.
C
C  Calls FAIL1, and FAIL2.
C
C  Input on call:
C
C    parameter list: MYR,KDX,IV,IK
C    common blocks:
C    /MAXIMA/ MAXYRS
C
C  Output on return:
C
C    function: FAIL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  ICY      I    See top of program.
C  IDX      I    See top of program.
C  JDX      I    See top of program.
C  MY       I    See top of program.
C
C
C  Notes.
C
C  FAIL was added for MOBILE4.1.
C  FAIL was modified in MOBILE5 to include ZML/DR.
C  FAIL was totally rewritten for MOBILE5 to call CHKFUN and
C  CHKATP. CHKATP indicates that ATP either exists or exists
C  without EVAP: this critical information is only temporary.
C  FAIL was modified to fix a bug in release MOBILE5, that was
C  causing evaporative emissions to go up when newer vehicles
C  were in use. This problem had been identified on
C  December 24, 1992.
C  FAIL was modified again for MOBILE5 corrected.  It was turned
C  into a driver subroutine and returns the value of FAIL1 or
C  FAIL2.
C
C
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
C
C  Assign values to MY, IDX, JDX, and ICY to decide based on
C  the value of ICY which function to return the value for.
C
      MY=MYR
      IDX=KDX
      JDX=MAXYRS+1-IDX
      ICY=JDX-1+MY
      IF(ICY.GT.1991) FAIL=FAIL1(MY,IDX,IV,IK)
      IF(ICY.LE.1991) FAIL=FAIL2(MY,IDX,IV,IK)
C
      RETURN
      END
