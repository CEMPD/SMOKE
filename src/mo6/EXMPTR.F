      SUBROUTINE EXMPTR(MY,IV,INERR)
C
C  EXMPTR selects the pointer to the non-methane fractions.
C
C  Called by BIGCFX.
C
C  Calls QUITER.
C
C  Changes : (Last change first)
C  
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  22 Jan  1999 @DynTel-ZK 2-000 Explicit Typing
C  25 Sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C
C  Input on call:
C
C    argument list: MY,IV
C    common blocks:
C    /NMETH2/ MYGMTH,MTHMYG
C    /REGION/ IREJN
C
C  Output on return:
C
C
C  Notes:
C
C  EXMPTR was modified for MOBILE4.1.
C
      IMPLICIT NONE 
C
      INCLUDE 'NMETH2.I'
      INCLUDE 'REGION.I'
C      
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER IG
C
      DO 10 IG=1,7
      IF(MY.LE.MYGMTH(IG,IV,IREJN)) GOTO 20
   10 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(FLOAT(MY),14,61,INERR)
C
   20 MTHMYG=IG
C
      RETURN
      END
