C Version @(#)$Id$  $Source$ $Date$

      SUBROUTINE TRANS( IERR )

C***************************************************************************
C subroutine body starts at line 57
C 
C    This program was taken from EPS2, which has the following header 
C    information:
C
C            SYSTEMS APPLICATIONS INTERNATIONAL
C            101 LUCAS VALLEY ROAD
C            SAN RAFAEL, CA  94903
C
C-----------------------------------------------------------------------
C
C    TRANSLATE THE CHARACTER STRINGS FROM UAMCNTL TO INTEGERS.
C    STORE IN COMMON BLOCK FROM MRGUAM.INC.
C
C    ARGUMENT DECLARATION.
C     OUTPUTS:
C       IERR  I ERROR FLAG
C
C-----------------------------------------------------------------------
C    LOG:
C-----------------------------------------------------------------------
C
C      11/10/91  --HJT--  ORIGINAL DEVELOPMENT
C
C-----------------------------------------------------------------------
C    INCLUDE FILES:
C-----------------------------------------------------------------------
C
      include "param.inc"
      include "mrguam.inc"
      include "uamcntl.inc"
C
C-----------------------------------------------------------------------
C    ARGUMENT DECLARATIONS:
C-----------------------------------------------------------------------
C
      INTEGER*4 IERR
C
C-----------------------------------------------------------------------
C    LOCAL VARIABLES:
C-----------------------------------------------------------------------
C
      CHARACTER*1 TEMPCH, CHCH
      INTEGER*4 ISPC, ICHAR
C
C-----------------------------------------------------------------------
C   ENTRY POINT:
C-----------------------------------------------------------------------
C
C   --- ASSUME FAILURE ---
C
      IERR = IFAIL
C
      DO 110 ISPC = 1, NSPECS
          DO 100 ICHAR = 1, 10
              TEMPCH = SPCNAM(ISPC)(ICHAR:ICHAR)
              READ(TEMPCH, 8000) INTNAM(ISPC,ICHAR)
              WRITE(CHCH, 9000) INTNAM(ISPC,ICHAR)
  100     CONTINUE
  110 CONTINUE
C
      IERR = ISUCES
      GOTO 9999
C
C-----------------------------------------------------------------------
C   FORMAT STATEMENTS:
C-----------------------------------------------------------------------
C
 8000 FORMAT(A1)
 9000 FORMAT(A1)
C
C-----------------------------------------------------------------------
C   RETURN POINT:
C-----------------------------------------------------------------------
C
 9999 CONTINUE
      RETURN
      END
