      SUBROUTINE SCFSLOW(IP,IV,MY,INERR)
C
C MSCFCO calculate the speed correction factor for speed
C 2.5 and 5 from old mobile5 method and modified to new method.
C
C Changes : (Last Change First)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C 17 Aug 99 @DynTel-ddj 2-698  Changed SCFSP from a subroutine CALL to
C               a function call.  Changed HC calculation from THC to
C               NMHC.
C 18 May 1999 @DynTel-HQ 2-698  Add IMPLICIT NONE statement
C  1 March 1999 @DynTel-HQ 2-698 calculate SCF for 2.5 and 5 m/h
C                          
C    Called By HRLOOP
C
C    Calls SCFSP.
C
C    argument list: IP,IV,MY,IDX,ISR
C    common blocks:
C    /SPEED6/    SALHCF
C    /SPEED9/  SCFF SCFAC SCFF7 SCFAC7
C
      IMPLICIT NONE
      INCLUDE 'SPEED6.I'
      INCLUDE 'SPEED9.I'
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER J
      INTEGER IPI
      REAL SCF2, SCF5, SCF7
      REAL SCFSP
C
C calculate the SCF for 2.5, 5, 7.1, for MOBILE5 method
C
      SCF2 = SCFSP(IP,IV,MY,2.5,INERR)
      SCF5 = SCFSP(IP,IV,MY,5.0,INERR)
      SCF7 = SCFSP(IP,IV,MY,7.1,INERR)
C
C calculate the new SCF
C
      IF (IP.EQ.1) THEN
         IPI = 4
      ELSE
         IPI = IP
      END IF
C
      DO J = 1, 3
         SCFF(J,1,IPI)  = SCF2 + SCFF7(J,IPI)  - SCF7
         SCFF(J,2,IPI)  = SCF5 + SCFF7(J,IPI)  - SCF7
         SCFAC(J,1,IPI) = SCF2 + SCFAC7(J,IPI) - SCF7
         SCFAC(J,2,IPI) = SCF5 + SCFAC7(J,IPI) - SCF7
      END DO
      RETURN
      END
