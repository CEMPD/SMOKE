      SUBROUTINE SCFTABMC(IP,IV,MY,IDX,ISR,INERR)
C
C SCFTABMC calculate the speed correction factor for speed
C 2.5-65 from old mobile5 method
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  18 May 1999 @DynTel-HQ 2-698  Add IMPLICIT NONE statement
C 12  May  1999 @Dyntel-HQ 2-698 simplify the code for caculation of SCF
C  1 March 1999 @DynTel-HQ 2-698 calculate SFC for 2.5-65 m/g
C                          
C                          
C    Called By HRLOOP
C
C    Calls  BIGSCMC.
C
C    argument list: IP,IV,MY,IDX,ISR
C    common blocks:
C    /SPEED6/  SALHCF
C    /SPEED9/  SCFL SCFAR M5SCF 
C
      IMPLICIT NONE
      INCLUDE 'SPEED6.I'
      INCLUDE 'SPEED9.I'
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER I, ISP, IDX, ISR
      REAL SP
C
      SCFL=1.0
      SCFR=1.0
      DO I=1,14
        M5SCF(I,IP)=1.0
      END DO
C
C Calculate the SCF for local road
C
      CALL BIGSCMC(IP,IV,MY,IDX,ISP,12.9,INERR)
      IF(INERR.GT.0) GOTO 99
      SCFL=SALHCF(IDX,IP,IV,ISR)
C
C Calculate the SCF for Ramp
C
      CALL BIGSCMC(IP,IV,MY,IDX,ISP,34.9,INERR)
      IF(INERR.GT.0) GOTO 99
      SCFR=SALHCF(IDX,IP,IV,ISR)
C
C Calculate the SCF For Freeway And AC
C
      CALL BIGSCMC(IP,IV,MY,IDX,ISP,2.5,INERR)
      IF(INERR.GT.0) GOTO 99
      M5SCF(1,IP)=SALHCF(IDX,IP,IV,ISR)
      SP=0.0
      DO I=2,14
      SP=SP+5.0
      CALL BIGSCMC(IP,IV,MY,IDX,ISP,SP,INERR)
      IF(INERR.GT.0) GOTO 99
      M5SCF(I,IP)=SALHCF(IDX,IP,IV,ISR)
      END DO
C
C
C
   99 RETURN
      END
