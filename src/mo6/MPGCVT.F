      REAL FUNCTION MPGCVT(ICV,INERR)
C
C  MPGCVT returns fuel ecomony values for composite vehicle types. Must
C  use harmonic means to compute composite MPGs.    
C
C  Called by OUTDSC
C
C  Calls EFIVT, QUITER.
C
C  Changes: (Last change first)
C
C  05 Apr 02: AIR Task Bug Fix 388, Modification of constants via QUITER.
C  03 Dec 01 AIR Task 27: New Routine
C
C  Input on call or from calls:
C
C     parameter list: ICV
C
      IMPLICIT NONE
C
      INCLUDE 'CO2COM.I'     ! MPG
      INCLUDE 'IVTYPE.I'     ! ILDGT12,ILDGT34,IHDGV,IHDDV,ILDGT
                             ! VTLDGT1,VTLDGT2,VTLDGT4,VTHDGV2B,VTHDGV8B
                             ! VTHDDV2B,VTHDDV8B,VTGASBUS,VTURBAN,VTSCHOOL,
                             ! VTLDDT12,VTLDDT34
      INCLUDE 'VMXCOM.I'
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: ICV      ! composite vehicle type
      INTEGER, INTENT(INOUT) :: INERR
C
C     Local variable
C
      INTEGER IV
      REAL    SUMD
      REAL    SUMN
C
C    Initialize
C
      MPGCVT=0.0
      SUMN=0.0
      SUMD=0.0
C
C  Begin calculations
C
      IF (ICV.EQ.ILDGT12) THEN
         DO IV=VTLDGT1,VTLDGT2
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
C
      ELSE  IF (ICV.EQ.ILDGT34) THEN
         DO IV=VTLDGT3,VTLDGT4
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
C
      ELSE  IF (ICV.EQ.IHDGV) THEN
         DO IV=VTHDGV2B,VTHDGV8B
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
         IF(MPG(VTGASBUS).GT.0.) THEN
         SUMN=SUMN+(1./MPG(VTGASBUS))*VMTMIX(VTGASBUS)
         SUMD=SUMD+VMTMIX(VTGASBUS)
         END IF
C
      ELSE IF (ICV.EQ.ILDDT) THEN
         IF(MPG(VTLDDT12).GT.0.) THEN
         SUMN=SUMN+(1./MPG(VTLDDT12))*VMTMIX(VTLDDT12)
         SUMD=SUMD+VMTMIX(VTLDDT12)
         END IF
         IF(MPG(VTLDDT34).GT.0.) THEN
         SUMN=SUMN+(1./MPG(VTLDDT34))*VMTMIX(VTLDDT34)
         SUMD=SUMD+VMTMIX(VTLDDT34)
         END IF
C
      ELSE IF (ICV.EQ.IHDDV) THEN
         DO IV=VTHDDV2B,VTHDDV8B
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
         DO IV=VTURBAN,VTSCHOOL
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
C
      ELSE  IF (ICV.EQ.ILDGT) THEN
         DO IV=VTLDGT1,VTLDGT4
         IF(MPG(IV).GT.0.) THEN
         SUMN=SUMN+(1./MPG(IV))*VMTMIX(IV)
         SUMD=SUMD+VMTMIX(IV)
         END IF
         END DO
      ELSE
         CALL QUITER(0.0,ICV,178,INERR)
      END IF
C
C  Compute harmonic mean MPG (Equation is inversted since
C  the inverse of the MPG is used above).
C
      IF(SUMN.GT.0.) MPGCVT=SUMD/SUMN
C
      RETURN
      END
