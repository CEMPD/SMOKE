      REAL FUNCTION ATEFCVT(ICV,IXV,IAT,INERR)
C
C  ATEFCVT returns AIR TOXIC emission factors for composite vehicle 
C  types.       
C
C  Called by ATOUTDSC
C
C  Calls QUITER.
C
C  Changes: (Last change first)
C
C  05 Apr 02: AIR Task Bug Fix 388, Modification of constants via QUITER.
C  09 Mar 02: AIR Task 28: Added USE ATHEAP due to expanded ADD HAPS changes.
C  20 Sep 01: AIR Tasks 15-21: Comments
C  17 Aug 01: AIR Tasks 15-21: AT version of EFCVT
C
C     Parameter list: ICV,IXV,IAT
C
C  Notes.
C
C
      USE ATHEAP, ONLY: EFAT
C
      IMPLICIT NONE
C
      INCLUDE 'IVTYPE.I'   ! ILDGT12,ILDGT34,IHDGV,IHDDV,ILDGT (IVTYPE.I)
      INCLUDE 'VMXCOM.I'   ! REGMIX,VMTMIX,VMTGT12,VMTGT34,VMTHDD,VMTHDG,VMLDGT
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: IAT      ! AT pollutant
      INTEGER, INTENT(IN) :: ICV      ! composite vehicle type
      INTEGER, INTENT(IN) :: IXV      ! Exhaust (1), Evap Type (2-6)
      INTEGER, INTENT(INOUT) :: INERR
C
C    Initialize
C
      ATEFCVT=0.0
C
C  Begin calculations
C
      IF (ICV.EQ.ILDGT12) THEN
         IF(VMTGT12>0.) 
     *     ATEFCVT=(EFAT(IXV,IAT,2)*VMTMIX(2) 
     *            + EFAT(IXV,IAT,3)*VMTMIX(3))/VMTGT12
C
      ELSE  IF (ICV.EQ.ILDGT34) THEN
         IF(VMTGT34>0.) 
     *     ATEFCVT=(EFAT(IXV,IAT,4)*VMTMIX(4)
     *            + EFAT(IXV,IAT,5)*VMTMIX(5))/VMTGT34
C
      ELSE  IF (ICV.EQ.IHDGV) THEN
         IF(VMTHDG>0.) 
     *     ATEFCVT=(EFAT(IXV,IAT, 6)*VMTMIX( 6)
     *            + EFAT(IXV,IAT, 7)*VMTMIX( 7)
     *            + EFAT(IXV,IAT, 8)*VMTMIX( 8)
     *            + EFAT(IXV,IAT, 9)*VMTMIX( 9)
     *            + EFAT(IXV,IAT,10)*VMTMIX(10) 
     *            + EFAT(IXV,IAT,11)*VMTMIX(11)
     *            + EFAT(IXV,IAT,12)*VMTMIX(12)
     *            + EFAT(IXV,IAT,13)*VMTMIX(13)
     *            + EFAT(IXV,IAT,25)*VMTMIX(25))/VMTHDG
C
      ELSE IF (ICV.EQ.ILDDT) THEN
         IF(VMTLDDT>0.) 
     &     ATEFCVT=(EFAT(IXV,IAT,15)*VMTMIX(15)
     &             +EFAT(IXV,IAT,28)*VMTMIX(28))/VMTLDDT
C
      ELSE IF (ICV.EQ.IHDDV) THEN
         IF(VMTHDD>0.) 
     *     ATEFCVT=(EFAT(IXV,IAT,16)*VMTMIX(16)
     *            + EFAT(IXV,IAT,17)*VMTMIX(17)
     *            + EFAT(IXV,IAT,18)*VMTMIX(18)
     *            + EFAT(IXV,IAT,19)*VMTMIX(19)
     *            + EFAT(IXV,IAT,20)*VMTMIX(20)
     *            + EFAT(IXV,IAT,21)*VMTMIX(21)
     *            + EFAT(IXV,IAT,22)*VMTMIX(22)
     *            + EFAT(IXV,IAT,23)*VMTMIX(23)
     *            + EFAT(IXV,IAT,26)*VMTMIX(26)
     *            + EFAT(IXV,IAT,27)*VMTMIX(27))/VMTHDD
C
      ELSE  IF (ICV.EQ.ILDGT) THEN
         IF(VMLDGT>0.) 
     *     ATEFCVT=(EFAT(IXV,IAT,2)*VMTMIX(2)
     *            + EFAT(IXV,IAT,3)*VMTMIX(3)
     *            + EFAT(IXV,IAT,4)*VMTMIX(4)
     *            + EFAT(IXV,IAT,5)*VMTMIX(5))/VMLDGT
      ELSE
         CALL QUITER(0.0,ICV,178,INERR)
      END IF
C
      RETURN
      END
