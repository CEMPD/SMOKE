      REAL FUNCTION EFCVT(ICV,IEF,IP) 
C 
C  EFCVT returns emission factors for composite vehicle  
C  types.        
C 
C  Called by OUTDT1, OUTDT3, EFCALX and HCCALX 
C 
C  Calls EFIVT, QUITER. 
C 
C  Changes: (Last change first) 
C 
C  26 fEB 01 @EPA-BAG Bug315 placed checks for zero VMT's in IF statements 
C  31 May 00 @EPA-ddj Bug117, Changed weighting from REGMIX to VMTMIX. 
C  16 May 00 @EPA-ddj BUG87, Removed conditional on STARFG. 
C  17 Apr 00 @EPA-djb Added code to avoid divide by zero. 
C   2 Feb 00 @DynTel-ddj 1-032  Modified LDDT to include mix of LDDT12 
C            and LDDT34. 
C  03 Nov 99 @DynTel-bag 1-006  Added units parameter GPS and include file 
C            UNIT_PARAM.I 
C  25 Mar 99 @DynTel-MLA 2-701  Added REGMIX weighting for evaporative 
C            components in grams. 
C  9 Nov. 98 @DynTel-HXQ 2-680 to delete the EFIVT(28,IEF,IP)*VMTMIX(28). 
C  06 Nov 98 @Dyntel-BAG 2-000 Renamed VMTGT1 and VMTGT2 to VMTGT12 and 
C            VMTGT34, respectively 
C 
C  18 SEP 98 @DynTel-ZK  2-000 Explicit Typing 
C 
C  Input on call or from calls: 
C 
C 
C     parameter list: ICV,IEF,IP 
C     common blocks: 
C     /CVNAMS/ ILDGT12,ILDGT34,IHDGV,IHDDV,ILDGT (IVTYPE.I) 
C     /FLAGS4/ STARFG 
C     /RESUL4/ ISTARTEF 
C     /UNIT_PARAM/ GPS 
C     /VMXCOM/ REGMIX,VMTMIX,VMTGT12,VMTGT34,VMTHDD,VMTHDG,VMLDGT 
C 
C  Notes. 
C 
C  10 Mar 98 @DynTel-RJD 2-622 Created this function for Mobile6 
C            to calculate emission factors for composite vehicle  
C            types.      
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C   4 Aug 98 @DynTel-MLA 2-000  Added call to QUITER for invalid value  
C            of ICV. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'RESUL4.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'UNIT_PARAM.I' 
      INCLUDE 'VMXCOM.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: ICV      ! composite vehicle type 
      INTEGER, INTENT(IN) :: IEF      ! emission type 
      INTEGER, INTENT(IN) :: IP       ! pollutant 
C 
C     Declare external functions 
C 
      REAL EFIVT 
C 
C    Initialize 
C 
      EFCVT=0.0 
C 
C  Begin calculations 
C 
      IF (ICV.EQ.ILDGT12) THEN 
         IF(VMTGT12>0.) EFCVT=(EFIVT(2,IEF,IP)*VMTMIX(2) + 
     *          EFIVT(3,IEF,IP)*VMTMIX(3))/VMTGT12 
C 
      ELSE  IF (ICV.EQ.ILDGT34) THEN 
         IF(VMTGT34>0.) EFCVT=(EFIVT(4,IEF,IP)*VMTMIX(4) + 
     *          EFIVT(5,IEF,IP)*VMTMIX(5))/VMTGT34 
C 
      ELSE  IF (ICV.EQ.IHDGV) THEN 
         IF(VMTHDG>0.) EFCVT=(EFIVT(6,IEF,IP)*VMTMIX(6) + 
     *          EFIVT(7,IEF,IP)*VMTMIX(7) + 
     *          EFIVT(8,IEF,IP)*VMTMIX(8) + 
     *          EFIVT(9,IEF,IP)*VMTMIX(9) + 
     *          EFIVT(10,IEF,IP)*VMTMIX(10) + 
     *          EFIVT(11,IEF,IP)*VMTMIX(11) + 
     *          EFIVT(12,IEF,IP)*VMTMIX(12) + 
     *          EFIVT(13,IEF,IP)*VMTMIX(13) + 
     *          EFIVT(25,IEF,IP)*VMTMIX(25))/VMTHDG 
C 
      ELSE IF (ICV.EQ.ILDDT) THEN 
         IF(VMTLDDT>0.) EFCVT = ( EFIVT(15,IEF,IP)*VMTMIX(15) + 
     &             EFIVT(28,IEF,IP)*VMTMIX(28) )/ VMTLDDT 
C 
      ELSE IF (ICV.EQ.IHDDV) THEN 
         IF(VMTHDD>0.) EFCVT=(EFIVT(16,IEF,IP)*VMTMIX(16) + 
     *          EFIVT(17,IEF,IP)*VMTMIX(17) + 
     *          EFIVT(18,IEF,IP)*VMTMIX(18) + 
     *          EFIVT(19,IEF,IP)*VMTMIX(19) + 
     *          EFIVT(20,IEF,IP)*VMTMIX(20) + 
     *          EFIVT(21,IEF,IP)*VMTMIX(21) + 
     *          EFIVT(22,IEF,IP)*VMTMIX(22) + 
     *          EFIVT(23,IEF,IP)*VMTMIX(23) + 
     *          EFIVT(26,IEF,IP)*VMTMIX(26) + 
     *          EFIVT(27,IEF,IP)*VMTMIX(27) )/VMTHDD 
C 
      ELSE  IF (ICV.EQ.ILDGT) THEN 
         IF(VMLDGT>0.) EFCVT=(EFIVT(2,IEF,IP)*VMTMIX(2) + 
     *          EFIVT(3,IEF,IP)*VMTMIX(3) + 
     *          EFIVT(4,IEF,IP)*VMTMIX(4) + 
     *          EFIVT(5,IEF,IP)*VMTMIX(5) )/VMLDGT 
      ELSE 
         CALL QUITER(0.0,ICV,178,1) 
      END IF 
C 
      RETURN 
      END 
