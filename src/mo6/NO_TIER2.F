      SUBROUTINE NO_TIER2() 
C 
C    NO_TIER2 sets the TIER2_FLAG to zero and modifies the LEVIMP 
C    array.  The modifications made to the LEVIMP array are to set 
C    all post-2001 TIER2 fraction to zero and all post-2001 LEV 
C    fractions to one. 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  30 Dec 99 @DynTel-ddj 1-021, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list:  
C    common blocks: 
C    /FLAGS5/      TIER2_FLAG 
C    /BASEQ9/      LEVIMP, MAXLYR 
C    /IOUCOM/      IOUREP 
C    /IVTYPE/      VTLDGT1, VTLDGT2, VTLDGT3, VTLDGT4, VTLDGV 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /BASEQ9/ LEVIMP          LEVIMP(ISTD, IYR, IV) 
C    /FLAGS5/ TIER2_FLAG 
C 
C  Common block variable / array dictionary: 
C 
C   Name     Type              Description 
C  ------    ----  ------------------------------------------------------- 
C TIER2_FLAG  I    Flag indicating whether Tier2 standard applicable: 
C                    1. Tier2 standard to be applied. 
C                    2. Tier2 standard bypassed. 
C LEVIMP      R    Array containing phase-in fractions for 10 emission 
C                  standards. 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER     IYR 
      INTEGER     IV 
C 
C       Set TIER2_FLAG to zero (bypassed EF_TIER2 function) 
C 
      TIER2_FLAG = 0 
C 
      WRITE(IOUREP,100) 
  100 FORMAT(/,'  User disabled Tier2 standard.  ', 
     &   'Tier2 standards are treated as LEV standards.') 
C 
C       Set phase-in fractions for Tier2 to zero and LEV to one 
C       for post-2001 model years LDGV, LDT1 and LDT2. 
C 
      DO IV = VTLDGV, VTLDGT2 
         DO IYR = 9, MAXLYR 
            LEVIMP(4,IYR,IV) = 0.0 
            LEVIMP(8,IYR,IV) = 1.0 
         END DO 
      END DO 
C 
C       Set phase-in fractions for Tier2 to zero and Tier1 to one 
C       for post-2001 model years LDT3 and LDT4. 
C 
      DO IV = VTLDGT3, VTLDGT4 
         DO IYR = 9, MAXLYR 
            LEVIMP(4,IYR,IV) = 0.0 
            LEVIMP(3,IYR,IV) = 1.0 
         END DO 
      END DO 
C 
      RETURN 
      END 
