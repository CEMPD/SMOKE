      REAL FUNCTION BIGSC1(IP,IV,SPD) 
C 
C  BIGSC1 computes SALHCF for IVB (Group B IV): HDGV, LDDV, LDDT, and HDDV. 
C  It also computes the Idle Speed Correction Factor Adjustment for this group. 
C 
C  Called by SCFTAB 
C 
C  Calls EXP. 
C 
C  Changes: (Last change first) 
C 
C  25 Sep 00 AIR Task 07: Bug 256, added initialization of function to 1.0. 
C  08 Jun 00 AIR Task 03: Removed type declarations for intrinsic functions. 
C  17 Aug 99 @DynTel-ddj 2-698  Changed from a subroutine to a function 
C                and simplified calculation for SPDGPB when 
C                VVHEAVY(IV) = 1. 
C  20 Apr 99 @DynTel-HQ 2-698 simplify the code for new speed 
C                correction 
C                           
C  Input on call: 
C 
C    argument list: IP,IV,IDX 
C    common blocks: 
C    /SPEED3/ SADJ 
C    /SPEED4/ IVB 
C    /SPEED7/ GPBSCO 
C    /VVNAMS/ VVHDGV, VVHEAVY, VVLDDSL  (IVTYPE.I) 
C 
C  Output on return: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  SB       R    vehicle speed from SPD, given IV, for IVB scf calculation 
C  SPDGPB   R    accumulator: stores scf calculations for subgroup IVB 
C 
C  Notes: 
C 
C  BIGSC1 added for MOBILE4.1 v4. 
C  BIGSC1 was modified for MOBILE4.1. 
C  BIGSC1 was modified for MOBILE5a in house version for new idle algorithm. 
C  15 July 1993 @ ARC-rwk Task 2-277 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  16 December 1997 Dyntel@MMS 2-645 modified multiple speed  to single speed. 
C                   SPD array reduced to a single vriable. 
C  20 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the SALHCF array 
C            to distinguish between start and running emissions. No speed 
C            correction is needed for start emissions. 
      IMPLICIT NONE  
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'SPEED3.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'SPEED7.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IP 
      REAL,    INTENT(IN) :: SPD 
C 
C     Declare local variables 
C 
      REAL SB 
      REAL SPDGPB 
C 
      SB=SPD 
      SPDGPB=1.0 
C 
C  Heavy Duty Vehicle classes 
C 
      IF(VVHEAVY(IV).EQ.1) THEN 
          SPDGPB = GPBSCO(1,IP,IVB) + SB*(GPBSCO(2,IP,IVB) + 
     &                                SB*GPBSCO(3,IP,IVB)) 
C 
C  Light Duty Diesel classes 
C 
      ELSEIF(VVLDDSL(IV).EQ.1) THEN 
        SPDGPB=GPBSCO(3,IP,IVB)*(SB**2-SADJ**2) 
     *        +GPBSCO(2,IP,IVB)*(SB-SADJ) 
      ENDIF 
C 
C  Unless HDDV and NOX, exponentiate the 2nd degree polynomial result. 
C 
      IF(VVHDGV(IV).EQ.0.OR.IP.NE.3) SPDGPB=EXP(SPDGPB) 
C 
        BIGSC1 = SPDGPB 
C 
      RETURN 
      END 
