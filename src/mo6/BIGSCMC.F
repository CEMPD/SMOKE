      SUBROUTINE BIGSCMC(IP,IV,MY,IDX,ISP,SPD) 
C 
C  BIGSCMC computes SALHCF for all Motorcycles  
C 
C  Called by SCFTABMC. 
C 
C  Calls ISPPTR and EXP. 
C 
C  Changes: (Last change first) 
C 
C  12 May 99 @DynTel-HQ 2-698 Add the code for new speed 
C                correction 
C 
C  Input on call: 
C 
C    argument list: IP,IV,MY,IDX,ISP 
C    common blocks: 
C    /SPEED1/ SK1 
C    /SPEED3/ SADJ,SCUNA1 
C    /SPEED4/ IVA 
C    /SPEED6/ SALHCF 
C    /SPEED8/ SIXT5 
C 
C  Output on return: 
C 
C    common blocks: 
C    /SPEED6/ SALHCF 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  BB       R    slope coefficient "BB" of high speed scf equ. 
C  IVA2     I    used to point to the LDGV for motorcycles in SIXT5. 
C  MYR      I    model year pointer for SIXT5 array 
C  SCF1     R    old scf equation computed using opmode adjusted speed value 
C  SP       R    Vehicle speed for high speed equation 
C 
C  Notes: 
C 
C  BIGSC2 added for MOBILE4.1 v4. 
C  BIGSC2 was modified for MOBILE4.1. 
C  BIGSC2 was modified for MOBILE5a in house version for new idle algorithm. 
C  15 July 1993 @ ARC-rwk Task 2-277 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  16 December 1997 Dyntel@MMS  modified multiple speed to single speed. 
C                   SPD array changed to a single variable. 
C  20 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the SALHCF array 
C            to distinguish between start and running emissions. No speed 
C            correction is needed for start emissions. 
C 
      IMPLICIT NONE  
      INCLUDE 'SPEED1.I' 
      INCLUDE 'SPEED3.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'SPEED6.I' 
      INCLUDE 'SPEED8.I' 
C 
      INTEGER ISR, ISP, ISPPTR, MY, IP, ICOEF5, IV, IDX 
C        
      REAL SCF1, SP, BB, SPD, TMPSPD, SCF 
C 
C  No speed correction is needed for start emissions (ISR=2) 
C 
      ISR = 1 
      SALHCF(IDX,IP,IV,ISR) = 1.0 
C 
      ISP=ISPPTR(MY,IP,0,2) 
C 
C  Use MOBILE3's 5th degree polynomial equation. 
 
      TMPSPD=SPD 
C 
      IF(TMPSPD.GT.48.0) TMPSPD=48.0 
C 
      SCF1=SK1(6,ISP,IP) 
      SCF=SK1(6,ISP,IP) 
      DO 10 ICOEF5=1,5 
      SCF1=SCF1*SADJ+SK1(6-ICOEF5,ISP,IP) 
      SCF=SCF*TMPSPD+SK1(6-ICOEF5,ISP,IP) 
   10 CONTINUE 
C 
      IF(IP.NE.3) THEN 
        SCF1=EXP(SCF1) 
        SCF=EXP(SCF) 
      ENDIF 
C 
C  Original SALHCF now used for low/mid speed CF only 
C 
        SALHCF(IDX,IP,IV,ISR)=SALHCF(IDX,IP,IV,ISR) 
     *                             *(SCF/SCF1) 
C 
C  Apply high speed correction when necessary.  For HC/CO and speeds less 
C  than or equal to 55 mph, use 48 mph value.  SCF is constant in this case. 
C 
      IF(IP.EQ.3.AND.SPD.GT.48.0.OR.IP.NE.3.AND.SPD.GT.55.0) THEN 
C 
C  New equation form for High Speed(.GT.48MPH) is: 
C  Speed Correction Factor * ( 1.0 + ( BB * speed-(set spd )) 
C  SP is either 48 for NOx or 55 for HC/CO 
C  Use 1981 model year LDGV SIXT5 value for motorcyles. 
C 
        SP=55.0 
        IF(IP.EQ.3) SP=48.0 
C 
C 
        BB=(SIXT5(1,IP,1)-SALHCF(IDX,IP,IV,ISR))/(65.0-SP) 
C 
        SALHCF(IDX,IP,IV,ISR)=SALHCF(IDX,IP,IV,ISR) 
     *                                   *(1.0+(BB*(SPD-SP) )) 
C 
      ENDIF 
C 
      RETURN 
      END 
