      REAL FUNCTION BIGSC3(IP,MY,SPD) 
C 
C  BIGSC3 calculates the freeway ramp and local roadway emissions 
C  as a function of hot running LA$ emissions for 2.5 and 5 mph using 
C  the MOBILE4 algorithm. 
C         For Vehicle 1,2,3,4,5 
C 
C  Called by SCFSP. 
C 
C  Calls ISPPTR and EXP 
C 
C  Changes: (Last change first) 
C  17 Aug 99 @DynTel-ddj 2-698 Changed to a function and removed IDX 
C                from parameter list. 
C  12 May 99 @DynTel-HQ 2-698 simplify the code for new speed 
C                correction 
C  20 Apr 99 @DynTel-HQ 2-698 simplify the code for new speed 
C                correction 
C 
C  Input on call: 
C 
C    argument list: IP,MY,ISP 
C    common blocks: 
C    /SPEED1/ SK2C,SK4C 
C    /SPEED3/ SADJ 
C    /SPEED4/ IVA 
C    /SPEED5/ MYHSPC 
C    /SPEED8/ SIXT5 
C 
C  Output on return: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  A        R    coefficient "A" of new scf equation: A/X+B or EXP(A+BX+CX**2) 
C  A2       R    SADJ based coefficient "A" of new scf equation: A/X+B 
C  B        R    coefficient "B" of new scf equation: A/X+B or EXP(A+BX+CX**2) 
C  B2       R    SADJ based coefficient "B" of new scf equation: A/X+B 
C  BB       R    slope coefficient "BB" of high speed scf equ. 
C  MYR      I    model year pointer for SIXT5 array 
C  S        R    Vehicle speed, held between 2.5 and 48 mph, used in equations. 
C  SP       R    Vehicle speed for high speed equation 
C  SCADJ2   R    new scf equation computed using opmode adjusted speed value 
C  SCUNA2   R    new scf equation computed using unadjusted user inputted speed 
C 
C  Notes: 
C 
C  BIGSC3 added for MOBILE4.1 v5. 
C  BIGSC3 was modified in MOBILE4.1. 
C  BIGSC3 was modified for MOBILE5a in house version for new idle algorithm. 
C  15 July 1993 @ ARC-rwk Task 2-277 
C  16 December 1997 Dyntel@MMS 2-645 modified multiple speed to a single speed. 
C                   SPD changed from an array to a single variable. 
C  20 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the SALHCF array 
C            to distinguish between start and running emissions. No speed 
C            correction is needed for start emissions. 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'SPEED1.I' 
      INCLUDE 'SPEED3.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'SPEED8.I' 
C 
      INTEGER ISP, ISPPTR, MY, IP 
      REAL A, B, SCUNA2, SCADJ2, A2, B2, S, C,SPD  
C 
C  Set the speed. 
C 
      S=SPD 
C  For NOx MY 1977 - 1979, IV 1, 2, 3, 4, 5 ( LDGV, 
C  LDGT12, and LDGT34 ); A, B, and C are constants. 
C  The equation for NOx is EXP(A+B*X+C*X**2) / EXP(A+B*Y+C*Y**2), 
C  where X & Y are defined as in the HC & CO case. 
C  Hard Code the 1977 - 1979 NOx Coeff's. 
C 
      IF(MY.LT.1980.AND.IP.EQ.3) THEN 
        A = 0.34665 
        B =-0.026082 
        C = 4.2835E-04 
C 
        SCUNA2 = EXP(A+S*(B+C*S)) 
        SCADJ2 = EXP(A+SADJ*(B+C*SADJ)) 
      ELSE 
C 
C  Adjust ISP to work with array SK2C and SK4C 
C 
        IF (MY.LT.1980) THEN 
           ISP = 1 
        ELSE 
C 
C  New equation form is (A/X+B) / (A/Y+B), where X = inputted speed 
C  and Y = operating mode adjusted FTP speed. 
C 
           ISP = ISPPTR(MY,IP,IVA,3) 
           IF (ISP.GT.13) THEN 
              ISP = 14 
           ELSE 
              IF (IVA.GT.1) THEN 
                 ISP = ISP + 1 
              END IF 
           END IF 
        END IF 
C 
        A=SK2C(1,ISP,IP,IVA) 
        B=SK2C(2,ISP,IP,IVA) 
C 
C     speed adjust 
C 
        A2=SK2C(1,ISP,IP,IVA) 
        B2=SK2C(2,ISP,IP,IVA) 
C 
C  Determine the SCFs and apply. 
C 
        SCUNA2 = (A/S) + B 
        SCADJ2 = (A2/SADJ) + B2 
C 
      END IF 
      BIGSC3 = SCUNA2/SCADJ2 
C 
      RETURN 
      END 
