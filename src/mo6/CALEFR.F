      SUBROUTINE CALEFR(INERR,RC) 
C 
C  CALEFR completes processing for the user-supplied effectiveness rates 
C    that were specified on the I/M EFFECTIVENESS labeled input record. 
C    This processing was deferred until the I/M, ATP, and purge/pressure 
C    inputs were complete. 
C 
C  Called by CHKINP. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  08 Jun 00 AIR Task 03: Removed non-block DO constucts. 
C  10 Feb 99 @DynTel-MLA 2-663  This code was moved from GETEFR in Mobile6. 
C 
C  
C  Input on call: 
C 
C    argument list : INERR 
C    common blocks: 
C     /ATPAR4/  EFF, EFFSAV 
C     /FLAGS2/  IMFLAG 
C     /FLAGS3/  ATPFLG 
C     /IMPAR2/  ITEST 
C     /IMPAR5/  DISCNT, ATPUEF, PPGUEF 
C     /IMPAR6/  INTYP 
C     /PRGCH2/  PRGFLG 
C     /PRSCH2/  PRSFLG 
C 
C  Output on return: 
C 
C    argument list : INERR,RC 
C    common blocks: 
C    /ATPAR4/  EFF 
C    /IMPAR5/  ATPUEF, PPGUEF, DISCNT 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  ---------------------------------------------------------- 
C  ATPADJ   R    Adjustment to ATP effectiveness due to presence of IM240. 
C  AVGIMF   R    Average of IM effectiveness rates = (HC + CO + NOx) / 3.0 
C  EFFREC   C    Character representation of Effectiveness Rec. (The actual 
C                record is variable length - some fields are optional.) 
C  IATP     I    One dimension of EFF array - see BLOCK DATA 25. 
C  ICOL     I    One dimension of EFF array - see BLOCK DATA 25. 
C  IROW     I    One dimension of EFF array - see BLOCK DATA 25. 
C 
C  Notes: 
C 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C 
      IMPLICIT NONE 
      INCLUDE 'ATPAR4.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'PRGCH2.I' 
      INCLUDE 'PRSCH2.I' 
C 
      REAL       ATPADJ 
      REAL       AVGIMF 
      INTEGER    IATP 
      INTEGER    ICOL 
      INTEGER    INERR 
      INTEGER    IROW 
      INTEGER    RC 
C         
C     Initializations for ATPUEF and PPGUEF were moved to INIFLG. 
C 
      RC = 0 
C 
C     The user-entered value for ATPUEF will be ignored unless an 
C     ATP program is explicitly modeled. The value for PPGUEF will 
C     be ignored unless a Purge or Pressure program is modeled. 
C 
      IF (ATPFLG.EQ.1) ATPUEF = 0.0 
C 
      IF (PRGFLG.EQ.1.AND.PRSFLG.EQ.1) PPGUEF = 0.50 
C 
C     Interpolate betwen the default values for test-only and test-and- 
C     repair programs, based on the value supplied by the user. If the  
C     user specified a value of 0.0 for ATPUEF (or if the user supplied 
C     no value and ATPUEF defaults to 0.0) the values we place in EFF 
C     will be the default values for test-and-repair ATP programs. First, 
C     make sure ATPUEF is between 0.0 and 1.0, inclusive. 
C 
      IF (ATPUEF.LT.0.0.OR.ATPUEF.GT.1.0) THEN 
        CALL QUITER(ATPUEF, 0, 167, INERR) 
        GOTO 90 
      ENDIF 
C 
C     Reset EFF (to default values if no alternate rates were entered). 
C 
      DO ICOL = 1,3 
        DO IATP = 1,2 
          DO IROW = 1,15 
            EFF(2,ICOL,IATP,IROW) = EFFSAV(2,ICOL,IATP,IROW) + ATPUEF* 
     & (EFFSAV(1,ICOL,IATP,IROW)-EFFSAV(2,ICOL,IATP,IROW)) 
          END DO 
        END DO 
      END DO 
C 
C     IM240 programs are required to cover three ATP inspection types. If  
C     the user is claiming greater than default effectiveness for their 
C     test-and-repair IM240 program, the effectiveness of these three ATP 
C     inspection types must be enhanced as well. Map the enhanced IM240 
C     effectiveness, which varies from 0.50 to 1.00 for three different 
C     pollutants, into enhanced effectiveness for ATP programs, which  
C     varies from 0.00 to 1.00. Skip the adjustment if the I/M program 
C     effectiveness is no greater than the default value, or if the I/M 
C     program effectiveness is less than the claimed ATP effectiveness, 
C     or if there are no test-and-repair IM240 programs.  

C 
      AVGIMF = (DISCNT(1,2) + DISCNT(2,2) + DISCNT(3,2)) / 3.0 
C 
      IF (AVGIMF.LE.0.50) GOTO 90 
C 
      ATPADJ = 2.0 * (AVGIMF - 0.50) 
C 
      IF (ATPADJ.LE.ATPUEF) GOTO 90 
C 
      IF (ITEST(1).EQ.4.AND.INTYP(1).GT.1) GOTO 32 
C 
      IF (IMFLAG.NE.3.AND.IMFLAG.NE.5) GOTO 90 
C 
      IF (ITEST(2).NE.4.OR.INTYP(2).EQ.1) GOTO 90 
C 
C     Adjust rows 3, 7, and 9 in EFF. 
C 
   32 CONTINUE 
C 
      DO ICOL = 1, 3 
C 
        DO IATP = 1, 2 
          EFF(2, ICOL, IATP, 3) = EFFSAV(2, ICOL, IATP, 3) + 
     & ATPADJ * (EFFSAV(1,ICOL,IATP,3)-EFFSAV(2,ICOL,IATP,3)) 
          EFF(2, ICOL, IATP, 7) = EFFSAV(2, ICOL, IATP, 7) + 
     & ATPADJ * (EFFSAV(1,ICOL,IATP,7)-EFFSAV(2,ICOL,IATP,7)) 
          EFF(2, ICOL, IATP, 9) = EFFSAV(2, ICOL, IATP, 9) + 
     & ATPADJ * (EFFSAV(1,ICOL,IATP,9)-EFFSAV(2,ICOL,IATP,9)) 
        END DO 
      END DO 
C 
C     Convert the 'effectiveness' rate, as entered by the 
C     user, to a 'discount' rate, as expected by the code. 
C 
      PPGUEF = 1.0 - PPGUEF 
C 
   90 RETURN 
      END      
