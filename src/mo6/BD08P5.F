      BLOCK DATA BD08P5 
C 
C  BLOCK DATA Subprogram 08: /IDLCOM/,/BSIDL/ 
C 
C  BLOCK DATA 08 initializes the base idle emission factors 
C  BLOCK DATA 08 also initializes deterioration factors for the idle 
C                emission facotrs. 
C 
C Changes: 
C 
C 07 May 01 AIR Task 11: New Block Data for MOBILE6.1. 
C 14 Nov 00 AIR Task 06: Added NIDLYR to common BSIDL2. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Oct-19-1993 @ ARC-bsg Subtask 322 must make the same on all platforms 
C  during the standardization caught this initialization error that was 
C  not identical on all platfoms. 
C 
C 
C  Common block array subscripts: 
C 
C  BASIDL(MAXYRS,MAXVEH)  - BASIDL(JDX,IV) 
C  EMSIDL(MAXVP1)      - EMSIDL(IV+1) 
C 
C  Common block dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /BSIDL/: 
C 
C  BASIDL   R    base idle emission rates, in g/hr 
C 
C  /BSIDL2/: 
C 
C  IDLYRS   I    base idle emission years 
C  NIDLYR   I    number of base idle years 
C 
C  /IDLCOM/: 
C 
C  IDLEPMFLG L   flag to print idle PM factors 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IDLE.I' 
      INCLUDE 'PARAME.I' 
C 
C  Oct-19-1993 @ ARC-bsg Subtask 322 must make the same on all platforms 
C  during the standardization caught this initialization error that was 
C  not identical on all platfoms. 
C 
      DATA IDLYRS /1987, 1990, 1993, NEYR/ 
      DATA NIDLYR /4/ 
C 
      DATA BASIDL / 
C     2BHDDV 
     *  5.37, 3.174, 1.86, 1.004, 
C 
C     LHDDV 
     *  5.37, 3.174, 1.86, 1.004, 
C 
C     MHDDV 
     *  5.37, 3.174, 1.86, 1.004, 
C 
C     HHDDV 
     *  5.37, 3.174, 1.86, 1.004, 
C 
C     BUSES 
     *  5.37, 3.174, 1.86, 1.004/ 
C 
      DATA IDLEPMFLG/.FALSE./ 
C 
      END 
