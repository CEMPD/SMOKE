      SUBROUTINE GETM5IM(ITEST,RC,INERR) 
C 
C 
C  GETM5IM Reads in the MOBILE5 ASM and IM240 Credits for 
C  ratio computation.  These are read from the file ASMDATA.D 
C  The ASMDATA.D file is MANDATORY if the user selects an ASM 
C  I/M program. 
C 
C 
C 
C  Called by: GETIPR 
C 
C 
C  Calls  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  11 Oct 00 AIR Task 02 Bug 272&275: Improved read error handling. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  29 Sep 00 AIR Task 07 Removed unused INERR in parameter list 
C  02 Oct 00 @EPA-djb Bug230 Changed hardcoded I/O unit to IOUIMD. 
C  30 Aug 00  @EPA-elg   Added ASM file echo to the output. 
C 
C  Input on call: 
C 
C    argument list: ITEST, INERR 
C 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list:  INERR 
C    common blocks: 
C 
C  Common block array Subscripts: 
C 
C      /LDGIMC/  ASMRAT(3,18,25) 
C 
C  Local variable Subscripts: 
C 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  --------------------------------------------------- 
C  ASMCRD      R    Contains the MOBILE5 ASM and IM240 Credits. 
C  ASMRAT      R    Contains the Ratio of the ASM to IM240 Credits. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'IOUCOM.I' 
C
      INTEGER, EXTERNAL :: JUNIT 
C
      INTEGER    RC 
C 
      INTEGER IMY, IP,IAGE,IBLK,ITEST,IT,INERR 
      REAL  ASMCRD(7,3,18,25) 
      CHARACTER*1  JUNK 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C 
C     Open the external data file. 
C 
ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )
      M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // 'ASMDATA.D'
      
      IOUIMD = JUNIT()      
      OPEN(IOUIMD, FILE=M6FILE, STATUS='OLD',ACTION='READ', 
     *     IOSTAT=RC, ERR=97) 
C 
      WRITE(IOUREP,120) 
  120 FORMAT(/'* Reading ASM I/M Test Credits from ASMDATA.D ') 
C 
C     IBLK = 1    ASM 2525 PHASE-IN 
C     IBLK = 2    ASM 2525 FINAL 
C     IBLK = 3    ASM 5015 PHASE-IN 
C     IBLK = 4    ASM 5015 FINAL 
C     IBLK = 5    ASM 2525/5015 PHASE-IN 
C     IBLK = 6    ASM 2525/5015 FINAL 
C     IBLK = 7    IM240 0.80 / 15.0 / 2.0 
C 
C     IBLK = 7  Contains the MOBILE5 IM240 credits for 
C               final cutpoints of 0.80 / 15.0 / 2.0 HC/CO/NOx. 
C               These are the basis for the ASM ratios. 
C 
      DO IBLK = 1,7 
        READ(IOUIMD,49) JUNK 
  49    FORMAT(A1) 
        IF(JUNK=='A') JUNK='B' !just to get rid of compiler warning 
C 
        DO IMY = 1,18 
          DO IP = 1,3 
            READ (IOUIMD,50,END=16,IOSTAT=RC,ERR=98) 
     &            (ASMCRD(IBLK,IP,IMY,IAGE),IAGE=2,25) 
  50        FORMAT(24F5.4) 
          END DO 
        END DO 
      END DO 
C 
C     The ASM I/M credit for IAGE=1 is zero. It is set in BD_IM1. 
C 
  16    DO IMY = 1,18 
          DO IP = 1,3 
            DO IAGE = 2,25 
              IT = ITEST - 5 
              ASMRAT(IP,IMY,IAGE) = ASMCRD(IT,IP,IMY,IAGE) / 
     *                              ASMCRD(7,IP,IMY,IAGE) 
            END DO 
          END DO 
        END DO 
C 
      RC = 0 
      GOTO 99 
C 
C     Errors 
C 
  97  CALL QUITER(0.0,0,564,INERR) 
      GOTO 99 
C 
  98  CALL QUITER(0.0,0,563,INERR) 
C 
  99  CLOSE(IOUIMD) 
      RETURN 
      END 
