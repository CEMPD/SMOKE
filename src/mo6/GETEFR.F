      SUBROUTINE GETEFR(INPREC,INERR,RC) 
C 
C  GETEFR reads the I/M EFFECTIVENESS labeled input record, which contains 
C    the user-supplied values for the effectivenes of test-and-repair I/M 
C    programs (three values). 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  14 Sep 00 @EPA-djb Bug 251, Eliminated read of ATPUEF and PPGUEF. 
C  26 Jun 00 @EPA-elg   Changed fixed read format to free format for 
C                       Exhaust I/M discount values. 
C  18 Nov 98 @DynTel-MLA 2-663  The old GETEFR routine was broken into 
C            three parts. Som initializations were moved to INIFLG. The 
C            input portions remain here. Post-processing was moved to 
C            the new CALEFR routine, after input processing for I/M, 
C            ATP, and purge/pressure programs is complete. 
C  08 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C  
C  Input on call: 
C 
C    argument list : INPREC 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list : INERR,RC 
C    common blocks: 
C    /IMPAR5/  DISCNT 
C    /IMPAR6/  IMEFFR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  ---------------------------------------------------------- 
C  IP       I    Pollutant index. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
C 
      CHARACTER  INPREC*(*) 
      INTEGER    IOS 
      INTEGER    RC 
      INTEGER    INERR 
      INTEGER    IP  
C 
C     The old GETEFR routine was broken into three parts. Some of the 
C     initializations were moved to INIFLG. Some post-processing has 
C     been moved to the new CALEFR routine, called by CHKINP. 
C 
      RC = 0 
C 
C     Read the user-supplied effectiveness rates. The effectiveness 
C     rate record should contain either three or five values. The first  
C     three are the rates for HC, CO, and NOx for test-and-repair I/M  
C     programs (use the same values for both computerized and manual  
C     I/M programs). The next two values are optional. They are the  
C     effectiveness rates for test-and-repair ATP and purge/pressure  
C     programs. (In prompting mode, ask the user to supply these values.) 
C 
      IMEFFR = 2 
      READ(INPREC(21:),FMT=*,ERR=95,END=96,IOSTAT=IOS) 
     &  (DISCNT(IP,2),IP=1,3) 
C 
C     Use the same values for manual and computerized test-and-repair 
C     I/M programs. 
C 
      DISCNT(1,3) = DISCNT(1,2) 
      DISCNT(2,3) = DISCNT(2,2) 
      DISCNT(3,3) = DISCNT(3,2) 
C 
C     Make sure all values entered by the user are within range. 
C 
      DO 20 IP = 1,3 
        IF (DISCNT(IP,2).LT.0.0.OR.DISCNT(IP,2).GT.1.0) THEN 
          CALL QUITER(DISCNT(IP,2), 0, 157, INERR) 
          RC = 1 
        ENDIF 
  20  CONTINUE 
C 
   90 RETURN 
C             
C  This is the exit for general FORTRAN I/O errors. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
      GOTO 99 
C 
   96 RC = -1 
      GOTO 99 
C 
   99 RETURN 
      END      
