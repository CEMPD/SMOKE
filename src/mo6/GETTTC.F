      SUBROUTINE GETTTC(INPREC,INERR,RC) 
C 
C  GETTTC processes the I/M INCLUDE TTC labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER 
C 
C  Changes: (Last change first) 
C 
C  27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads. 
C  02 Oct 00 @EPA-djb Bug230, Removed external file I/O (IOUTTC) statements. 
C  15 Sep 00 AIR Task 08: Expanded number of I/M programs to handle EVAP 
C  04 Oct 99 @Dyntel-jwr 2-705 Commented out code that reads external  
C            tech training cred data file 
C 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C  11 Feb 99 @DynTel-MLA 2-663  This routine was substantially rewritten 
C            for Mobile6. 
C  29 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /TTCINF/  ITTCFG 
C    /IOUCOM/ 
C    /CHKLST/  MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C     /CHKLST/ IMCHECK 
C     /TTCINF/ TTCDAT,TTCSML 
C 
C  Notes: 
C 
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C  Subroutine was developed in accordance with this task request.   
C  
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'TTCINF.I' 
C 
      INTEGER    INERR 
      CHARACTER  INPREC*(*) 
      INTEGER    RC 
C 
      CHARACTER  DATAFIELD*100 
      INTEGER    IMPGM 
      INTEGER    IOS 
      INTEGER    PTR1 
      INTEGER    PTR2 
      REAL       X 
C 
      RC = 0 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,IOSTAT=IOS,ERR=98) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program. 
C 
C     Reset ITTCFG to 1 for untrained techncicians 
      IMCHECK(1,IMPGM) = 2 
      ITTCFG(IMPGM) = 1          
C 
C     Look for the name of a file containing user-supplied TTC credits. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
C 
C     If no file name was supplied, the default values for 
C     the TTC credits will be used. 
C 
      IF (DATAFIELD.EQ.' ') THEN 
        TTCSML(IMPGM) = ' ' 
        GOTO 99 
      ENDIF 
C 
  98  RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
C 
  99  RETURN 
      END 
       
