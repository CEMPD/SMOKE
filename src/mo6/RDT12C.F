      SUBROUTINE RDT12C(INERR,RC) 
C 
C  RDT12C reads in the TECH1&2 I/M credits if the altitude is "new" and 
C  there is an I/M program. 
C 
C  Called by PRCSCN. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  02 Oct 00 @EPA-djb Bug230 Replaced hardcoded I/O unit with IOUIMD. 
C  11 Apr 00 AIR Task 02: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  22 Jun 99 @DynTel-MLA 2-000  Corrected a comment. 
C  21 Dec 1998 @DynTel-ZK 2-000 Explicit Typing  
C  19 Nov 98 @DynTel-MLA 2-663  Because the RSDFLG has become an array 
C            rather than a simple integer, and because the application 
C            of RSD to TECH12 vehicles in Mobile6 is uncertain at this 
C            point, all references to RSDFLG have been removed from the 
C            code which interprets the smart line. 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list:  
C    common blocks: 
C    /FLAGS2/ IMFLAG 
C    /IMPAR6/ IFILET 
C    /IMPAR7/ TNAME,IREAD,SVALTC 
C    /IOUCOM/ IOUOUT 
C    /REGION/ IREJN 
C    /RSDINF/ RSDFLG 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /IMPAR4/ SMRTL1  
C    /IMPAR7/ IREAD,SVALTC 
C    /IM12HC/ CR12HC 
C    /IM12CO/ CR12CO 
C    /REGION/ IREJN 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  BUFFER    C    To keep "Smart line" record  
C**CSMART(1) C   Tells if Hybrid program   
C  CSMART(2) C   Tells if RSD program 
C  CSMART(3) C   Reserv  
C  DTNAME   C    Holds the description of TECH1&2 
C  ICTECH   I    Used to check the TECH1&2 validity 
C  IMSKIP   I    Counts how many lines to skip in the I/M file 
C 
C  Notes: 
C 
C  18 Mar 98 @DynTel-MLA 2-659  This routine is new for Mobile6. This 
C            code was located in GETSC1, which is obsolete in Mobile6. 
C 
      IMPLICIT NONE 
      CHARACTER*120  BUFFER 
      CHARACTER*1    CSMART(3) 
      CHARACTER*80   DTNAME 
C 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'IM12HC.I' 
      INCLUDE 'IM12CO.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'RSDINF.I' 
C 
      INTEGER, EXTERNAL :: JUNIT
C
      INTEGER  RC, ITALD, ICTECH, INERR, IAGE1, ISTR, ITECH, IMSKIP, 
     *ISMART, IBY, IOS 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
      RC = 0 
C 
C  Do not read any I/M credits unless IMFLAG is greater than 1. Read new 
C  credits whenever alternate credits are requested. Also read new credits  
C  whenever region (high or low altitude) changes, and when standard  
C  credits are requested for the current scenario but alternate credits  
C  were used for the last scenario. The only time we will skip reading new  
C  credits is when standard credits were used for the last scenario and  
C  standard credits (for the same region) are requested for the current  
C  scenario. 
C 
      IF(IMFLAG.EQ.1) GOTO 98 
      IF (IFILET.EQ.1.AND.SVALTC.EQ.1.AND.IREAD.EQ.IREJN) GOTO 98 
C 
      IREAD = IREJN 
      SVALTC = IFILET 
C 
C     If the user supplied an alternate credit file name, it will be in 
C     the variable TNAME. Otherwise, we can use TNAME to hold the standard 
C     credit file name 'TECH12.D.' 
C 
      IF(IFILET.EQ.1) TNAME='TECH12.D' 
ccs   All M6 input files live in SMK_M6PATH
      CALL GETENV( "SMK_M6PATH", M6DIR )

      M6FILE = TRIM( M6DIR ) // '/' // TNAME
      
      IOUIMD = JUNIT()            
      OPEN(IOUIMD,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     * IOSTAT=IOS,ERR=96) 
C 
   20 READ(IOUIMD,141,END=95,ERR=96,IOSTAT=IOS) BUFFER 
      READ(UNIT=BUFFER,FMT=140,END=95,ERR=96,IOSTAT=IOS) ITALD, 
     *  ICTECH,(CSMART(ISMART), ISMART=1,3), DTNAME 
C      
C Check MOBILE5b or MOBILE5a format 
C 
      IF ((CSMART(1) .EQ. '1') .OR. (CSMART(1) .EQ. '2')) THEN  
C 
C     Check the tech type and altitude and if they are not right then 
C     jump back and get the next one after skipping 400 lines. 
C 
C 
C  26-Feb-1997 @DynTel-yc Req.2-637 The check for hybrid credits in TECH12.D 
C               was disabled to avoid error. 
C 
         IF( (ITALD.NE.IREJN) .OR.(ICTECH.NE.1) ) THEN 
            GOTO 21 
         ELSE    
            GOTO 30 
         ENDIF 
      ELSE      
C MOBILE5a format 
         READ(UNIT=BUFFER,FMT=142,END=95,ERR=96) ITALD,ICTECH,DTNAME 
         IF (ITALD.EQ.IREJN .AND. ICTECH.EQ.1) GOTO 30 
      ENDIF   
C      
   21 DO 25 IMSKIP=1,400 
        READ(IOUIMD,150) 
   25 CONTINUE 
  140 FORMAT(I1,1X,I1,T21,3A1,1X,A80) 
  141 FORMAT(120A) 
  142 FORMAT(I1,1X,I1,17X,A80)  
  150 FORMAT(1X) 
      GOTO 20 
C 
C     Reading the data in this manner accomplishes the same task as 
C     Block Data Subprograms 06, and 07 once did. 
C 
C     Loop for the two (2) ITECHs to be stored in the record. 
C 
   30 WRITE(IOUOUT,240) DTNAME 
  240 FORMAT(' Reading TECH1&2 I/M Credit Data File...'/A80) 
      WRITE(IOUREP,245) TNAME 
  245 FORMAT('*** I/M credits for Tech1&2 vehicles were read from', 
     &       ' the following external'/,'    data file: ',A) 
      SMRTL1=DTNAME 
      DO 45 ITECH=1,2 
C 
C       Loop for the five (5) stringency values in the record. 
C 
        DO 40 ISTR=1,5 
C 
C         Loop for the values in the record 
C 
          DO 35 IAGE1=1,20 
             READ(IOUIMD,160,END=95,ERR=96) 
     *           (CR12HC(IBY,IAGE1,ISTR,ITECH),IBY=1,19) 
   35     CONTINUE 
C 
   40   CONTINUE 
   45 CONTINUE 
C 
C     Loop for the two (2) ITECHs to be stored in the record. 
C 
      DO 60 ITECH=1,2 
C 
C       Loop for the five (5) stringency values in the record. 
C 
        DO 55 ISTR=1,5 
C 
C         Loop for the values in the record 
C 
          DO 50 IAGE1=1,20 
            READ(IOUIMD,160,END=95,ERR=96) 
     *          (CR12CO(IBY,IAGE1,ISTR,ITECH),IBY=1,19) 
  160       FORMAT(19F4.3) 
   50     CONTINUE 
C 
   55   CONTINUE 
   60 CONTINUE 
C 
      CLOSE(IOUIMD) 
      GOTO 98 
C 
   95 CALL QUITER(0.,0,54,INERR) 
      CLOSE(IOUIMD) 
      RC = IOS 
      GOTO 98 
   96 CALL QUITER(0.,0,55,INERR) 
      CLOSE(IOUIMD) 
      RC = IOS 
      GOTO 98 
C 
   98 RETURN 
      END 
