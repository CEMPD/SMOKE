      SUBROUTINE GETIPR(INPREC,INERR,RC) 
C 
C  GETIPR processes the I/M PROGRAM labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls UCCOMP,NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C 25 May 01 @EPA-ELG Added loaded idle 
C 11 Oct 00 AIR Task 02: Bug Fix 275: Improved read error handling in GETM5IMC; put 
C           INERR back in parameter list. 
C 29 Sep 00 AIR Task 07 Fixed parameters in call to GETM5IMC; remove INERR 
C 27 Sep 00 AIR Task 07: Bug Fix 255: Free format reads 
C 25 Aug 00 AIR Task 08:       Added processing for EVAP I/M 
C 17 Jul 00 @EPA-elg           Routine now calls GETM5IM for ASM I/M programs. 
C 16 Mar 99 @DynTel-RJD 2-694  Added OBD I/M to the list of test types. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /CHKLST/ IMCHECK 
C    /FLAGS2/ IMFLAG 
C    /IMPAR1/ ICYIM 
C    /IMPAR2/ ITEST 
C    /IMPAR6/ IFREQ,INTYP 
C    /LDGIMC/ IMSTRICT 
C    /TTCINF/ IASM 
C    /YEARS4/ IYSTART,IYEND,IY1960 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  IMPGM      I   The I/M program number. 
C  NDATA      I   Loop control variable. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  TESTFIELD  C   Upper case version of the program type input field. 
C  X          R   Temporary real to integer variable 
C 
C  Notes: 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'TTCINF.I' 
      INCLUDE 'YEARS4.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        IMPGM 
      INTEGER        IOS 
      INTEGER        NDATA 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
      CHARACTER*30   TESTFIELD 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C   IMFLAG = 2  IMFLAG is now set in TESTIM. 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program and check off 
C     the other data items that this subroutine is supposed to read. 
C 
      DO NDATA = 1,6 
        IMCHECK(NDATA,IMPGM) = 2 
      END DO 
C 
C     Look for the program start and end dates. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      ICYIM(1,IMPGM)=INT(X) 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      ICYIM(2,IMPGM)=INT(X) 
      IF (ICYIM(1,IMPGM).LT.IY1960.OR.ICYIM(1,IMPGM).GT.IYEND) THEN 
        CALL QUITER(0.0,ICYIM(1,IMPGM),35,INERR) 
        RC = RC + 2 
      ELSEIF (ICYIM(2,IMPGM).LT.IY1960.OR. 
     &        ICYIM(2,IMPGM).GT.IYEND) THEN 
        CALL QUITER(0.0,ICYIM(2,IMPGM),36,INERR) 
        RC = RC + 3 
      ELSEIF (ICYIM(1,IMPGM).GT.ICYIM(2,IMPGM)) THEN 
        CALL QUITER(0.0,0,37,INERR) 
        RC = RC + 4 
      ENDIF 
C 
C     Look for the program frequency. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) X 
      IFREQ(IMPGM)=INT(X) 
      IF (IFREQ(IMPGM).LT.1.OR.IFREQ(IMPGM).GT.3) THEN 
        CALL QUITER(0.,IFREQ(IMPGM),101,INERR) 
        RC = RC + 5 
      ENDIF 
C 
C     Look for the program type (test-and-repair or test-only). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
      IF (TESTFIELD.EQ.'T/O') THEN 
        INTYP(IMPGM) = 1 
      ELSEIF (TESTFIELD.EQ.'TRC') THEN 
        INTYP(IMPGM) = 2 
      ELSEIF (TESTFIELD.EQ.'TRM') THEN 
        INTYP(IMPGM) = 3 
      ELSE 
        CALL QUITER(0.0,0,82,INERR) 
        RC = RC + 6 
      ENDIF 
C 
C     Finally, look for the test type (Idle, ASM, IM240, etc.). 
C     Look at all of the characters from the end of the last data field 
C     up to the end of the record. (Note: IASM is used only by OUTIMP.) 
C     Set IMSTRICT based on the program type. Assume lenient cutpoints 
C     for IDLE and ASM PHASE-IN programs and strict cutpoints for ASM 
C     FINAL programs. For IM240 programs, the user must enter cutpoints 
C     which may vary by model year. Cutpoints will be judged as strict 
C     or lenient at the appropriate point in the code for each model 
C     year by examining the NOX cutpoint entered by the user for that 
C     model year. 
C 
      IASM(IMPGM) = 1 
      IMSTRICT(IMPGM) = 0 
      CALL UCCOMP(INPREC(PTR2:),TESTFIELD) 
      IF (TESTFIELD.EQ.'IDLE') THEN 
        ITEST(IMPGM) = 1 
        IMSTRICT(IMPGM) = 1 
      ELSE IF (TESTFIELD.EQ.'2500/IDLE') THEN 
        ITEST(IMPGM) = 2 
        IMSTRICT(IMPGM) = 1 
      ELSEIF (TESTFIELD.EQ.'LOADED/IDLE') THEN 
        ITEST(IMPGM) = 2 
        IMSTRICT(IMPGM) = 1 
      ELSE IF (TESTFIELD.EQ.'I/M 240'.OR. 
     &         TESTFIELD.EQ.'IM240') THEN 
        ITEST(IMPGM) = 4 
        IMSTRICT(IMPGM) = 0 
C     Note: ASM TEST has disappeared... 
CC    ELSE IF (TESTFIELD.EQ.'RSD ONLY') THEN 
CC      ITEST(IMPGM) = 5 
      ELSE IF (TESTFIELD.EQ.'ASM 2525 PHASE-IN') THEN 
        ITEST(IMPGM) = 6 
        IMSTRICT(IMPGM) = 1 
      ELSE IF (TESTFIELD.EQ.'ASM 2525 FINAL') THEN 
        ITEST(IMPGM) = 7 
        IMSTRICT(IMPGM) = 2 
      ELSE IF (TESTFIELD.EQ.'ASM 5015 PHASE-IN') THEN 
        ITEST(IMPGM) = 8 
        IMSTRICT(IMPGM) = 1 
      ELSE IF (TESTFIELD.EQ.'ASM 5015 FINAL') THEN 
        ITEST(IMPGM) = 9 
        IMSTRICT(IMPGM) = 2 
      ELSE IF (TESTFIELD.EQ.'ASM 2525/5015 PHASE-IN') THEN 
        ITEST(IMPGM) = 10 
        IMSTRICT(IMPGM) = 1 
      ELSE IF (TESTFIELD.EQ.'ASM 2525/5015 FINAL') THEN 
        ITEST(IMPGM) = 11 
        IMSTRICT(IMPGM) = 2 
      ELSE IF (TESTFIELD.EQ.'OBD I/M') THEN 
        ITEST(IMPGM) = 12 
        IMSTRICT(IMPGM) = 1 
C 
C     Look for the four EVAP I/M test types: EVAP OBD,  
C     EVAP OBD & GC, FP & GC, and GC. Note: starting 
C     the evap-type I/M at 20 allows for expansion of 
C     exhaust-type programs. 
C 
      ELSE IF (TESTFIELD.EQ.'EVAP OBD') THEN 
        ITEST(IMPGM) = 20 
      ELSE IF (TESTFIELD.EQ.'EVAP OBD & GC') THEN 
        ITEST(IMPGM) = 21 
      ELSE IF (TESTFIELD.EQ.'FP & GC') THEN 
        ITEST(IMPGM) = 22 
      ELSE IF (TESTFIELD.EQ.'GC') THEN 
        ITEST(IMPGM) = 23 
C 
      ELSE 
        CALL QUITER(0.,INTYP(IMPGM),162,INERR) 
        RC = RC + 7 
      ENDIF 
C 
      IF (ITEST(IMPGM).GE.6.AND.ITEST(IMPGM).LE.11) THEN 
        IASM(IMPGM) = 2 
        CALL GETM5IM(ITEST(IMPGM),RC,INERR) 
      ENDIF 
C 
      IF (RC.NE.0) THEN 
        INERR = INERR + 1 
      ENDIF 
C 
      GOTO 99 
C 
C Error or EOF on any attempted read. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
C       
   99 RETURN 
      END 
