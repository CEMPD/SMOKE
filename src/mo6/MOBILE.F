      SUBROUTINE MOBILE(MINPUT,NRUN,NFILE,INERR,RC) 
C 
C  Called by user's program (example: Air Quality Analysis 
C  System). 
C 
C  Calls BDSAVE, CHKINP, CHKLAB, DSPTCH, FNDLAB, INIFLG, 
C        OUTHDR, OUTNEW, PRCSCN, REINIT, and RESTOR. 
C 
C  Changes: (Last change first) 
C 
C  21 Feb 01 AIR Task 11: Included code to print header for PM runs. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  24 Aug 99 @DynTel-MLA 2-684  Renamed some variables. Avoid using MIN 
C            (the name of an intrinsic function) as a variable name, etc. 
C            Removed call to INILAB. (It is now called in DRIVER.) 
C  28 Jul 99 @DynTel-MLA 2-694  Added call to DALSCN. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Added call to INILAB, rewrote calls to 
C            CHKLAB because the calling sequence was changed. 
C  24 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  The following variables are passed in argument lists or common 
C  blocks and in MOBILE: 
C 
C  Input on call: 
C 
C    argument lists: MINPUT, NRUN, NFILE 
C 
C    common blocks: 
C    /EVAL/   MEVAL,MEVOLD 
C    /FLAGS2/ NEWFLG 
C    /IOUCOM/ IOUERR 
C    /REGION/ IREJN,IROLD 
C    /SAVE01/ JCALL 
C    /YEARS4/ IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  ENDRUN   I    Token number for the "END OF RUN" input label. 
C  ICY      I    Calendar year for the current scenario. 
C  ICYOLD   I    Calendar year from the previous scenario. 
C  INPREC   C    The next labeled input record from the input file. 
C  LABEL    C    The label field from the current labeled input record. 
C  LABTYPE  I    Token number for the label on the current input rec. 
C  NSCEN    I    Number of scenarios in the current test run. 
C  SCENARIO I    Token number for the "SCENARIO RECORD" input label. 
C  SECTION  I    Input section where a record belongs (1=global only, 
C                2=global or scenario, 3=scenario only) 
C  STDLAB   C    The full standard text for an input label. 
C 
C  Notes: 
C 
C  For MOBILE4.1 the subroutine MOBILE was created from MAIN. 
C  The /BYMYCd/ CBs were added. 
C  May-18-1993 @ ARC-bsg Subtask 238 the common block IMPAR8 
C  was removed, it contained the variable ICALC 
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5 
C  was modified, the HDGV I/M credit was expanded to include 
C  an initial value for NOx. 
C  24-August-1994 @CSC-pme request 2-446 Include file BASE12.I was removed 
C  (24-May-1996) @DynTel-yc Request 2-621 The calendar year was expanded 
C                 from year 2020 to year 2051 
C  (7-April-1997)@Dyntel-rjd (after gjr) Request 2-204 Changed calculations  
C                 for annual mileage and average accumulation rates when  
C                 calculating July emissions.  July emissions will be  
C                 calculated by aging the VMT mix by 6 months rather that  
C                 interpolating July emissions from two consecutive  
C                 Januaries.  Subroutines adjust and saver and function  
C                 zcalc no longer needed. Interpolation loop removed. 
C   16 Mar 98 @DynTel-MLA 2-659  This subroutine was completely rewritten 
C             for Mobile6. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'SAVE01.I' 
C 
C     Parameter declarations. 
C 
      INTEGER     INERR 
      CHARACTER   MINPUT*(*) 
      INTEGER     NFILE 
      INTEGER     NRUN 
      INTEGER     RC 
C 
C     Local variable declarations. 
C 
      INTEGER     ENDRUN 
      INTEGER     ICY 
      INTEGER     ICYOLD 
      CHARACTER   INPREC*150 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      INTEGER     NSCEN 
      INTEGER     SCENARIO 
      INTEGER     SECTION 
      CHARACTER   STDLAB*40 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Subroutine call, setup run. 
C 
      CALL REINIT 
C 
      IF(JCALL.EQ.1) THEN  
        CALL BDSAVE 
      ELSE  
        CALL RESTOR 
      ENDIF 
      JCALL=2 
C 
C     Initialize non-COMMON and local variables. 
C 
      ICYOLD = -1 
      INERR = 0 
      NSCEN = 0 
C 
C     Find the token numbers for the "SCENARIO RECORD" and the  
C     "END OF RUN" record labels. 
C 
      LABEL = 'SCENARIO RECORD' 
      CALL CHKLAB(LABEL,STDLAB,SCENARIO,SECTION) 
      LABEL = 'END OF RUN' 
      CALL CHKLAB(LABEL,STDLAB,ENDRUN,SECTION) 
C 
C     Set the internal flags (Control Section flags  
C     in Mobile5) to default values. 
C 
      CALL INIFLG 
C 
C     Check for end-of-file on the input file. Ignore blank lines at EOF. 
C 
      INPREC = ' ' 
   15 READ(IOUGEN,110,IOSTAT=RC,ERR=9900,END=9900) INPREC 
  110 FORMAT(A) 
      IF (INPREC.EQ.' ') GOTO 15 
      BACKSPACE(IOUGEN) 
C 
C     Write the run header to the output files and a progress message  
C     to the console. 
C 
      CALL OUTHDR(IOUREP,MINPUT,NFILE,NRUN) 
C 
      IF(PARTFLG) CALL OUTHDR(IOUPM,MINPUT,NFILE,NRUN) 
C 
      IF(IOUOUT.EQ.6.OR.IOUOUT.EQ.9) WRITE(IOUOUT,120) 
  120 FORMAT(' Reading information.') 
C 
C     Pick up the next labeled input record from the input file. 
C 
   20 CALL FNDLAB(INPREC,STDLAB,LABTYPE,SECTION,INERR,RC) 
C 
C     Special processing may be required for Scenario and End of Run recs. 
C 
      IF (LABTYPE.EQ.SCENARIO.OR.LABTYPE.EQ.ENDRUN) THEN 
C 
C       Check the global or One-Time data (any data entered before  
C       the first scenario record) for consistency and completeness. 
C 
        IF (LABTYPE.EQ.SCENARIO.AND.NSCEN.EQ.0) THEN 
          CALL CHKINP(INERR,RC) 
        ENDIF 
C 
C       Finish processing the previous scenario, if any. 
C 
        IF (NSCEN.GT.0) THEN 
C 
          CALL PRCSCN(ICY,ICYOLD,INERR,RC) 
C 
C         Return any storage dynamically allocated for the last scenario. 
C 
          CALL DALSCN() 
C 
        ENDIF 
C 
      ENDIF 
C 
C     Process the current input record. 
C 
      CALL DSPTCH(INPREC,STDLAB,SECTION,NFILE,NSCEN,NRUN,ICY, 
     &            INERR,RC) 
C 
C     Get the next input record. 
C 
      IF (LABTYPE.NE.ENDRUN) THEN 
        GOTO 20 
      ENDIF 
C 
C     Echo the alternate BER rates entered by the user even if some 
C     error prevented processing. 
C 
      IF ((NSCEN.EQ.0).AND. 
     *   (NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)) THEN 
        CALL OUTNEW() 
      ENDIF 
C 
 9900 RETURN 
      END 
