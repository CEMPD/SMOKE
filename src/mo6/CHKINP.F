      SUBROUTINE CHKINP(INERR,RC) 
C 
C  CHKINP performs some checks on the input data before the first scenario   
C  is processed. 
C 
C  Called by MOBILE 
C 
C  Calls CALCUM, CALEFR, CEVTPD, CHKMRN, GETCUM, INITIM, SAVEP, 
C        and TESTIM. 
C 
C  Changes: (Last change first) 
C 
C  20 Jun 01 @EPA-djb Bug326, Removed call to CHKPGM (I/M checks). 
C  25 Sep 00 @EPA-djb Bug238, Removed call to DBWRUN (database output). 
C  07 Jul 00 Removed unsed LPOD from CALCUM call. 
C  24 Aug 99 @DynTel-MLA 2-684 Call database output routines. 
C  18 May 99 @DynTel-HXQ 2-650  Delete call to CEVTPD 
C  11 Feb 99 @DynTel-MLA 2-663  Added call to INITIM 
C  19 Nov 98 @DynTel-MLA 2-663  Removed call to GETIMC. That routine is 
C            obsolete in Mobile6. Removed call to GETTTC. GETTTC has 
C            been modified and is now called directly from DSPTCH. 
C  12 Nov 98 @DynTel-MLA 2-663  Added call to TESTIM and replaced a call 
C            to GETEFR with a call to CALEFR. 
C  26 Aug 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    /FLAGS2/ IMFLAG,NEWFLG 
C    /FLAGS4/ COLDFG 
C    /IMPAR6/ IMEFFR 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
C  This subroutine is new for Mobile6. 
      USE DATABASE, ONLY: DBFLAG 
      IMPLICIT NONE 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IMPAR6.I' 
 
      INTEGER INERR 
      INTEGER    RC 
C 
      RC = 0 
C 
C  Check for consistency between the mileage accumulation rates and the 
C  registration distributions. JULMYR is normalized even if no 
C  inconsistencies are found. 
C 
      CALL CHKMRN(INERR,RC) 
      IF (RC.NE.0.AND.INERR.EQ.0) INERR = 1 
C 
C     Make sure I/M program inputs are complete and consistent and 
C     initialize some I/M variables. 
C 
      CALL TESTIM(INERR,RC) 
      IF (INERR.EQ.0) THEN 
        CALL INITIM() 
      ENDIF 
C 
C     CHKPGM checks for consistency between the I/M, ATP, Pressure, and 
C     Purge check program paramters. 
C 
C     These checks have *not* been updated for MOBILE6 and have been 
C     dropped. 
C 
C     CALL CHKPGM(INERR) 
C 
C     If IMEFFR is 2, an I/M EFFECTIVENESS labeled input record was 
C     read. Some processing was deferred until the I/M, ATP, and 
C     purge/pressure program inputs were complete. Call CALEFR here 
C     to complete processing for user-supplied effectiveness rates. 
C 
      IF (IMEFFR.EQ.2) THEN 
        CALL CALEFR(INERR,RC) 
        IF (RC.NE.0) GOTO 90 
      ENDIF 
C 
C     TAMFLG, MYMRFG, IMFLAG, and ATPFLG control inputs affecting the lead 
C     phaseout date tampering feature. Set up the lead phaseout save area. 
C 
      CALL CALCUM() 
      CALL GETCUM 
      CALL SAVEPB 
C 
C 
      IF(NEWFLG.GE.5) THEN 
         COLDFG=1 
         CALL QUITER(0.,0,135,INERR) 
      ENDIF 
C 
   90 RETURN 
      END 
