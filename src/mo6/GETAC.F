      SUBROUTINE GETAC(INPREC,NSCEN,INERR,RC) 
C 
C  GETALH reads in and validates the optional USER supplied parameters 
C  for the AC correction factor. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER and UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  21 Sep 00 AIR Task 07: Bug Fix 192: Added Quiter messages for user-supplied AC input. 
C            Bug Fix 256: Added storage of CLD_CVR into GLCD_CVR. Bug Fix 260: Applied 
C            storage of global values by case, rather than all at once. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  13 Jan 00 @DynTel-ddj 1-013, Added code to handle the CLOUD COVER, 
C            SUNRISE/SUNSET, PEAK SUN, NO AC, ABSOLUTE HUMIDITY and 
C            AC FRACTIN labels. 
C  06 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    argument list: INERR,INPREC 
C    common blocks: 
C    /FLAGS2/ ALHFLG 
C    /IOUCOM/ IOUGEN,IOUOUT 
C    /MAXIMA/ MAXLTW 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C    common blocks: 
C    /ALUHIN/ AC,ABSHUM,GLAHUM,GLAC,GSUNRISE,GSUNSET,GPEAKSUN,GCLD_CVR 
C 
C Local variables dictionary 
C 
C   Name        Type               Description 
C  ------      ------      -------------------------------------- 
C 
C  Notes: 
C 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. 
C Oct-27-97 DynTel@gjr 2-645 Removed code for interactive input. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER, INTENT(INOUT)    :: INERR 
      INTEGER, INTENT(IN)       :: NSCEN 
      INTEGER, INTENT(INOUT)    :: RC 
C 
      CHARACTER, INTENT(IN)     :: INPREC*(*) 
C 
      REAL DUMMY 
C 
      INTEGER IDUM(2) 
C 
      SELECT CASE (RC) 
C 
        CASE(1)  ! CLOUD COVER label 
          READ(INPREC(21:),*,IOSTAT=RC,END = 96,ERR=95) DUMMY 
          IF (DUMMY .LT. 0.0 .OR. DUMMY  .GT. 1.0) THEN 
             CALL QUITER(DUMMY, 7, 45, INERR) 
             CLD_CVR = 0.0 
          ELSE 
             CLD_CVR = DUMMY 
             CALL QUITER(CLD_CVR, 0, 617, INERR) 
             IF(NSCEN.EQ.0) GCLD_CVR = CLD_CVR 
          ENDIF 
C 
        CASE(2)  ! SUNRISE/SUNSET label 
          READ(INPREC(21:),*,IOSTAT=RC,END = 96,ERR=95) IDUM 
          IF (IDUM(1) .LT. 5 .OR. IDUM(1) .GT. 9) THEN 
             CALL QUITER(7.0, IDUM(1), 56, INERR) 
          ELSE IF(IDUM(2) .LT. 5 .OR. IDUM(2) .GT. 9) THEN 
               CALL QUITER(7.0, IDUM(2), 56, INERR) 
          ELSE 
C 
             CALL QUITER(FLOAT(IDUM(1)), IDUM(2), 618, INERR)  
C 
C         Set sunrise and sunset index for MOBILE6.  Index 1 is 
C         equivalent to 6 AM and index 24 is equivalent to 5 AM. 
C 
             SUNRISE = IDUM(1) - 5 
             IF(SUNRISE .LT. 1) SUNRISE = SUNRISE + 24 
             SUNSET  = IDUM(2) + 7 
             IF(NSCEN.EQ.0) THEN 
               GSUNRISE=SUNRISE 
               GSUNSET =SUNSET 
             ENDIF 
          ENDIF 
C 
C 
        CASE(3)  ! PEAK SUN label 
          READ(INPREC(21:),*,IOSTAT=RC,END = 96,ERR=95) IDUM 
          IF (IDUM(1) .LT. 9 .AND. IDUM(1) .GT. 1) THEN 
             CALL QUITER(7.0, IDUM(1), 57, INERR) 
          ELSE IF(IDUM(1) .GT. 12 ) THEN 
             CALL QUITER(7.0, IDUM(1), 57, INERR) 
          ELSE IF(IDUM(2) .LT. 1 .OR. IDUM(2) .GT. 5) THEN 
             CALL QUITER(7.0, IDUM(2), 57, INERR) 
          ELSE 
C 
             CALL QUITER(FLOAT(IDUM(1)), IDUM(2), 619, INERR) 
C 
C         Set PEAKSUN indices for MOBILE6.  Index 1 is 
C         equivalent to 6 AM and index 24 is equivalent to 5 AM. 
C 
             PEAKSUN(1) = IDUM(1) - 5 
             IF(PEAKSUN(1) .LT. 1) PEAKSUN(1) = PEAKSUN(1) + 12 
             PEAKSUN(2) = IDUM(2) + 7 
             IF(NSCEN.EQ.0) GPEAKSUN = PEAKSUN 
          ENDIF 
C 
C 
        CASE(4)  ! NO AC label 
          AC_FLAG = 0 
          RC = 0 
          CALL QUITER(0., 0, 620, INERR) 
C 
        CASE(5)  ! ABSOLUTE HUMIDITY label 
          READ(INPREC(21:),*,IOSTAT=RC,END = 96,ERR=95) DUMMY 
          IF(DUMMY .LT. 20.0 .OR. DUMMY .GT. 528.0) THEN 
             CALL QUITER(DUMMY, 0 , 47 , INERR) 
          ELSE 
             ABSHUM = DUMMY 
             IF(NSCEN.EQ.0) GLAHUM = ABSHUM 
          END IF 
C 
        CASE(6)  ! AC FRACTION label 
          READ(INPREC(21:),*,IOSTAT=RC,END = 96,ERR=95) DUMMY 
          IF(DUMMY .LT. 0.0 .OR. DUMMY .GT. 1.0) THEN 
             CALL QUITER(DUMMY,0,44,INERR) 
          ELSE 
             AC = DUMMY 
             IF(AC .GT. 0.0) CALL QUITER(AC, 0, 44, INERR) 
             IF(AC .EQ. 0.0) AC_FLAG = 0 
             IF(NSCEN.EQ.0) GLAC = AC 
          END IF 
C 
        CASE DEFAULT 
C 
      END SELECT 
C 
      RETURN 
C 
C  Error or EOF on any attempted read. 
C 
   95 RC = 1 
      GOTO 99 
C 
   96 RC = -1 
C 
   99 CALL QUITER(0.0,0,133,INERR) 
C 
      RETURN 
      END 
