C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETIM2(INERR,*)
C
C  GETIM2 reads in and processes the parameter record for I/M program #2.
C
C  Calls QUITER and YRTEST.
C
C  Input on call:
C
C     common blocks:
C     /FLAGS1/ PROMPT
C     /IOUCOM/ IOUGEN,IOUASK
C     /YEARS4/ IYSTART,IY1960,IYEND
C
C  Output on return:
C
C     parameter list: INERR
C     common blocks:
C     /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C     /IMPAR2/ ILDT
C     /IMPAR4/ CUTHC,CUTCO,CUTNO
C     /IMPAR6/ IFREQ,INTYP,IMEFFR
C     /TTCINF/ IHYB,IASM
C
C  Local variable / array dictionary:
C
C   Name   Type                      Description
C  ------  ----  ---------------------------------------------------------------
C  DCTFLG   C    Flag - user will enter effectiveness ratios for I/M pgms.
C  IVIM     I    See top of program
C  JCUTFL   I    Analagous to the ICUTFL from GETIMP
C  ECUTHC,ECUTCO,ECUTNO
C           R    Temporary variables to hold cutpoints
C  JFILET   I    Temporary variables to hold alternate I/M credit flags
C                from the second I/M program parameter record. 
C
C
C  Notes:
C
C  GETIM2 was added to MOBILE5.  It allows for a second I/M program.
C  /IMPAR4/ was included here to allow GETIM2 to gather non-default
C  cutpoints
C
C  GETIM2 was modified for request 2-364 to permit the specification of the
C  alternate I/M credit flags on either the first or the second I/M 
C  program parameter record. TMM@ARC 12/22/93
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Modified the "Program type" and "Test type" checking.
C  Set values of IHYB and IASM.
C 
C  8 Feb 96 @DynTel-MLA 2-611 Added DCTFLG to indicate whether user will
C     enter alternate discount rates for test-and-repair I/M programs.
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'YEARS4.I'
      INCLUDE 'TTCINF.I'
C
C     Define EFFFLG and define and initialize temporary I/M credit flags.
C   
      INTEGER    JFILET(2)
      CHARACTER  DCTFLG
C
C
      JFILET(1) = 1
      JFILET(2) = 1 
C
C     Initialize JCUTFL to indicate use default cutpoints
C
      JCUTFL=1
      IF(PROMPT.EQ.2) THEN
         WRITE(IOUASK,200)
  200    FORMAT(' ','Enter the I/M Program #2 parameters record:'/
     *          '0','Program start year, pre-81 stringency,'/
     *          ' ','first and last model year getting benefits,'/
     *          ' ','old tech waiver rate, new tech waiver rate,'/
     *          ' ','compliance rate,'/
     *          ' ','inspection type,'/
     *          ' ','flag for effectiveness rates (optional),'/
     *          ' ','frequency of I/M Program #2 inspection,'/
     *          ' ','the format is:'/
     *          ' ','(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1)')
C
         READ(IOUGEN,100,ERR=99,END=99)
     *   ICYIM(2),ISTRIN(2),MODYR(1,2),MODYR(2,2),
     *   WAIVER(1,2),WAIVER(2,2),CRIM(2),INTYP(2),DCTFLG,IFREQ(2)
  100    FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1)
C
         WRITE(IOUASK,210)
  210    FORMAT(' ','Enter the rest of the I/M Program #2 parameters'
     *             ,' record:'/
     *          '0','Vehicle classes covered, test type,'/
     *          ' ','cutpoint flag (optional),'/
     *          ' ','and any alternate cutpoints (contingent)'/
     *          ' ','the format is:'/
     *          ' ','4I1,1X,I1,I1 (optional),3X,2(F4.2,1X),F4.2'
     *             ,' (contingent)')
C
         READ(IOUGEN,110,ERR=99,END=99)
     *   (ILDT(IVIM,2),IVIM=1,4),ITEST(2),
     *   JCUTFL,ECUTHC,ECUTCO,ECUTNO
  110    FORMAT(4I1,1X,2I1,3(1X,F4.2))
         GOTO 10
      ENDIF
C
C MOBILE5 uses a different format for interactive input because
C the format description for interactive input is too long to
C on one line.
C
      READ(IOUGEN,120,ERR=99,END=99)
     * ICYIM(2),ISTRIN(2),MODYR(1,2),MODYR(2,2),WAIVER(1,2),WAIVER(2,2),
     * CRIM(2),INTYP(2),DCTFLG,IFREQ(2),(ILDT(IVIM,2),IVIM=1,4),
     * ITEST(2),JCUTFL,JFILET(1),JFILET(2),ECUTHC,ECUTCO,ECUTNO
  120 FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1,1X,4I1,1X,4I1
     *     ,3(1X,F4.2))
C
C     If either JFILET(1) or JFILET(2) is 2, set the corresponding
C       IFILET to 2, otherwise leave the IFILET values alone.
C     TMM@ARC 12/22/93 2-364
C
      IF (JFILET(1).EQ.2) IFILET(1)=2
          IF (JFILET(2).EQ.2) IFILET(2)=2      
C
C
C  Check I/M Program #2 inputs.
C
   10 CALL YRTEST(ICYIM(2),32,IY1960,IYEND,INERR)
      IF(ISTRIN(2).LT.10.OR.ISTRIN(2).GT.50)
     *   CALL QUITER(0.,ISTRIN(2),33,INERR)
C
C little is known what in fact these subroutines do
C
      CALL YRTEST(MODYR(1,2),35,IYSTART,IYEND,INERR)
      CALL YRTEST(MODYR(2,2),36,IYSTART,IYEND,INERR)
      IF(MODYR(1,2).GT.MODYR(2,2)) CALL QUITER(0.,0,37,INERR)
      IF(WAIVER(1,2).LT.0.0.OR.WAIVER(1,2).GT.50)
     *   CALL QUITER(WAIVER(1,2),0,100,INERR)
      IF(WAIVER(2,2).LT.0.0.OR.WAIVER(2,2).GT.50)
     *   CALL QUITER(WAIVER(2,2),0,100,INERR)
      IF(CRIM(2).LT.0.0.OR.CRIM(2).GT.100.)
     *   CALL QUITER(CRIM(2),0,34,INERR)
      IF(INTYP(2).LT.1.OR.INTYP(2).GT.4)
     *   CALL QUITER(0.,INTYP(2),82,INERR)
C     Don't reset DCTFLG to 1. It may have been set in GETIMP!
      IF(DCTFLG.EQ.'2')THEN
        IMEFFR = 2
      ELSEIF(DCTFLG.NE.' '.AND.DCTFLG.NE.'1') THEN
        CALL QUITER(0.0, 0, 102, INERR)
      ENDIF
      IF(IFREQ(2).LT.1.OR.IFREQ(2).GT.2)
     *   CALL QUITER(0.,IFREQ(2),101,INERR)
C
C  Convert the waiver rate from percent to percentage for PCLEFT.
C
      WAIVER(1,2)=WAIVER(1,2)/100.
      WAIVER(2,2)=WAIVER(2,2)/100.
C
      DO 40 IVIM=1,4
        IF(ILDT(IVIM,2).LT.1.OR.ILDT(IVIM,2).GT.2)
     *    CALL QUITER(0.,ILDT(IVIM,2),72,INERR)
   40 CONTINUE
C     
      IF(ITEST(2).LT.1.OR.ITEST(2).GT.5)
C
     *    CALL QUITER(0.,ITEST(2),73,INERR)
      IF(JCUTFL.EQ.2) THEN
         CUTHC(2)=ECUTHC
         CUTCO(2)=ECUTCO
         CUTNO(2)=ECUTNO
      ENDIF
C
C  User specifies a Retest Based Hybrid Program Design.
C
        IHYB(2)=1
        IF(INTYP(2).EQ.4) THEN
          INTYP(2)=1
          IHYB(2)=2
        ENDIF
C
C  User specifies the ASM test procedure.
C
        IASM(2)=1
        IF(ITEST(2).EQ.5) THEN
          ITEST(2)=3
          IASM(2)=2
        ENDIF
C      
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 CALL QUITER(0.,0,129,INERR)
C
      RETURN1
      END
