C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETIMP(INERR,*)
C
C  GETIMP reads in and processes the parameter record for I/M Program #1
C
C  Called by ONESEC.
C
C  Calls QUITER and YRTEST.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IY1941,IY1960,IY2020
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT,ITEST
C    /IMPAR4/ CUTHC,CUTCO,CUTNO
C    /IMPAR6/ IFREQ,INTYP,IFILET
C
C  Local variable / array dictionary:
C
C   Name   Type                      Description
C  ------  ----  ---------------------------------------------------------------
C  IVIM     I    See top of program
C  ICUTFL   I    Indicates alternate cutpoints when equal to 2
C  DCUTHC,DCUTCO,DCUTNO
C           R    Temporary variables to hold cutpoints
C
C  Notes:
C
C  GETIMP was modified in MOBILE4.1 to exclude the IM240 program.
C  and again in MOBILE5 v2.0 to check for new user inputs in the I/M parameter
C  these new inputs allow the user to specify their own cutpoints, and filenames
C  or they can mix and match between their own values and default MOBILE5 values.
C  Note TEST TYPE is right next to new flags in MOBILE5, this is done to allow
C  the I/M parameter record to be more compatible with MOBILE4.1, recall that
C  there once were 2 (two) indicators of alternate I/M credits located at the
C  end of the I/M parameter record and they were in the following FORTRAN
C  format 1X, 2I1
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR4/ CUTHC(2),CUTCO(2),CUTNO(2)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C
      INTEGER UNIT
C
C     Variables for the my VS cut-points test
C
      DIMENSION CH(15),CHT(10)
C
      CHARACTER*1 CH,CHT
      CHARACTER*2 CM24YR,CPRGYR,CPRSYR
C
      DATA CM24YR /'  '/
      DATA CPRGYR /'  '/
      DATA CPRSYR /'  '/
      DATA CHT /'1','2','3','4','5','6','7','8','9','0'/
C
      DO 5 I=1 ,15
         CH(I)='$'
 5    CONTINUE
C
      UNIT=IOUGEN
C
C     Initialize ICUTFL to indicate use default cutpoints.
C
      ICUTFL=1
C
C     User Imputed Values
      IF(PROMPT.EQ.2) THEN
C
         WRITE(IOUASK,200)
  200    FORMAT(' ','Enter the I/M Program parameters record:'/
     *          '0','Program start year, stringency,'/
     *          ' ','first and last model year getting benefits,'/
     *          ' ','old tech waiver rate, new tech waiver rate,'/
     *          ' ','compliance rate,'/
     *          ' ','inspection type,'/
     *          ' ','frequency of I/M inspection,'/
     *          ' ','the format is:'/
     *          ' ','(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,1X,I1)')
C
         READ(IOUGEN,100,ERR=99,END=99)
     *   ICYIM(1),ISTRIN(1),MODYR(1,1),MODYR(2,1),
     *   WAIVER(1,1),WAIVER(2,1),CRIM(1),INTYP(1),IFREQ(1)
  100    FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,1X,I1)
C
         WRITE(IOUASK,210)
  210    FORMAT('0','Enter the rest of the I/M Program parameters '
     *             ,'record:'/
     *          '0','vehicle classes covered, test type,'/
     *          ' ','cutpoint flag (optional),'/
     *          ' ','alternate filename flags (optional)'/
     *          ' ','and any alternate cutpoints (contingent)'/
     *          ' ','the format is:'/
     *          ' ','4I1,1X,I1,I1(optional),I1(optional),'
     *             ,'I1(optional),3(1X,F4.2)(contingent)')
C
C
         READ(IOUGEN,110,ERR=99,END=99)
     *   (ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     *   ICUTFL,IFILET(1),IFILET(2),DCUTHC,DCUTCO,DCUTNO
  110    FORMAT(4I1,1X,4I1,3(1X,F4.2))
C
         GOTO 25
       ENDIF
C
C MOBILE5 uses a different format for interactive input because
C the format description forinteractive input is too long to
C on one line.
C
        READ(IOUGEN,120,ERR=99,END=99)
     * ICYIM(1),ISTRIN(1),MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),
     * CRIM(1),INTYP(1),IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     * ICUTFL,IFILET(1),IFILET(2)
  120 FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,1X,I1,1X,4I1,1X,4I1)
C
C
C     Try and find The *OLD* M4.1 YEARS
C
      BACKSPACE(UNIT)
      READ(IOUGEN,101) CM24YR,CPRGYR,CPRSYR
 101  FORMAT(T36,3(1X,A2))
C
C     Read the character array CH
C
      BACKSPACE(UNIT)
      READ(IOUGEN,102) (CH(J),J=1,15)
 102  FORMAT(T36,15(A1))
C
C     try and find a "." period for floating point number
C     if true go and reread the I/M input line
      DO 15 J=1,15
         IF (CH(J).EQ.'.') GOTO 20
 15   CONTINUE
C
C     check for The three yrs and /or 15 blanks in a row
C
      IF((CM24YR.EQ.'  '.AND.
     *      CPRGYR.EQ.'  '.AND.
     *      CPRSYR.EQ.'  ').OR.(
     *       CH(1).EQ.' '.AND.CH(2).EQ.' '.AND.CH(3).EQ.' '
     *.AND.CH(4).EQ.' '.AND.CH(5).EQ.' '.AND.CH(6).EQ.' '
     *.AND.CH(7).EQ.' '.AND.CH(8).EQ.' '.AND.CH(9).EQ.' '
     *.AND.CH(10).EQ.' '.AND.CH(11).EQ.' '.AND.CH(12).EQ.' '
     *.AND.CH(13).EQ.' '.AND.CH(14).EQ.' '.AND.CH(15).EQ.' '
     *)) GOTO 25
C
C     Check for a character string or numbers
C     if only characters *NO* M152 Warning...
C     1st M4.1 yr: CM24YR
C
      IF (CH(1).EQ.' '.AND.CM24YR.NE.'  ') THEN
           DO 60 I=1,10
           IF (CH(2).EQ.CHT(I)) THEN
              DO 61 J=1,10
              IF (CH(3).EQ.CHT(I)) THEN
                 CALL QUITER(0.,0.,152,INERR)
                 GOTO 25
              ENDIF
 61           CONTINUE
          ENDIF
 60       CONTINUE
      ENDIF
C
C     2nd M4.1 yr: CPRGYR
C
      IF (CH(4).EQ.' '.AND.CPRGYR.NE.'  ') THEN
           DO 64 I=1,10
           IF (CH(5).EQ.CHT(I)) THEN
              DO 65 J=1,10
              IF (CH(6).EQ.CHT(I)) THEN
                 CALL QUITER(0.,0.,152,INERR)
                 GOTO 25
              ENDIF
 65           CONTINUE
          ENDIF
 64       CONTINUE
      ENDIF
C
C     3rd M4.1 yr: CPRSYR
C
      IF (CH(7).EQ.' '.AND.CPRSYR.NE.'  ') THEN
           DO 68 I=1,10
           IF (CH(8).EQ.CHT(I)) THEN
              DO 69 J=1,10
              IF (CH(9).EQ.CHT(I)) THEN
                 CALL QUITER(0.,0.,152,INERR)
                 GOTO 25
              ENDIF
 69           CONTINUE
          ENDIF
 68       CONTINUE
      ENDIF
C
C     Skip the Read
C
      GOTO 25
C
C     Read in the M5.0 F4.2 Cut-Points
C
 20   CONTINUE
C
      BACKSPACE (UNIT)
      READ(IOUGEN,121,ERR=99,END=99) DCUTHC,DCUTCO,DCUTNO
 121  FORMAT(T36,3(1X,F4.2))
C
C
C  Check I/M inputs.
C
 25   CONTINUE
        CALL YRTEST(ICYIM(1),32,IY1960,IY2020,INERR)
        IF(ISTRIN(1).LT.10.OR.ISTRIN(1).GT.50)
     *   CALL QUITER(0.,ISTRIN(1),33,INERR)
        CALL YRTEST(MODYR(1,1),35,IY1941,IY2020,INERR)
        CALL YRTEST(MODYR(2,1),36,IY1941,IY2020,INERR)
        IF(MODYR(1,1).GT.MODYR(2,1)) CALL QUITER(0.,0,37,INERR)
        IF(WAIVER(1,1).LT.0.0.OR.WAIVER(1,1).GT.50.)
     *             CALL QUITER(WAIVER(1,1),0,100,INERR)
        IF(WAIVER(2,1).LT.0.0.OR.WAIVER(2,1).GT.50.)
     *             CALL QUITER(WAIVER(2,1),0,100,INERR)
        IF(CRIM(1).LT.0.0.OR.CRIM(1).GT.100.)
     *             CALL QUITER(CRIM(1),0,34,INERR)
        IF(INTYP(1).LT.1.OR.INTYP(1).GT.3)
     *             CALL QUITER(0.,INTYP(1),82,INERR)
        IF(IFREQ(1).LT.1.OR.IFREQ(1).GT.2)
     *             CALL QUITER(0.,IFREQ(1),101,INERR)
C
C  Convert the 2 waiver rates from percents to percentages for PCLEFT.
C
        WAIVER(1,1)=WAIVER(1,1)/100.
        WAIVER(2,1)=WAIVER(2,1)/100.
C
        DO 70 IVIM=1,4
        IF(ILDT(IVIM,1).LT.1.OR.ILDT(IVIM,1).GT.2)
     *    CALL QUITER(0.,ILDT(IVIM,1),72,INERR)
   70 CONTINUE
C
        IF(ITEST(1).LT.1.OR.ITEST(1).GT.4)
     *    CALL QUITER(0.,ITEST(1),73,INERR)
C
C  Check ICUTFL and if it indicates alternate cutpoints MOBILE5 copies
C  these values into CUTHC(1),CUTCO(1),and CUTNO(1).
C
      IF(ICUTFL.EQ.2) THEN
C
         CUTHC(1)=DCUTHC
 
         CUTCO(1)=DCUTCO
 
         CUTNO(1)=DCUTNO
C
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
