C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /FLAGS2/ IMFLAG 
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IYSTART,IY1960,IYEND
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT,ITEST
C    /IMPAR4/ CUTHC,CUTCO,CUTNO
C    /IMPAR6/ IFREQ,INTYP,IFILET, IMEFFR
C    /TTCINF/ ITTCFG,TTCEFF,IHYB,IASM
C    /RSDINF/ RSDFLG
C
C  Local variable / array dictionary:
C
C   Name   Type                      Description
C  ------  ----  ---------------------------------------------------------------
C  DCTFLG   C    Flag - user will enter effectiveness rates for I/M pgms.
C  IVIM     I    See top of program
C  ICUTFL   I    Indicates alternate cutpoints when equal to 2
C  DCUTHC,DCUTCO,DCUTNO
C           R    Temporary variables to hold cutpoints
C  NUMIMP   I    Number I/M programs
C  TI1FLG   I    Yes/No Tier1 program
C  RDTTCF   I    Read or don't read TTC.D file
C
C  Notes:
C
C  GETIMP was modified in MOBILE4.1 to exclude the IM240 program.
C  and again in MOBILE5 v2.0 to check for new user inputs in the I/M parameter
C  these new inputs allow the user to specify their own cutpoints, and filenames
C  or they can mix and match between their own values and default MOBILE5 values
C  Note TEST TYPE is right next to new flags in MOBILE5, this is done to allow
C  the I/M parameter record to be more compatible with MOBILE4.1, recall that
C  there once were 2 (two) indicators of alternate I/M credits located at the
C  end of the I/M parameter record and they were in the following FORTRAN
C  format 1X, 2I1
C  December-14-1994 @CSC-zp Request 2-468 Changed the second parameter to 
C  integer value in call to QUITER subroutine according to parameter's format
C  in the subroutine itself.
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added reading new records when IMFLAG=6 (TTC information).  
C  Modified the "Program type" and "Test type" checking.
C  Set values of IHYB and IASM.
C  Initialize flags for processing TTC and RSD programs
C
C  8 Feb 96 @DynTel-MLA 2-611 Added DCTFLG to indicate whether user will
C     enter alternate discount rates for test-and-repair I/M programs.
C 
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'RSDINF.I' 
      INCLUDE 'TTCINF.I'
      INCLUDE 'YEARS4.I'
C
C
      CHARACTER DCTFLG*1
      INTEGER UNIT
      INTEGER NUMIMP,RDTTCF,TI1FLG
C
C Variables for the my VS cut-points test
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
C Initialize ICUTFL to indicate use default cutpoints.
C
      ICUTFL=1
C
C Initialize flags for processing TTC and RSD programs
C
      RSDFLG=1
      ITTCFG=1
      TTCEFF=1.0
C
C Read records with TTC, RSD control information
C
      IF (IMFLAG .EQ. 6) THEN
        RDTTCF=1
        READ(IOUGEN,105,ERR=98,END=98) NUMIMP,TI1FLG,ITTCFG,RSDFLG
 105    FORMAT(4(I1,1X))
C
C Check for the valid values
C
        IF ( (NUMIMP .LT. 1) .OR. (NUMIMP .GT. 2) ) THEN
          CALL QUITER(22.,NUMIMP,1,INERR)
          GOTO 998
        ENDIF
C
        IF ( (TI1FLG .LT. 1) .OR. (TI1FLG .GT. 2) ) THEN
          CALL QUITER(23.,TI1FLG,1,INERR)
          GOTO 998
        ENDIF
C
        IF ( (ITTCFG .LT. 1) .OR. (ITTCFG .GT. 3) ) THEN
          CALL QUITER(24.,ITTCFG,1,INERR)
          GOTO 998
        ENDIF
C
        IF ( (RSDFLG .LT. 1) .OR. (RSDFLG .GT. 2) )  THEN
          CALL QUITER(25.,RSDFLG,1,INERR)
          GOTO 998
        ENDIF
C
        IF (ITTCFG .EQ. 3 ) THEN
          READ(IOUGEN,106,ERR=98,END=98) RDTTCF,TTCEFF
 106      FORMAT(I1,1X,F5.2)
C
          IF ( (RDTTCF .LT. 1) .OR. (RDTTCF .GT. 2) )  THEN
            CALL QUITER(26.,RDTTCF,1,INERR)
            GOTO 998
          ENDIF
C
          IF ( (TTCEFF .LT. 0.0) .OR. (TTCEFF .GT. 1.0) )  THEN
            CALL QUITER(TTCEFF,0,165,INERR)
            TTCEFF=1.0
          ENDIF
C
C ITTCFG=4 use the user's percent of the TTC effectiveness and default
C          values for TTC program.
C ITTCFG=5 use the user's percent of the TTC effectiveness and read
C          values for TTC program from the TTC.D file.
C
          IF (RDTTCF .EQ. 1) THEN
            ITTCFG=4
          ELSE
            ITTCFG=5
          ENDIF
C
        ENDIF
C        
C Create IMFLAG for further processing
C
        IMFLAG=1+NUMIMP
C
        IF (TI1FLG .EQ. 2) THEN
          IMFLAG=IMFLAG+2
        ENDIF
C
      ENDIF      
C
C User Inputted Values
C
      IF(PROMPT.EQ.2) THEN
C
        WRITE(IOUASK,200)
  200   FORMAT(' ','Enter the I/M Program parameters record:'/
     *         '0','Program start year, stringency,'/
     *         ' ','first and last model year getting benefits,'/
     *         ' ','old tech waiver rate, new tech waiver rate,'/
     *         ' ','compliance rate,'/
     *         ' ','inspection type,'/
     *         ' ','flag for effectiveness rates (optional),'/
     *         ' ','frequency of I/M inspection,'/
     *         ' ','the format is:'/
     *         ' ','(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1)')
C
        READ(IOUGEN,100,ERR=99,END=99)
     * ICYIM(1),ISTRIN(1),MODYR(1,1),MODYR(2,1),
     * WAIVER(1,1),WAIVER(2,1),CRIM(1),INTYP(1),DCTFLG,IFREQ(1)
  100   FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1)
C
        WRITE(IOUASK,210)
  210   FORMAT('0','Enter the rest of the I/M Program parameters '
     *            ,'record:'/
     *         '0','vehicle classes covered, test type,'/
     *         ' ','cutpoint flag (optional),'/
     *         ' ','alternate filename flags (optional)'/
     *         ' ','and any alternate cutpoints (contingent)'/
     *         ' ','the format is:'/
     *         ' ','4I1,1X,I1,I1(optional),I1(optional),'
     *            ,'I1(optional),3(1X,F4.2)(contingent)')
C
C
        READ(IOUGEN,110,ERR=99,END=99)
     * (ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     * ICUTFL,IFILET(1),IFILET(2),DCUTHC,DCUTCO,DCUTNO
  110   FORMAT(4I1,1X,4I1,3(1X,F4.2))
C
        GOTO 25
      ENDIF
C
C MOBILE5 uses a different format for interactive input because
C the format description for interactive input is too long to fit
C on one line.
C
      READ(IOUGEN,120,ERR=99,END=99)
     * ICYIM(1),ISTRIN(1),MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),
     * CRIM(1),INTYP(1),DCTFLG,IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),
     * ITEST(1),ICUTFL,IFILET(1),IFILET(2)
  120 FORMAT(4(I2,1X),2(F2.0,1X),F3.0,1X,I1,A1,I1,1X,4I1,1X,4I1)
C
C Try and find The *OLD* M4.1 YEARS
C
      BACKSPACE(UNIT)
      READ(IOUGEN,101) CM24YR,CPRGYR,CPRSYR
 101  FORMAT(T36,3(1X,A2))
C
C Read the character array CH
C
      BACKSPACE(UNIT)
      READ(IOUGEN,102) (CH(J),J=1,15)
 102  FORMAT(T36,15(A1))
C
C Try and find a "." period for floating point number;
C if true, go and reread the I/M input line
C
      DO 15 J=1,15
C
        IF (CH(J).EQ.'.') GOTO 20
C
 15   CONTINUE
C
C Check for three years and/or 15 blanks in a row
C
      IF ( (CM24YR.EQ.'  '.AND.CPRGYR.EQ.'  '.AND.CPRSYR.EQ.'  ').OR.
     *     (CH(1).EQ.' '.AND.CH(2).EQ.' '.AND.CH(3).EQ.' '
     * .AND.CH(4).EQ.' '.AND.CH(5).EQ.' '.AND.CH(6).EQ.' '
     * .AND.CH(7).EQ.' '.AND.CH(8).EQ.' '.AND.CH(9).EQ.' '
     * .AND.CH(10).EQ.' '.AND.CH(11).EQ.' '.AND.CH(12).EQ.' '
     * .AND.CH(13).EQ.' '.AND.CH(14).EQ.' '.AND.CH(15).EQ.' ') ) GOTO 25
C
C Check for a character string or numbers;
C if only characters 'NO' then M152 Warning ...
C
C 1st M4.1 yr: CM24YR
C
      IF (CH(1).EQ.' '.AND.CM24YR.NE.'  ') THEN
C
        DO 60 I=1,10
C
          IF (CH(2).EQ.CHT(I)) THEN
C
            DO 61 J=1,10
C
              IF (CH(3).EQ.CHT(I)) THEN
                CALL QUITER(0.,0,152,INERR)
                GOTO 25
              ENDIF
C
 61         CONTINUE
C
          ENDIF
C
 60     CONTINUE
C
      ENDIF
C
C 2nd M4.1 yr: CPRGYR
C
      IF (CH(4).EQ.' '.AND.CPRGYR.NE.'  ') THEN
C
        DO 64 I=1,10
C
          IF (CH(5).EQ.CHT(I)) THEN
C
            DO 65 J=1,10
C
              IF (CH(6).EQ.CHT(I)) THEN
                CALL QUITER(0.,0,152,INERR)
                GOTO 25
              ENDIF
C
 65         CONTINUE
C
          ENDIF
C
 64     CONTINUE
C
      ENDIF
C
C 3rd M4.1 yr: CPRSYR
C
      IF (CH(7).EQ.' '.AND.CPRSYR.NE.'  ') THEN
C
        DO 68 I=1,10
C
          IF (CH(8).EQ.CHT(I)) THEN
C
            DO 69 J=1,10
C
              IF (CH(9).EQ.CHT(I)) THEN
                CALL QUITER(0.,0,152,INERR)
                GOTO 25
              ENDIF
C
 69         CONTINUE
C
          ENDIF
C
 68     CONTINUE
C
      ENDIF
C
C Skip the Read
C
      GOTO 25
C
C Read in the M5.0 F4.2 Cut-Points
C
 20   CONTINUE
C
      BACKSPACE (UNIT)
      READ(IOUGEN,121,ERR=99,END=99) DCUTHC,DCUTCO,DCUTNO
 121  FORMAT(T36,3(1X,F4.2))
C
C Check I/M inputs.
C
 25   CONTINUE
C
      CALL YRTEST(ICYIM(1),32,IY1960,IYEND,INERR)
C
      IF(ISTRIN(1).LT.10.OR.ISTRIN(1).GT.50)
     *CALL QUITER(0.,ISTRIN(1),33,INERR)
C
      CALL YRTEST(MODYR(1,1),35,IYSTART,IYEND,INERR)
      CALL YRTEST(MODYR(2,1),36,IYSTART,IYEND,INERR)
      IF(MODYR(1,1).GT.MODYR(2,1)) CALL QUITER(0.,0,37,INERR)
C
      IF(WAIVER(1,1).LT.0.0.OR.WAIVER(1,1).GT.50.)
     *CALL QUITER(WAIVER(1,1),0,100,INERR)
C
      IF(WAIVER(2,1).LT.0.0.OR.WAIVER(2,1).GT.50.)
     *CALL QUITER(WAIVER(2,1),0,100,INERR)
C
      IF(CRIM(1).LT.0.0.OR.CRIM(1).GT.100.)
     *CALL QUITER(CRIM(1),0,34,INERR)
C
      IF(INTYP(1).LT.1.OR.INTYP(1).GT.4)
     *CALL QUITER(0.,INTYP(1),82,INERR)
C
      IF(DCTFLG.EQ.' '.OR.DCTFLG.EQ.'1') THEN
        IMEFFR = 1
      ELSEIF(DCTFLG.EQ.'2')THEN
        IMEFFR = 2
      ELSE
        CALL QUITER(0.0, 0, 102, INERR)
      ENDIF
C
      IF(IFREQ(1).LT.1.OR.IFREQ(1).GT.2)
     *CALL QUITER(0.,IFREQ(1),101,INERR)
C
C Convert the 2 waiver rates from percents to percentages for PCLEFT.
C
      WAIVER(1,1)=WAIVER(1,1)/100.
      WAIVER(2,1)=WAIVER(2,1)/100.
C
      DO 70 IVIM=1,4
        IF(ILDT(IVIM,1).LT.1.OR.ILDT(IVIM,1).GT.2)
     * CALL QUITER(0.,ILDT(IVIM,1),72,INERR)
   70 CONTINUE
C
      IF(ITEST(1).LT.1.OR.ITEST(1).GT.5)
     *CALL QUITER(0.,ITEST(1),73,INERR)
C
C Check ICUTFL and if it indicates alternate cutpoints MOBILE5 copies
C these values into CUTHC(1),CUTCO(1),and CUTNO(1).
C
      IF(ICUTFL.EQ.2) THEN
        CUTHC(1)=DCUTHC
        CUTCO(1)=DCUTCO
        CUTNO(1)=DCUTNO
      ENDIF
C
C User specifies a Retest Based Hybrid Program Design
C
      IHYB(1)=1
C
      IF(INTYP(1).EQ.4) THEN
        INTYP(1)=1
        IHYB(1)=2
      ENDIF
C
C User specifies the ASM test procedure.
C
      IASM(1)=1
C
      IF(ITEST(1).EQ.5) THEN
        ITEST(1)=3
        IASM(1)=2
      ENDIF
C
      RETURN
C
C EOF on any attempted read => take alternate return 1 => run aborts.
C
   98 CALL QUITER(0.,0,164,INERR)
      RETURN1
C      
   99 CALL QUITER(0.,0,129,INERR)
      RETURN1
C
  998 CALL QUITER(0.,0,53,INERR)
      RETURN1
C
      END
