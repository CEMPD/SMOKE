C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTIMP
C
C  OUTIMP prints the I/M program parameters.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT,ITEST
C    /IMPAR4/ CUTHC,CUTCO,CUTNO
C    /IMPAR6/ IFREQ,INTYP
C    /IOUCOM/ IOUREP
C    /STRING/ NOYES
C
C  Local array subscripts:
C
C  TFREQ(2,2)   -  TFREQ ( ICH, IFREQ(IMPNUM) )
C  TINTYP(7,3)  -  TFREQ ( ICH, INTYP(IMPNUM) )
C  TTYPE(4,3)   -  TTYPE ( ICH, ITEST(IMPNUM) )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  TFREQ    C*4  I/M inspection frequency description
C  TINTYP   C*4  I/M program type description
C  TTYPE    C*4  I/M test type description
C
C  Notes:
C
C  Two paths this can go on.  If IMFLAG is 2 (two) then just 1 (one) program
C  is printed if IMFLAG is 3 (three)  then there are 2 (two) programs output.
C
C
      CHARACTER*4 NOYES,COMMA,PERIOD
      CHARACTER*1 COLON
C
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR4/ CUTHC(2),CUTCO(2),CUTNO(2)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
C
      CHARACTER*4 TFREQ(2,2),TINTYP(7,3),TTYPE(4,4)
C
      DATA TTYPE/'Idle','    ','    ','    ',
     2           '2500',' rpm',' / I','dle ',
     3           'Load','ed /',' Idl','e   ',
     4           'IM24','0 te','st  ','    '/
C
      DATA TFREQ/'Annu','al  ',
     *           'Bien','nial'/
C
      DATA TINTYP/'Test',' Onl','y   ','    ','    ','    ','    ',
     2            'Comp','uter','ized',' Tes','t an','d Re','pair',
     3            'Manu','al T','est ','and ','Repa','ir  ','    '/
C
      WAIVER(1,1)=WAIVER(1,1)*100.
      WAIVER(2,2)=WAIVER(2,2)*100.
      WAIVER(2,1)=WAIVER(2,1)*100.
      WAIVER(1,2)=WAIVER(1,2)*100.
C
      IF(IMFLAG.EQ.3) GOTO 60
      WRITE(IOUREP,200)
  200 FORMAT('0','I/M program selected:'/)
      WRITE(IOUREP,215) ICYIM(1),ISTRIN(1),
     *     MODYR(1,1), MODYR(2,1),
     *     WAIVER(1,1),WAIVER(2,1),CRIM(1),
     *                 (TINTYP(J,INTYP(1)),J=1,7),
     *                 (TFREQ(I,IFREQ(1)),I=1,2)
  215 FORMAT('0','    Start year (January 1):             ',I4/
     *       ' ','    Pre-1981 MYR stringency rate:        ',I2,'%'/
     *       ' ','    First model year covered:           ',I4/
     *       ' ','    Last model year covered:            ',I4/
     *       ' ','    Waiver rate (pre-1981):             ',F3.0,'%'/
     *       ' ','    Waiver rate (1981 and newer):       ',F3.0,'%'/
     *       ' ','    Compliance Rate:                   ',F4.0,'%'/
     *       ' ','    Inspection type:                    ',7A4,/
     *       ' ','    Inspection frequency                ',2A4)
      WRITE(IOUREP,216) (NOYES(ILDT(IVIM,1)),IVIM=1,4)
  216 FORMAT(' ','    Vehicle types covered:              LDGV - ',A4/
     * ' ',T41,'LDGT1 - ',A4/' ',T41,'LDGT2 - ',A4/' ',T42,'HDGV - ',A4)
      WRITE(IOUREP,220) (TTYPE(ICH,ITEST(1)),ICH=1,4)
  220 FORMAT(' ','    1981 & later MYR test type:         ',4A4)
      WRITE(IOUREP,223) CUTHC(1),CUTCO(1),CUTNO(1)
  223 FORMAT(' ','    Cutpoints, HC: ',F8.3,'   CO: ',
     * F8.3, '   NOx: ', F8.3)
C
      GOTO 70
C
   60 WRITE(IOUREP,225)
  225 FORMAT('0','I/M program #1 selected:                          ',
     *           'I/M program #2 selected:'/)
      WRITE(IOUREP,230) ICYIM(1),ICYIM(2),
     *   ISTRIN(1),ISTRIN(2),
     *   MODYR(1,1),MODYR(1,2),
     *   MODYR(2,1),MODYR(2,2),
     *   WAIVER(1,1),WAIVER(1,2),
     *   WAIVER(2,1),WAIVER(2,2),
     *   CRIM(1),CRIM(2),
     *   (TINTYP(J,INTYP(1)),J=1,7),(TINTYP(J,INTYP(2)),J=1,7),
     *   (TFREQ(I,IFREQ(1)),I=1,2),(TFREQ(I,IFREQ(2)),I=1,2)
  230 FORMAT('0', 'Start year (Jan 1):  ',I4,'                         '
     *           ,'Start year (Jan 1):  ',I4/
     *   ' ','Pre-1981 stringency: ',I2,'%','                          '
     *      ,'Pre-1981 stringency: ',I2,'%'/
     *        ' ','First MYR covered:   ',I4,'                         '
     *           ,'First MYR covered:   ',I4/
     *        ' ','Last MYR covered:    ',I4,'                         '
     *           ,'Last MYR covered:    ',I4/
     *  ' ','Waiver (pre-1981):   ',F3.0,'%','                         '
     *     ,'Waiver (pre-1981):   ',F3.0,'%'/
     *  ' ','Waiver (1981+):      ',F3.0,'%','                         '
     *     ,'Waiver (1981+):      ',F3.0,'%'/
     *   ' ','Compliance Rate:    ',F4.0,'%','                         '
     *      ,'Compliance Rate:    ',F4.0,'%'/
     *     ' ','Inspection type:     ',29X,
     *         'Inspection type:     '/1X,7A4,22X,7A4/
     *           ' ','Inspection frequency: ',2A4,'                    '
     *              ,'Inspection frequency: ',2A4)
      WRITE(IOUREP,240) ((NOYES(ILDT(IVIM,IMF)),IMF=1,2),IVIM=1,4)
  240 FORMAT(' ','I/M program #1 vehicle types',
     *           22X,'I/M program #2 vehicle types'/
     *  T7, 'LDGV - ', A4, T57,  'LDGV - ',A4/
     * T6, 'LDGT1 - ', A4, T56, 'LDGT1 - ',A4/
     * T6, 'LDGT2 - ', A4, T56, 'LDGT2 - ',A4/
     *  T7, 'HDGV - ', A4, T57,  'HDGV - ',A4)
      WRITE(IOUREP,250) (TTYPE(ICH,ITEST(1)),ICH=1,4),
     *                  (TTYPE(ICH,ITEST(2)),ICH=1,4)
  250 FORMAT(' ','1981 & later MYR test type:', 23X,
     *           '1981 & later MYR test type:'/6X,4A4,34X,4A4)
      WRITE(IOUREP,260) CUTHC,CUTCO,CUTNO
  260 FORMAT(' ','Cutpoints, HC:  ',F8.3,26X,
     *           'Cutpoints, HC:  ',F8.3/
     *       ' ','Cutpoints, CO:  ',F8.3,26X,
     *           'Cutpoints, CO:  ',F8.3/
     *       ' ','Cutpoints, NOx: ',F8.3,26X,
     *           'Cutpoints, NOx: ',F8.3/)
C
   70 WAIVER(1,1)=WAIVER(1,1)/100.
      WAIVER(2,2)=WAIVER(2,2)/100.
      WAIVER(2,1)=WAIVER(2,1)/100.
      WAIVER(1,2)=WAIVER(1,2)/100.
C
   90 CONTINUE
C
   99 RETURN
      END
