C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /IMPAR4/ CUTHC,CUTCO,CUTNO,SMARTL,SMRTL1
C    /IMPAR6/ IFREQ,INTYP
C    /IOUCOM/ IOUREP
C    /TTCINF/ IHYB,IASM,ITTCFG,TTCEFF,TTCSML
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
C  ITST1    I    to index an array variables 
C
C  Notes:
C
C  Two paths this can go on.  If IMFLAG is 2 (two) then just 1 (one) program
C  is printed if IMFLAG is 3 (three)  then there are 2 (two) programs output.
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added printing "smart line I/M credit description"
C  Added printing "smart line Tech1&2 credit description"
C  Added printing "smart line TTC.D file"
C
C  9 Feb 96 @DynTel-MLA 2-611 Print the discounts for HC, CO, and NOx
C      that appply to the selected I/M program. The discounts are always
C      printed, whether they are defaults or user-supplied.
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'STRING.I'
      INCLUDE 'TTCINF.I'
      INTEGER     IP
      CHARACTER*4 TFREQ(2,2),TINTYP(7,4),TTYPE(4,5)
C
      DATA TTYPE/'Idle','    ','    ','    ',
     2           '2500',' rpm',' / I','dle ',
     3           'Load','ed /',' Idl','e   ',
     4           'IM24','0 te','st  ','    ',
     5           'ASM ','Test','    ','    '/
C
      DATA TFREQ/'Annu','al  ',
     *           'Bien','nial'/
C
      DATA TINTYP/'Test',' Onl','y   ','    ','    ','    ','    ',
     2            'Comp','uter','ized',' Tes','t an','d Re','pair',
     3            'Manu','al T','est ','and ','Repa','ir  ','    ',
     4            'Rete','st B','ased',' Hyb','rid ','Prog','ram '/
C
      WAIVER(1,1)=WAIVER(1,1)*100.
      WAIVER(2,2)=WAIVER(2,2)*100.
      WAIVER(2,1)=WAIVER(2,1)*100.
      WAIVER(1,2)=WAIVER(1,2)*100.
C
      IF(IMFLAG.EQ.3) GOTO 60
C      
C  One I/M program      
C
      WRITE(IOUREP,200)
  200 FORMAT('0','I/M program selected:'/)
C
      ITST1=INTYP(1)
      IF(IHYB(1).EQ.2) ITST1=4
C
      WRITE(IOUREP,215) ICYIM(1),ISTRIN(1), MODYR(1,1), MODYR(2,1),
     *     WAIVER(1,1), WAIVER(2,1), CRIM(1), (TINTYP(J,ITST1),J=1,7)
  215 FORMAT('0','    Start year (January 1):             ',I4/
     *       ' ','    Pre-1981 MYR stringency rate:        ',I2,'%'/
     *       ' ','    First model year covered:           ',I4/
     *       ' ','    Last model year covered:            ',I4/
     *       ' ','    Waiver rate (pre-1981):             ',F3.0,'%'/
     *       ' ','    Waiver rate (1981 and newer):       ',F3.0,'%'/
     *       ' ','    Compliance Rate:                   ',F4.0,'%'/
     *       ' ','    Inspection type:                    ',7A4)
C
C     Print the program effectiveness rates if the program is either
C     a manual or a computerized test-and-repair program.
C
      IF (INTYP(1).EQ.2.OR.INTYP(1).EQ.3) THEN
        WRITE(IOUREP,216) (DISCNT(IP, INTYP(1)), IP=1,3)
  216   FORMAT(' ','    Effectiveness - HC: ', F4.2, '  CO: ', F4.2, 
     *             '  NOx: ', F4.2)
      ENDIF
C
      WRITE(IOUREP,218) (TFREQ(I,IFREQ(1)),I=1,2)
  218 FORMAT(' ','    Inspection frequency                ',2A4)
C
      WRITE(IOUREP,219) (NOYES(ILDT(IVIM,1)),IVIM=1,4)
  219 FORMAT(' ','    Vehicle types covered:              LDGV - ',A4/
     * ' ',T41,'LDGT1 - ',A4/' ',T41,'LDGT2 - ',A4/' ',T42,'HDGV - ',A4)
C
      ITST1=ITEST(1)
      IF(IASM(1).EQ.2) ITST1=5
      WRITE(IOUREP,220) (TTYPE(ICH,ITST1),ICH=1,4)
  220 FORMAT(' ','    1981 & later MYR test type:         ',4A4)
      WRITE(IOUREP,223) CUTHC(1),CUTCO(1),CUTNO(1)
  223 FORMAT(' ','    Cutpoints, HC: ',F8.3,'   CO: ',
     * F8.3, '   NOx: ', F8.3)
C     
C  Print "smart line I/M and Tech1&2  credit description"
C
      IF (SMRTL1 .NE. ' ') WRITE(IOUREP,224) SMRTL1
      WRITE(IOUREP,224) SMARTL(1)                  
  224 FORMAT(' ',A80)
C
      GOTO 70
C
C Two I/M programs
C
   60 WRITE(IOUREP,225)
  225 FORMAT('0','I/M program #1 selected:                          ',
     *           'I/M program #2 selected:'/)
C
      ITST1=INTYP(1)
      IF(IHYB(1).EQ.2) ITST1=4
      ITST2=INTYP(2)
      IF(IHYB(2).EQ.2) ITST2=4
C
      WRITE(IOUREP,230) ICYIM(1),ICYIM(2),
     *   ISTRIN(1),ISTRIN(2),
     *   MODYR(1,1),MODYR(1,2),
     *   MODYR(2,1),MODYR(2,2),
     *   WAIVER(1,1),WAIVER(1,2),
     *   WAIVER(2,1),WAIVER(2,2),
     *   CRIM(1),CRIM(2),
     *   (TINTYP(J,ITST1),J=1,7),(TINTYP(J,ITST2),J=1,7)
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
     *         'Inspection type:     '/1X,7A4,22X,7A4)
C
C     Print effectiveness ratios if either program is test-and-repair.
C
      IF (INTYP(1).EQ.2.OR.INTYP(1).EQ.3.OR.
     *    INTYP(2).EQ.2.OR.INTYP(2).EQ.3) THEN
        IF (INTYP(1).EQ.2.OR.INTYP(1).EQ.3) THEN
C         First program is test-and-repair...
          IF (INTYP(2).EQ.2.OR.INTYP(2).EQ.3) THEN
C           Both programs are test-and-repair.
            WRITE(IOUREP,232) (DISCNT(IP,INTYP(1)),IP=1,3),
     *                        (DISCNT(IP,INTYP(2)),IP=1,3)
  232       FORMAT(' ','Effectiveness - HC: ', F4.2, '  CO: ', F4.2,
     *                 '  NOx: ', F4.2, '     ',
     *                 'Effectiveness - HC: ', F4.2, '  CO: ', F4.2, 
     *                 '  NOx: ', F4.2)
          ELSE
C           First program is test-and-repair, second is not.
            WRITE (IOUREP,234) (DISCNT(IP,INTYP(1)),IP=1,3)
  234       FORMAT(' ','Effectiveness - HC: ', F4.2, '  CO: ', F4.2,
     *                 '  NOx: ', F4.2, '     ')
          ENDIF
        ELSE
C         First program is not test-and-repair, but second is.
          WRITE (IOUREP,236) (DISCNT(IP,INTYP(2)),IP=1,3)
  236     FORMAT(T52,'Effectiveness - HC: ', F4.2, '  CO: ', F4.2,
     *               '  NOx: ', F4.2, '     ')
        ENDIF
      ENDIF
C
      WRITE(IOUREP,238) (TFREQ(I,IFREQ(1)),I=1,2),
     *                  (TFREQ(I,IFREQ(2)),I=1,2)
  238 FORMAT(' ','Inspection frequency: ',2A4,'                    '
     *          ,'Inspection frequency: ',2A4)
C
      WRITE(IOUREP,240) ((NOYES(ILDT(IVIM,IMF)),IMF=1,2),IVIM=1,4)
  240 FORMAT(' ','I/M program #1 vehicle types',
     *           22X,'I/M program #2 vehicle types'/
     *  T7, 'LDGV - ', A4, T57,  'LDGV - ',A4/
     * T6, 'LDGT1 - ', A4, T56, 'LDGT1 - ',A4/
     * T6, 'LDGT2 - ', A4, T56, 'LDGT2 - ',A4/
     *  T7, 'HDGV - ', A4, T57,  'HDGV - ',A4)
C
      ITST1=ITEST(1)
      IF(IASM(1).EQ.2) ITST1=5
      ITST2=ITEST(2)
      IF(IASM(2).EQ.2) ITST2=5
      WRITE(IOUREP,250) (TTYPE(ICH,ITST1),ICH=1,4),
     *                  (TTYPE(ICH,ITST2),ICH=1,4)
C
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
C  Print "smart line I/M and Tech1&2   credit description"
C
      IF (SMRTL1 .NE. ' ') WRITE(IOUREP,224) SMRTL1
      WRITE(IOUREP,224) SMARTL(1)                  
      IF ((SMARTL(2) .NE. ' ') .AND. (SMARTL(1) .NE. SMARTL(2)))
     *  WRITE (IOUREP,224) SMARTL(2) 
C
C  Print "smart line" from TTC.D file and percent of effectiveness. 
C
   70 IF (ITTCFG .GT. 1) WRITE(IOUREP,270) TTCEFF*100
 270  FORMAT (1X,'With ',F5.1,'% Technician Training and Certification',
     *        ' Credit')
      IF (ITTCFG .EQ. 5 .AND. TTCSML .NE. ' ')
     *       WRITE(IOUREP,224) TTCSML   
C
      WAIVER(1,1)=WAIVER(1,1)/100.
      WAIVER(2,2)=WAIVER(2,2)/100.
      WAIVER(2,1)=WAIVER(2,1)/100.
      WAIVER(1,2)=WAIVER(1,2)/100.
C
   90 CONTINUE
C
   99 RETURN
      END
