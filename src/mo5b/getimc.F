C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETIMC(INERR,*)
C
C  GETIMC obtains and verifies the I/M credits, depending on IMFLAG.
C  If IMFLAG is 2 then GETIMC reads in I/M program #1 and if IMFLAG is
C  3 then it reads in I/M program #1 credits and I/M program #2 credits.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS2/ IMFLAG
C    /IMPAR2/ ITEST
C    /IMPAR4/ CUTHC,CUTCO,CUTNO,SMARTL
C    /IMPAR6/ IFREQ,ICUTFL,IFILET,INTYP
C    /IOUCOM/ IOUIMD,IOUASK
C    /RSDINF/ RSDFLG
C    /TTCINF/ IHYB
C
C  Output on return:
C
C    common blocks:
C    /ICR4VA/ CRD4VA
C    /ICR4VB/ CRD4VB
C    /IMPAR7/ TNAME
C
C  Local array subscripts:
C
C  SMART(4)     -    SMART( ISMART )
C
C  Local variable / array dictionary:
C
C   Name   Type                      Description
C  ------  ----  ---------------------------------------------------------------
C  BUFFER   C    To keep "Smart line" record 
C  CSMART   C    See middle of subroutine
C  ITECH    I    See top of program
C  IP       I    See top of program
C  IBY      I    See top of program
C  IMPNO1   I    Indicates having already found first I/M credits
C                1=not found 2=found
C  IMPNO2   I    Indicates having already found second I/M credits
C                1=not found 2=found
C  IMSKIP   I    Counts how many lines to skip in the I/M file
C  SMART    I    See middle of subroutine
C
C  Notes:
C
C  IMPAR7 was added to pass to GETSC1 the non-default filename for the
C  TECH1&2 I/M credit data.  IMPAR4 was added to MOBILE5 version 2 to
C  pass user updated cutpoints from GETIMP and GETIM2.  The case of
C  identical test types, frequency  and cutpoints is handled
C  seperately and the two credit arrays are then meant to be identical
C
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added checking "smart line" if Hybrid or RSD program.
C
      CHARACTER*40 NAMIMC 
      CHARACTER*80 IMNAME
      INTEGER SMART(7)
      CHARACTER*1 CSMART(3)
      CHARACTER*120 BUFFER
c
cmh  LMOS modification from Mobile5a
cmh  Add definition of ENVL, ENVG, and NMIMC
c
      CHARACTER*132 ENVL,ENVG,NMIMC
cmh end
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'ICR4VA.I'
      INCLUDE 'ICR4VB.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IMPAR7.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'RSDINF.I'
      INCLUDE 'TTCINF.I' 

cmh add declarations
      LOGICAL WFLAG
cmh
cmh add commons
      COMMON /WARNING/ WFLAG
cmh

C
      IMPNO1 = 1
      IMPNO2 = 1
      SMARTL(1)=' '
      SMARTL(2)=' '
C
C  IMTIER flag signal the application or absence of the
C  I/M Special Credit for Tier 1 1994 - 1999+ LDV's
C
      IMTIER=0
      IF(IMFLAG.GT.3) THEN
         IMTIER=1
c
cmh  LMOS modification from Mobil5a
cmh  Comment out writing of message to IOUREP
c
c        IF(IMFLAG.EQ.4.OR.IMFLAG.EQ.5) WRITE(IOUERR,205)
c 205    FORMAT('0',
c    *   '(12/1/92) Guidance Memo Credits being used ',
c    *   'for Tier 1 vehicles')
         IF(IMFLAG.EQ.4) IMFLAG=2
         IF(IMFLAG.EQ.5) IMFLAG=3
      ENDIF
C
C     IF(IMFLAG.EQ.2) IMPNO2=2
C
C  User supplies filename when IFILET Flag says to
C
      IF(PROMPT.EQ.2.AND.IFILET(1).EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ',
     *  'Please enter the alternate TECH1&2 I/M credits filename:')
      IF(IFILET(1).EQ.2) READ(IOUGEN,100) TNAME
      IF(PROMPT.EQ.2.AND.IFILET(2).EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ',
     *  'Please enter the alternate TECH4+ I/M credits filename:')
C
      IF(IFILET(2).EQ.2) THEN
         READ(IOUGEN,100) NAMIMC
  100 FORMAT(A40)
c
cmh  LMOS Modification from Mobile5a
cmh  Read the user-supplied, alternate TECH4+ credit file from EMS_LOC.
cmh  Otherwise get the default file, IMDATA.D, from EMS_LOC.
c
         CALL GETENV("EMS_LOC",ENVL)
         NMIMC=ENVL(1:LNBLNK(ENVL))//'/'//NAMIMC
         OPEN(IOUIMD,FILE=NMIMC,STATUS='OLD')
cmh      OPEN(IOUIMD,FILE=NAMIMC,STATUS='OLD')
         GOTO 10
      ENDIF

cmh  Continue LMOS modification
cmh
      CALL GETENV("EMS_LOC",ENVG)
      NMIMC=ENVG(1:LNBLNK(ENVG))//'/IMDATA.D'
      OPEN(IOUIMD,FILE=NMIMC,STATUS='OLD')

cmh   NAMIMC='IMDATA.D'
cmh   OPEN(IOUIMD,FILE=NAMIMC,STATUS='OLD')
cmh end
C
C  Read the Smart control strip, cut point values and description line
C  Check technology and if it is incorrect then skip 48 lines
C
C
C  SMART(1) = tells the altitude 1=Low 2=High
C  SMART(2) = tells the frequency 1=Annual 2=Biennial of I/M program #1.
C  SMART(3) = tells the technology 1=TECH1&2 2=TECH4+
C  SMART(4) = tells the test type 1=Idle 2=2500/Idle 3=Loaded/Idle
C                                                    4=IM240 Test
C  CCUTHC   = Cutpoint for HC
C  CCUTCO   = Cutpoint for CO
C  CCUTNO   = Cutpoint for NOx
C  CSMART(1) or SMART(5) = tells if Hybrid program (seSET 
C  CSMART(2) or SMART(6)= tells if RSD program
C  CSMART(3) or SMART(7)= reserv 
C  IMNAME   = Description line
C
   
   10 READ(IOUIMD,109,END=95,ERR=96) BUFFER
  109 FORMAT (A120)
      READ(UNIT=BUFFER,FMT=110,ERR=96) (SMART(ISMART),ISMART=1,4),
     *    CCUTHC,CCUTCO,CCUTNO,(CSMART(ISMART),ISMART=1,3),IMNAME
  110 FORMAT(4I1,1X,F5.1,1X,F4.1,1X,F5.1,1X,3A1,1X,A80)
C  
C Check MOBILE5b or MOBILE5a format
C
      IF ( (CSMART(1) .EQ. '1') .OR. (CSMART(1) .EQ. '2') ) THEN
C
        DO 11 ISMART=1,3
C
          IF (CSMART(ISMART) .EQ. '2') THEN
            SMART(ISMART+4) = 2
          ELSE
            SMART(ISMART+4) = 1
          ENDIF
C
  11    CONTINUE
C
      ELSE
C
C Use the Mobile5a format
C
        READ(UNIT=BUFFER,FMT=111,ERR=96) (SMART(ISMART),ISMART=1,4),
     * CCUTHC,CCUTCO,CCUTNO,IMNAME
  111   FORMAT(4I1,1X,F5.1,1X,F4.1,1X,F5.1,1X,A80)
        SMART(5)=0
        SMART(6)=0
        SMART(7)=0
      ENDIF
C
      IF(SMART(3).NE.2) GOTO 70
C
C  Checks if I/M Program #1 has the correct frequency, test type and cutpoints
C  as read.  If so I/M Program #1 credit data are read and if not 96 lines are
C  skipped
C
      IF(IFREQ(1).NE.SMART(2).OR.
     *   ITEST(1).NE.SMART(4).OR.
     *   CUTHC(1).NE.CCUTHC.OR.
     *   CUTCO(1).NE.CCUTCO.OR.
     *   CUTNO(1).NE.CCUTNO.OR.
C  Check if Hybrid program
     *   ((IHYB(1) .EQ. 2) .AND. (SMART(5) .NE. 2)) .OR.
C  Check if RSD program     
     *   ((RSDFLG .EQ.2) .AND. (SMART(6) .NE. 2))) THEN
            IF(IMFLAG.EQ.3) GOTO 40
            GOTO 70
      ENDIF
C
cmh add conditional
      IF( WFLAG ) THEN
          WRITE(IOUASK,220) IMNAME
      ENDIF

  220 FORMAT( 7X,' Reading I/M credits information'/ 8X,A80)
C
C     Now that the data are being read in set the flag to indicate.
C
      IMPNO1=2
      SMARTL(1)=IMNAME
C
C     Loop for the 18 ITECH to be stored in the record.
C
      DO 20 ITECH=1,18
C
C       Loop for the three (3) pollutants in the record.
C
        DO 15 IP=1,3
          READ(IOUIMD,120,END=95,ERR=96)
     *        (CRD4VA(IBY,IP,ITECH),IBY=1,24)
   15   CONTINUE
   20 CONTINUE
C
C  First checks if there is a second I/M program and then if there is
C  I/M Program #2 is checked for frequency, test type and cutpoints.
C  If they are correct I/M Program #2 credit data is read and if not
C  96 lines are skipped
C
C  If there are two I/M programs, their test types, frequency and
C  cutpoints are identical then the IMPNO2 flag is turned on and the
C  arrays are madeto be identical and then GETIMC returns to ONESEC.
C
       IF(IMFLAG.EQ.2) GOTO 97
C
C  If there is not a 2nd I/M program and credits for the first program
C  have been read (which has happened if the program has reached this
C  spot) then the next thing to be done is to close the file and
C  return to ONESEC.
C
      IF(IFREQ(1).EQ.IFREQ(2).AND.CUTHC(1).EQ.CUTHC(2).AND.
     *   CUTCO(1).EQ.CUTCO(2).AND.ITEST(1).EQ.ITEST(2)) THEN
C
C          Loop for the 16 ITECH to be stored in the record.
C
           DO 35 ITECH=1,18
C
C            Loop for the three (3) pollutants in the record.
C
             DO 30 IP=1,3
C
C              Loop for the 24 years in the record.
C
               DO 25 IBY=1,24
C
C                Assign the second set of credits the same
C
                 CRD4VB(IBY,IP,ITECH)=CRD4VA(IBY,IP,ITECH)
   25          CONTINUE
   30        CONTINUE
C
   35      CONTINUE
           GOTO 97
      ENDIF
C
C  By this point the second I/M program credits may or may not have been read.
C  And if the have the only thing left is to close the file and return.  Or
C  they haven't and then the flow goes back to the read line.
C
      IF(IMPNO2.EQ.2) GOTO 97
      GOTO 10
   40 IF(IFREQ(2).NE.SMART(2).OR.
     *   ITEST(2).NE.SMART(4).OR.
     *   CUTHC(2).NE.CCUTHC.OR.
     *   CUTCO(2).NE.CCUTCO.OR.
     *   CUTNO(2).NE.CCUTNO.OR. 
C  Check if Hybrid program
     *   ((IHYB(2) .EQ. 2) .AND. (SMART(5) .NE. 2)) .OR.
C  Check if RSD program     
     *   ((RSDFLG .EQ.2) .AND. (SMART(6) .NE. 2)) ) GOTO 70
C
C  Set the read indicator flag
C
      IMPNO2=2
      SMARTL(2)=IMNAME

cmh add conditional
      IF( WFLAG ) THEN
          WRITE(IOUASK,230) IMNAME
      ENDIF

  230 FORMAT
     *( 7X,' Reading second I/M program credits information'/ 8X,A80)
C
C     Loop for the 18 ITECH to be stored in the record.
C
      DO 60 ITECH=1,18
C
C       Loop for the three (3) pollutants in the record.
C
        DO 50 IP=1,3
          READ(IOUIMD,120,END=95,ERR=96)
     *        (CRD4VB(IBY,IP,ITECH),IBY=1,24)
   50   CONTINUE
   60 CONTINUE
C
C  Check to see if first I/M credit have been read and if not go to
C  the read. Or if they have then return.
C
      IF(IMPNO1.EQ.1) GOTO 10
      GOTO 97
C
C  Skip through data file till next set of records.
C
   70 DO 80 IMSKIP=1,54
        READ(IOUIMD,130,END=95,ERR=96)
   80 CONTINUE
      GOTO 10
  120 FORMAT(24F5.3)
  130 FORMAT(1X)
   95 CALL QUITER(0.,0,54,INERR)
      GOTO 99
   96 CALL QUITER(0.,0,55,INERR)
      GOTO 99
   97 CLOSE(IOUIMD)
      RETURN
   99 RETURN1
      END
