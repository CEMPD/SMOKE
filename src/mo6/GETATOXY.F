      SUBROUTINE GETATOXY(INPREC,INERR,RC) 
C 
C  GETATOXY processes the Air Toxic OXYGENATE labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls UCCOMP,NXTTOK,NXTREC,QUITER 
C 
C  Changes: (Last change first) 
C 
C 26 Jun 01 AIR Task 15: New subroutine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /ATOX1/ ATOXY 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  NDATA      I   Loop control variable. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  TESTFIELD  C   Upper case version of the program type input field. 
C  SUM        R   Temporary summation variable 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX2.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'OXY1.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      CHARACTER*30   TESTFIELD 
      CHARACTER*150  M6REC 
      INTEGER        I 
      INTEGER        IOS 
      INTEGER        ITYPE 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           MKTSH 
      REAL           OXY 
      REAL           SUM 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Initialize the ATOXY array. 
C 
      ATOXY=0. 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 is set by the NXTTOK subroutine.  
C 
C     Loop for the 4 oxygenate types (MTBE, ETBE, Ethanol, TAME). 
C 
      DO I=1,4 
C 
      IF(I.EQ.1) THEN 
        M6REC=INPREC      ! First record already read via OXYGENATE command 
      ELSE 
        CALL NXTREC(IOUGEN,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) GOTO 85 
      ENDIF 
C 
      PTR2 = 21 
C 
C     Obtain and verify the oxygenate type 
C 
      CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
C 
      IF (TESTFIELD.EQ.'MTBE') THEN 
        ITYPE = 1 
      ELSEIF (TESTFIELD.EQ.'ETBE') THEN 
        ITYPE = 2 
      ELSEIF (TESTFIELD.EQ.'ETHANOL' .OR. 
     *        TESTFIELD.EQ.'ETOH') THEN 
        ITYPE = 3 
      ELSEIF (TESTFIELD.EQ.'TAME') THEN 
        ITYPE = 4 
      ELSE 
        CALL QUITER(0.0,0,571,INERR) 
        GOTO 95 
      ENDIF 
C 
C     Obtain and verify the vol% oxygen content. 
C 
      CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) OXY 
C 
      IF(OXY.LT.0. .OR. OXY.GT.5.) THEN 
        CALL QUITER(OXY,0,572,INERR) 
        GOTO 98 
      ENDIF 
C 
C     Obtain and verify the market share. 
C 
      CALL NXTTOK(M6REC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,IOSTAT=IOS) MKTSH 
C 
      IF(MKTSH.LT.0. .OR. MKTSH.GT.1.) THEN 
        CALL QUITER(MKTSH,0,573,INERR) 
        GOTO 98 
      ENDIF 
C 
C     Store the inputs 
C 
      ATOXY(3,ITYPE) = OXY     ! vol% 
      ATOXY(2,ITYPE) = MKTSH 
C 
      END DO 
C 
C  Determine the market share for no oxygenate fuels. 
C  Catch error that sum of 4 market share exceeds 100% 
C 
      SUM=0. 
      DO I=1,4 
        SUM=SUM+ATOXY(2,I) 
      END DO 
C 
      IF(SUM.GT.1.0) THEN 
        CALL QUITER(SUM,0,574,INERR) 
        GOTO 98 
      ENDIF 
C 
      ATOXY(2,5)=1.0 - SUM 
C 
C     Show that the 4 required OXYGENATE data were entered  
C     for this Air Toxic program. 
C 
      TOXCHK(6) = 2 
C 
C     Since the OXYGENATE command supercedes the OXYGENATED FUELS 
C     command, replace the market share and average oxygenate values 
C     which would have been read in. The input vol% oxygenate values 
C     must be changed to wt%, then to a fractional weight. 
C     Note that Ethers = MTBE + ETBE + TAME. 
C 
      SHRMKT(2) = ATOXY(2,1) + ATOXY(2,2) + ATOXY(2,4) !Ethers 
      SHRMKT(3) = ATOXY(2,3)                           !Alcohols 
C 
      MKTSH = SHRMKT(2) + SHRMKT(3) 
      IF(MKTSH .GT. 1.0) CALL QUITER(MKTSH,0,4,INERR) 
      SHRMKT(1)=1.0-SHRMKT(2)-SHRMKT(3) 
C 
      SUM = 0. 
      DO I=1,4 
        ATOXY(1,I) = ATOXY(3,I) * VOLTOWT(I) ! Change from vol% to wt% 
        IF(I.NE.3) SUM = SUM + ATOXY(1,I) * ATOXY(2,I)  ! wt% 
      END DO 
C 
      IF(SHRMKT(2) .GT. 0.) THEN 
        OXYCNT(1) = SUM / (SHRMKT(2)*100.) !Weighted Ether blend oxygen content, change wt% to fraction 
      ELSE 
        OXYCNT(1) = 0.0 
      ENDIF 
      OXYCNT(2) = ATOXY(1,3) / 100. !Alcohol oxygen content, change from wt% to fraction 
C 
C     Verify values using OXYGENATED FUELS command criteria. If 
C     outside range, warn user and reset only the OXYGENATED FUELS 
C     values, NOT the OXYGENATE values. 
C 
      IF(OXYCNT(1).LT.0..OR.OXYCNT(1).GT..027001) THEN 
        CALL QUITER(OXYCNT(1),0,5,INERR) 
        OXYCNT(1)=MAX(0.,MIN(.027,OXYCNT(1))) 
      ENDIF 
      
      IF(OXYCNT(2).LT.0..OR.OXYCNT(2).GT..035001) THEN 
        CALL QUITER(OXYCNT(2),0,6,INERR) 
        OXYCNT(2)=MAX(0.,MIN(.035,OXYCNT(2))) 
      ENDIF 
C 
      IF(SHRMKT(2).GT.0..AND.OXYCNT(1).EQ.0. .OR. 
     *   SHRMKT(2).EQ.0..AND.OXYCNT(1).GT.0.) 
     *   CALL QUITER(0.,0,8,INERR) 
      IF(SHRMKT(3).GT.0..AND.OXYCNT(2).EQ.0. .OR. 
     *   SHRMKT(3).EQ.0..AND.OXYCNT(2).GT.0. .OR. 
     *   SHRMKT(3).EQ.0..AND.IGASHW.EQ.2) 
     *   CALL QUITER(0.,0,9,INERR) 
C 
C     The Air Toxic OXYGENATE command is at the Run, rather than at 
C     the Scenario level. Therefore, set all flags and store values 
C     accordingly. 
C 
      OXYFLG = 2 
      GLOXFLAG = 2 
C 
      GLSHARE(1) = SHRMKT(1) 
      GLSHARE(2) = SHRMKT(2) 
      GLSHARE(3) = SHRMKT(3) 
C 
      GLOXCNT(1) = OXYCNT(1) 
      GLOXCNT(2) = OXYCNT(2) 
C 
      GOTO 99 
C 
C Error or EOF on any attempted read. 
C 
   85 WRITE(IOUERR,185) M6REC 
  185 FORMAT(/'*** ERROR reading or parsing the following ', 
     &        'OXYGENATE input record:', 
     &       /'*** ',A) 
      GOTO 98 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
C 
   98 INERR = INERR + 1 
      RC = 1 
C       
   99 RETURN 
      END 
