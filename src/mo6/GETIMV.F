      SUBROUTINE GETIMV(INPREC,INERR,RC) 
C 
C  GETIMV processes the I/M VEHICLES labeled input record. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,QUITER 
C 
C  Changes: (Last change first) 
C 
C  27 Sep 00 AIR Task 07: Bug Fix 255: skip over QUITER call for bug 183 
C  03 Oct 00 @EPA-bg Task 07: Bug Fix 183, added vehicle type info. for 
C            error message 
C  22 Sep 00 AIR Task 07: Bug Fix 183, trap incorrect number of flags. 
C  25 Aug 00 AIR Task 08: Allowed for up to MAXIMPGM programs. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /CHKLST/ IMCHECK 
C    /IMPAR2/ ILDT 
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  IMPGM      I   The I/M program number which uses this data. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  X          R   Temporary real to integer variable 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'IMPAR2.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        I 
      INTEGER        IMPGM 
      INTEGER        IOS 
      INTEGER        IVIM 
      INTEGER        IVEHICLEGROUP 
      INTEGER        PTR1 
      INTEGER        PTR2 
      REAL           X 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Look for the first data field following the delimiter in column 20. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*,ERR=95,END=95,IOSTAT=IOS) X 
      IMPGM=INT(X) 
      IF (IMPGM.LT.1.OR.IMPGM.GT.MAXIMPGM) THEN 
        CALL QUITER(0.0,IMPGM,153,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Show that some data was entered for this I/M program and check off 
C     the I/M VEHICLES data item for this I/M program. 
C 
      IMCHECK(1,IMPGM) = 2 
      IMCHECK(9,IMPGM) = 2 
C 
C     Read the flags for the vehicle types covered by this I/M program. 
C     There should be a group of five flags (for LDGV/T), a group of eight 
C     flags (for HDGV), and a single flag (for GAS BUS) at the end of the 
C     record. If the data length is not as expected, quit read and warn user. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IVEHICLEGROUP=1 
      IF(LEN(TRIM(DATAFIELD)).NE.5) GOTO 98 
      READ (DATAFIELD,110,ERR=95,IOSTAT=IOS) (ILDT(I,IMPGM),I=1,5) 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IVEHICLEGROUP=2 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 98 
      READ (DATAFIELD,110,ERR=95,IOSTAT=IOS) (ILDT(I,IMPGM),I=6,13) 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IVEHICLEGROUP=3 
      IF(LEN(TRIM(DATAFIELD)).NE.1) GOTO 98 
      READ (DATAFIELD,110,ERR=95,IOSTAT=IOS) ILDT(14,IMPGM) 
  110 FORMAT(8I1) 
C 
C     Check the data. 
C 
      DO IVIM = 1,14 
        IF(ILDT(IVIM,IMPGM).LT.1.OR.ILDT(IVIM,IMPGM).GT.2) THEN 
          CALL QUITER(0.,ILDT(IVIM,IMPGM),72,INERR) 
          RC = 2 
        ENDIF 
      END DO 
C 
      GOTO 99 
C 
C Error or EOF on any attempted read. 
C 
   95 RC = IOS 
      CALL QUITER(0.0,IOS,172,INERR) 
      GOTO 99 
C 
   98 CALL QUITER(0.,IVEHICLEGROUP,129,INERR) 
      RC = 2 
C 
C 
   99 RETURN 
      END 
