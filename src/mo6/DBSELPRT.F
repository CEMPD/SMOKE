      SUBROUTINE DBSELPRT(INPREC,INERR,RC) 
C 
C  DBSELPRT selects the particulate elements that will appear in the database 
C  *AND* descriptive outputs. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C  30 Jan 01 AIR Task 11: New routine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELPART 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /IOUCOM/ IOUERR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  NFOUND     I   Number of valid pollutant names found on the input record. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  TESTFIELD  C   Upper case version of the pollutant type input field. 
C 
C  Notes: 
C 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'PART1.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        NFOUND 
      INTEGER        PTR1 
      INTEGER        PTR2 
      CHARACTER*10   TESTFIELD 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set all the particulate types to 'not selected', using array init. 
C 
      DBSELPART = 1 
C 
      PARTFLG=.FALSE. 
C 
C     Scan the data field of the input record, looking for the names of 
C     the pollutants that should be included in the database output. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
      PTR2 = 21 
      NFOUND = 0 
C 
   10 CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
C 
      IF (DATAFIELD.EQ.' ') GOTO 70 
C 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
C 
C     Make allowance for O (letter) or 0 (number) in SO4, SOF, SO2 
C 
      IF (TESTFIELD.EQ.'SO4' .OR. 
     *    TESTFIELD.EQ.'S04') THEN 
        DBSELPART(1) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'GASPM') THEN 
        DBSELPART(4) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'OCARBON') THEN 
        DBSELPART(2) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'ECARBON') THEN 
        DBSELPART(3) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'LEAD') THEN 
        DBSELPART(5) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'SO2' .OR. 
     *        TESTFIELD.EQ.'S02') THEN 
        DBSELPART(6) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'BRAKE') THEN 
        DBSELPART(8) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'TIRE') THEN 
        DBSELPART(9) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'NH3') THEN 
        DBSELPART(7) = 2 
        NFOUND = NFOUND + 1 
      ELSE 
        WRITE(IOUOUT,110) TESTFIELD 
        WRITE(IOUERR,110) TESTFIELD 
  110   FORMAT('*** ERROR: ',A10,' is not a valid particulate type.') 
        INERR = INERR + 1 
        RC = 1 
      ENDIF 
C 
      GOTO 10 
C 
C     If no particulate elements are listed, then ALL are selected, using 
C     array init. 
C  
   70 IF (NFOUND.EQ.0 .AND. RC.EQ.0) THEN 
        DBSELPART = 2 
        NFOUND = 9 
      ENDIF 
C 
      IF (NFOUND.EQ.0) THEN 
        WRITE(IOUOUT,140) 
        WRITE(IOUERR,140) 
  140   FORMAT(/'*** ERROR: No particulate types were selected by the', 
     &         ' PARTICULATES record.') 
        INERR = INERR + 1 
        RC = 2 
      ELSE 
        PARTFLG=.TRUE. 
      ENDIF 
C 
   99 RETURN 
      END 
