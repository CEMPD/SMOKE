      SUBROUTINE DBSELTX(INPREC,INERR,RC) 
C 
C  DBSELTOX selects the air toxics elements that will appear in the database 
C  *AND* descriptive outputs. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C  12 Jul 11 AIR Task 21: Added Acrolein 
C  10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C  30 Jan 01 AIR Task 11: New routine. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELTOX 
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
      USE DATABASE, ONLY : DBSELTOX 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'IOUCOM.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      CHARACTER*10   TESTFIELD 
      INTEGER        I 
      INTEGER        NFOUND 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set the Air Toxic flag to false so that errors abort. 
C 
      TOXFLG = .FALSE. 
C 
C     Set all the air toxic types (exhaust and evap) to 'not selected.'  
C 
      DO I = 1,MAXIAT 
        DBSELTOX(I) = 1 
      END DO 
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
C     Make allowance for O (letter) or 0 (number) in COH2 and COHX 
C 
      IF     (TESTFIELD.EQ.'C6H6'  .OR. 
     *        TESTFIELD.EQ.'BENZ') THEN 
        DBSELTOX(1) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'MTBE') THEN 
        DBSELTOX(2) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'BUTA') THEN 
        DBSELTOX(3) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'COH2' .OR. 
     *        TESTFIELD.EQ.'C0H2' .OR. 
     *        TESTFIELD.EQ.'FORM') THEN 
        DBSELTOX(4) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'COHX' .OR. 
     *        TESTFIELD.EQ.'C0HX' .OR. 
     *        TESTFIELD.EQ.'ACET') THEN 
        DBSELTOX(5) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'ACRO' .OR. 
     *        TESTFIELD.EQ.'ACROLEIN') THEN 
        DBSELTOX(6) = 2 
        NFOUND = NFOUND + 1 
      ELSE 
        WRITE(IOUOUT,110) TESTFIELD 
        WRITE(IOUERR,110) TESTFIELD 
  110   FORMAT('*** ERROR: ',A10,' is not a valid Air Toxic type.') 
        INERR = INERR + 1 
        RC = 1 
      ENDIF 
C 
      GOTO 10 
C 
C     If no particulate elements are listed, then ALL are selected. 
C  
   70 IF (NFOUND.EQ.0) THEN 
      DO I = 1,MAXIAT 
        DBSELTOX(I) = 2 
      END DO 
      NFOUND = MAXIAT 
      ENDIF 
C 
      IF (NFOUND.EQ.0) THEN 
        WRITE(IOUOUT,140) 
        WRITE(IOUERR,140) 
  140   FORMAT(/'*** ERROR: No air toxics were selected by the', 
     &         ' AIR TOXICS record.') 
        INERR = INERR + 1 
        RC = 2 
      ELSE 
        TOXFLG = .TRUE. 
      ENDIF 
C 
   99 RETURN 
      END 
