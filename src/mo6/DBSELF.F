      SUBROUTINE DBSELF(INPREC,INERR,RC) 
C 
C  DBSELF selects the facility types that will appear in the database 
C  output. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C 13 Sep 00 @EPA-elg bug228 changed error handling conditions. 
C 10 Feb 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /IOUCOM/ IOUERR 
C 
C  Output on return: 
C 
C    argument list: INERR,RC 
C    common blocks: 
C    /DATABASE/  DBSELFAC 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type                      Description 
C  ------  ----  --------------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  NFOUND     I   Number of valid faclity types found on the input record. 
C  PTR1       I   Pointer used in NXTTOK to parse the input record. 
C  PTR2       I   Pointer used in NXTTOK to parse the input record. 
C  TESTFIELD  C   Upper case version of the facility type input field. 
C 
C  Notes: 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        I 
      INTEGER        NFOUND 
      INTEGER        PTR1 
      INTEGER        PTR2 
      CHARACTER*10   TESTFIELD 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set all the facility types to 'not selected.' They are all 
C     initialized to 'selected' so they will be printed by default 
C     unless this routine is called to restrict the output. 
C 
      DO I = 1,5 
        DBSELFAC(I) = 1 
      END DO 
C 
C     Scan the data field of the input record, looking for the names of 
C     the facility types that should be included in the database output. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
      PTR2 = 21 
      NFOUND = 0 
   10 CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF (DATAFIELD.EQ.' ') GOTO 70 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
      IF (TESTFIELD.EQ.'FREEWAY') THEN 
        DBSELFAC(1) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'ARTERIAL') THEN 
        DBSELFAC(2) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'LOCAL') THEN 
        DBSELFAC(3) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'RAMP') THEN 
        DBSELFAC(4) = 2 
        NFOUND = NFOUND + 1 
      ELSEIF (TESTFIELD.EQ.'NONE') THEN 
        DBSELFAC(5) = 2 
        NFOUND = NFOUND + 1 
      ELSE 
        WRITE(IOUOUT,110) TESTFIELD 
        WRITE(IOUERR,110) TESTFIELD 
  110   FORMAT(' *** ERROR: ',A10,' is not a valid facility type.') 
        INERR = INERR + 1 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
      GOTO 10 
C 
   70 IF (NFOUND.EQ.0) THEN 
        WRITE(IOUOUT,140) 
        WRITE(IOUERR,140) 
  140   FORMAT(/' *** ERROR: No facility types were selected by', 
     &         ' the DATABASE FACILITIES record.') 
        INERR = INERR + 1 
        RC = 2 
      ENDIF 
C 
   99 RETURN 
      END 
