      SUBROUTINE DBSELP(INPREC,INERR,RC) 
C 
C  DBSELP selects the pollutants that will appear in the database 
C  *AND* descriptive outputs. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 08 Mar 01 AIR Task 11: Change to allow blank input (no HC,CO,NOx output). 
C 02 Oct 00 @EPA-djb Bug230 Changed write to * to IOUOUT. 
C 13 Sep 00 @EPA-elg bug228 changed error handling conditions. 
C 10 Apr 00 AIR Task 02: Changed processing to POLLUTANTS so that 
C           database and descriptive outputs can use this routine. 
C 10 Feb 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELPOL 
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
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        I 
      INTEGER        PTR1 
      INTEGER        PTR2 
      CHARACTER*10   TESTFIELD 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set all the pollutant types to 'not selected.' They are all 
C     initialized to 'selected' so they will be printed by default 
C     unless this routine is called to restrict the output. 
C 
      DO I = 1,3 
        DBSELPOL(I) = 1 
      END DO 
C 
C     Scan the data field of the input record, looking for the names of 
C     the pollutants that should be included in the database output. 
C     Initialize by setting PTR2 to 21 (the position of the blank following 
C     the delimiter). PTR1 is set by the NXTTOK subroutine. 
C 
      PTR2 = 21 
   10 CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF (DATAFIELD.EQ.' ') GOTO 99 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
      IF (TESTFIELD.EQ.'HC') THEN 
        DBSELPOL(1) = 2 
      ELSEIF (TESTFIELD.EQ.'CO') THEN 
        DBSELPOL(2) = 2 
      ELSEIF (TESTFIELD.EQ.'NOX') THEN 
        DBSELPOL(3) = 2 
      ELSE 
        WRITE(IOUOUT,110) TESTFIELD 
        WRITE(IOUERR,110) TESTFIELD 
  110   FORMAT('*** ERROR: ',A10,' is not a valid pollutant type.') 
        INERR = INERR + 1 
        RC = 1 
      ENDIF 
C 
      GOTO 10 
C 
   99 RETURN 
      END 
