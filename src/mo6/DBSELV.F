      SUBROUTINE DBSELV(INPREC,INERR,RC) 
C 
C  DBSELV selects the vehicles that will appear in the database output. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 11 Oct 00 AIR Task 02 Bug Fix 275: Improved read error handling. 
C 02 Oct 00 @EPA-bag Bug183  remove QUITER call for error message 
C 22 Sep 00 AIR Task 07 Bug Fix 183: Trap incorrect number of flags. 
C 10 MAR 00 @DynTel-MLA 1-013  Replaced vehicle numbers with MAXVEH. 
C 10 Feb 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELVEH 
C    /IVTYPE/    MAXVEH 
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
C 
C  Notes: 
C 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10 :: DATAFIELD 
      CHARACTER*19 :: VEHICLEGROUP 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
      INTEGER        I 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set all the vehicle types to 'not selected.' They are all 
C     initialized to 'selected' so they will be printed by default 
C     unless this routine is called to restrict the output. 
C 
      DO I = 1,MAXVEH 
        DBSELVEH(I) = 1 
      END DO 
C 
C     Read the flags for the vehicle types selected. There should be a 
C     group of five flags (for LDGV/T), a group of eight (for HDGV), a 
C     group of 2 (for LDDV/T), another group of 8 (for HDD), a single 
C     flag (for MC), and finally a group of 3 (for BUS). If the data 
C     length is not as expected, quit read and warn user. 
C 
      PTR2 = 21 
C 
      VEHICLEGROUP = 'LIGHT DUTY GAS    ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.5) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELVEH(I),I=1,5) 
C 
      VEHICLEGROUP = 'HEAVY DUTY GAS    ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELVEH(I),I=6,13) 
C 
      VEHICLEGROUP = 'MOTOROCYCLE       ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.1) GOTO 95 
      READ (DATAFIELD,110,ERR=96)  DBSELVEH(24) 
C 
      VEHICLEGROUP = 'LIGHT DUTY DIESEL ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.3) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELVEH(I),I=14,15) 
      READ (DATAFIELD(3:),110,ERR=96) DBSELVEH(28) 
C 
      VEHICLEGROUP = 'HEAVY DUTY DIESEL ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.8) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELVEH(I),I=16,23) 
C 
      VEHICLEGROUP = 'BUS               ' 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.3) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELVEH(I),I=25,27) 
  110 FORMAT(8I1) 
C 
C     Check the data. 
C 
      DO I = 1,MAXVEH 
        IF(DBSELVEH(I).LT.1.OR.DBSELVEH(I).GT.2) THEN 
          WRITE(IOUOUT,10) DBSELEFT(I),IOOFBD 
          WRITE(IOUERR,10) DBSELEFT(I),IOOFBD 
   10     FORMAT(13X,I4,5A4,'a database selection flag (1 or 2).') 
          IF(RC.EQ.0) INERR = INERR + 1 
          RC = 2 
        ENDIF 
      END DO 
      GOTO 99 
C 
C     Data length error encountered 
   95 WRITE(IOUOUT,15) VEHICLEGROUP 
      WRITE(IOUERR,15) VEHICLEGROUP 
   15 FORMAT( 
     *'*** Error:  The number of vehicle flags entered',/ 
     *13X,'for the DATABASE VEHICLES command in the ',/ 
     *13X, A19,/ 
     *13X,'vehicle group is incorrect.'/) 
      GOTO 98 
C 
C     Read errors. 
C 
   96 WRITE(IOUOUT,25) TRIM(DATAFIELD) 
      WRITE(IOUERR,25) TRIM(DATAFIELD) 
   25 FORMAT(/'*** ERROR: The following DATABASE VEHICLES input line'/ 
     &        '    cannot be read: ',A) 
C      
   98 RC = 2 
      INERR = INERR + 1 
C 
   99 RETURN 
      END 
