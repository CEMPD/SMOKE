      SUBROUTINE DBSELE(INPREC,INERR,RC) 
C 
C  DBSELE selects the ETypes that will appear in the database output. 
C 
C  Called by DBOPTS, PROCHDR. 
C 
C  Calls LJCOMP,NXTTOK,QUITER,UCCOMP. 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 06 Mar 01 AIR Task 11: Added brake and tire wear flags. 
C 11 Oct 00 AIR Task 02 Bug Fix 275: Improved read error handling. 
C 02 Oct 00 @EPA-bag Bug183  Call QUITER to write error message 
C 22 Sep 00 AIR Task 07 Bug Fix 183: Trap incorrect number of flags. 
C 10 Feb 99 @DynTel-MLA 2-684  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: INPREC 
C    common blocks: 
C    /DATABASE/  DBSELETP 
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
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
      CHARACTER*(*)  INPREC 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*10   DATAFIELD 
      INTEGER        I 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C ..................................................................... 
C 
      RC = 0 
C 
C     Set all the ETypes to 'not selected.' They are all initialized 
C     to 'selected' so they will be printed by default unless this 
C     routine is called to restrict the output. If the data length is 
C     not as expected, quit read and warn user. 
C 
      DO I = 1,10 
        DBSELEFT(I) = 1 
      END DO 
C 
C     Read the flags for the ETypes selected. There are two groups 
C     of four flags. 
C 
      PTR2 = 21 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.4) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELEFT(I),I=1,4) 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.4) GOTO 95 
      READ (DATAFIELD,110,ERR=96) (DBSELEFT(I),I=5,8) 
C 
C     Brake and Tire Wear flags are optional, so no error message 
C     occurs if they are missing. Default is to not to report them. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      IF(LEN(TRIM(DATAFIELD)).NE.2) GOTO 5 
      READ (DATAFIELD,110,ERR=96) (DBSELEFT(I),I=9,10) 
C 
  110 FORMAT(4I1) 
C 
C     Check the data. 
C 
    5 DO I = 1,10 
        IF(DBSELEFT(I).LT.1.OR.DBSELEFT(I).GT.2) THEN 
          WRITE(IOUOUT,10) DBSELEFT(I),IOOFBD 
          WRITE(IOUERR,10) DBSELEFT(I),IOOFBD 
   10     FORMAT(13X,I4,5A4,'a database selection flag (1 or 2).') 
          IF(RC.EQ.0) INERR = INERR + 1 
          RC = 2 
        ENDIF 
      END DO 
      GOTO 99 
C 
C     Data length error encountered. 
   95 WRITE(IOUOUT,15) 
      WRITE(IOUERR,15) 
   15 FORMAT( 
     *'*** Error:  The number of emission type flags'/ 
     *13X,'entered for the DATABASE EMISSIONS command'/ 
     *13X,'is incorrect.'/) 
      GOTO 98 
C 
C     Read errors. 
C 
   96 WRITE(IOUOUT,25) TRIM(DATAFIELD) 
      WRITE(IOUERR,25) TRIM(DATAFIELD) 
   25 FORMAT(/'*** ERROR: The following DATABASE EMISSIONS input line'/ 
     &        '    cannot be read: ',A) 
C      
   98 RC = 2 
      INERR = INERR + 1 
C 
   99 RETURN 
      END 
