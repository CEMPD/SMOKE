      SUBROUTINE DBWRUN(INERR,RC) 
C 
C  DBWRUN writes the run-level data to the database output tables. 
C 
C  Called by CHKINP 
C 
C  Calls QUITER 
C 
C  Changes: (Last change first) 
C 
C 
C 
C  Input on call: 
C 
C  Output on return: 
C 
C    argument list: INERR, RC 
C 
C  Modules used: 
C 
C    DATABASE - DBBATCH, DBFLAG, DBNFILE, DBNRUN, DBSTRNUM, DBUNIT, 
C               TABCHAR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  I        I    Loop control variable (implicit DO loops).  
C  ITAB     I    Database table number (for error messages). 
C 
C 
C 
C  Notes: 
C 
C  This subroutine is new for Mobile6. 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'VDATA.I' 
C 
      INTEGER    INERR 
      INTEGER    RC 
C 
      INTEGER    I 
      INTEGER    ITAB 
C 
      RC = 0 
C 
C     Return if no database output is requested or if the database was 
C     not opened correctly. 
C 
      IF (DBFLAG /= 2) GOTO 90 
C 
C     Write the batch file name to Table 3, if necessary. The batch file 
C     always has zeroes for the file number, run number, and scenario 
C     number fields and a file type of 'BATCH.' 
C 
      IF (DBNFILE == 1 .AND. DBNRUN == 1 .AND. DBBATCH /= ' ') THEN 
        WRITE(DBUNIT(3),120,IOSTAT=RC,ERR=85) (0,TABCHAR,I=1,3), 
     &        'BATCH',TABCHAR,DBBATCH 
      ENDIF 
  120 FORMAT(3(I4,A1),A8,A1,A80) 
C 
C     Write the name of the current input file to the DATAFILES table, 
C     if necessary. Also save the VERSION string in the MESSAGES table. 
C     The VERSION string will be saved once for each input file. Unless 
C     we are running in batch mode, it will be saved exactly once. In 
C     batch mode, the keys are still unique because the file number is 
C     different in each VERSION record. 
C 
      IF (DBNRUN == 1) THEN 
C 
        ITAB = 3 
        WRITE(DBUNIT(3),120,IOSTAT=RC,ERR=85) DBNFILE,TABCHAR, 
     &        DBNRUN,TABCHAR,0,TABCHAR,'INPUT',TABCHAR,DBINPUT 
C 
        ITAB = 4 
        DBSTRNUM = DBSTRNUM + 1 
        WRITE(DBUNIT(4),130,IOSTAT=RC,ERR=85) DBNFILE, TABCHAR, 
     &        0,TABCHAR,0,TABCHAR,'VERSION',TABCHAR,DBSTRNUM,TABCHAR, 
     &        1,TABCHAR,VERSION 
  130   FORMAT(3(I4,A1),A8,A1,2(I4,A1),A80) 
C 
      ENDIF 
C 
C     Store the run parameters in the RUNDATA table. 
C 
      ITAB = 5 
      WRITE(DBUNIT(5),140,IOSTAT=RC,ERR=85) DBNFILE,TABCHAR, 
     *       DBNRUN,(TABCHAR,DBRUNFLAGS(I),I=1,14) 
  140 FORMAT(I4,A1,I4,A1,14(A1,A1)) 
C 
      GOTO 90 
C 
   85 CONTINUE 
      CALL QUITER(0.0, ITAB, 525, INERR) 
      GOTO 90 
C 
   90 RETURN 
      END 
