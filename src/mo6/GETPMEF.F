      SUBROUTINE GETPMEF(INPREC,INERR,RC) 
C 
C  GETPMEF obtains the gas carbon and diesel PM ZML, DR1 and DR2 filenames 
C  from the PARTICULATE EF command and then calls the routines to read in the data. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,UCCOMP,GETPMGZML,GETPMGDR1,GETPMGDR2, 
C        GETPMDZML,GETPMDDR1,GETPMDDR2 
C 
C  Changes: (Last change first) 
C 
C  20 Oct 01  AIR Task 25: Now reads Gas Carbon and Diesel Total PM separately 
C  28 Aug 01: Air Task 23c: DR1 and DR2 implementation. 
C  25 May 01: AIR Task 17: New Subroutine. 
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
C 
C  Local variable / array dictionary: 
C 
C   Name    Type                      Description 
C  ------   ----  --------------------------------------------------------- 
C  DATAFIELD  C   One blank-delimited data item from the input record.  
C  FILENAME   C   Input filename for PM ZML or DR values.  
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'PART1.I' 
C 
      CHARACTER*(*)  INPREC 
      INTEGER        INERR 
      INTEGER        RC 
C 
      CHARACTER*80   DATAFIELD 
      CHARACTER*80   FILENAME 
      INTEGER        PTR1 
      INTEGER        PTR2 
C 
C ..................................................................... 
C 
C    Initially set PRTCHK for PARTICULATE EF to abort in routine PMCALX. 
C 
      PRTCHK(3) = 1 
C 
C    Look for the first data field (gas carbon ZML filename) following the delimiter 
C    in column 20. Initialize by setting PTR2 to 21 (the position of the blank 
C    following the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETPMGZML. 
C 
      CALL GETPMGZML(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the second data field (gas carbon DR1 filename). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETPMGDR1. 
C 
      CALL GETPMGDR1(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the third data field (gas carbon DR2 filename). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETPMGDR2. 
C 
      CALL GETPMGDR2(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the fourth data field (diesel total PM ZML filename). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETPMDZML. 
C 
      CALL GETPMDZML(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the fifth data field (diesel total PM DR1 filename)). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETDPMDDR1. 
C 
      CALL GETPMDDR1(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the sixth data field (diesel total PM DR2 filename)). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)
C 
C     Process this filename via GETPMDDR2. 
C 
      CALL GETPMDDR2(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Gas Carbon and Diesel Total PM ZML, DR1 and DR2 data have been  
C     properly read in, so update PRTCHK 
C 
      PRTCHK(3) = 2 
C 
   99 RETURN 
      END 
        
