      SUBROUTINE GETPMEF(INPREC,INERR,RC) 
C 
C  GETPMEF obtains the PM ZML, DR1 and DR2 filenames from the PARTICULATE EF command 
C  and then calls the routines to read in the data. 
C 
C  Called by DSPTCH. 
C 
C  Calls NXTTOK,UCCOMP,GETPMZML,GETPMDR1,GETPMDR2 
C 
C  Changes: (Last change first) 
C 
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
C    Look for the first data field (ZML filename) following the delimiter 
C    in column 20. Initialize by setting PTR2 to 21 (the position of the blank 
C    following the delimiter). PTR1 does not need to be initialized. 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 

ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)

C 
C     Process this filename via GETPMZML. 
C 
      CALL GETPMZML(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the second data field (DR1 filename). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 

ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)

C 
C     Process this filename via GETPMDR1. 
C 
      CALL GETPMDR1(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     Look for the third data field (DR2 filename). 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 

ccs      CALL UCCOMP(DATAFIELD,FILENAME) 
      CALL LJCOMP(DATAFIELD,FILENAME)

C 
C     Process this filename via GETPMDR2. 
C 
      CALL GETPMDR2(TRIM(FILENAME),INERR,RC) 
C 
C     If any error occurs, abort 
C 
      IF(INERR.GT.0.OR.RC.NE.0) RETURN 
C 
C     PM ZML, DR1 and DR2 data have been properly read in, so update PRTCHK 
C 
      PRTCHK(3) = 2 
C 
   99 RETURN 
      END 
        
