C 
C  Changes: (Last change first) 
C 
C 28 Sep 01: AIR TASK 25: Changed ATSSFLAG to SSFLAG to accomodate PM output. 
C 06 Sep 01: AIR Task 15: Changed TOXCHK back to 6, removed RVP OXY WAIVER. 
C 30 Aug 01: AIR Task 22: Added ATSSFLAG (air toxic spreadsheet) 
C 26 Jun 01: AIR Task 15: Expanded TOXCHK to 7. Added AT Waiver. 
C 10 Jul 01: AIR Task 15: Added MAXIAT 
C 26 Jun 01: AIR Task 15: Expanded TOXCHK to 6. Added ATOXY. 
C 30 Jan 01: AIR Task 11: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      INTEGER MAXIAT 
      INTEGER TOXCHK 
      LOGICAL SSFLAG 
      LOGICAL TOXFLG 
      REAL    AROM 
      REAL    ATOXY 
      REAL    BENZ 
      REAL    E200 
      REAL    E300 
      REAL    OLEF 
C 
      COMMON /ATOX1/ TOXCHK(6),E200,E300,AROM,OLEF,BENZ, 
     &               ATOXY(3,5),MAXIAT,TOXFLG,SSFLAG 
