      REAL FUNCTION PMFRACTS(MY,ICY,IV,INUM,INERR) 
C 
C  PMFRACTS calculates fractions for use in any number of functions 
C 
C Changes: 
C 
C 13 Feb 01 AIR Task 11: Changed code to use MOBILE6-based vehicle 
C           classes and technology fractions. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by CEFPM, CXEFGF. 
C 
C  Calls LEADED, PMLEDPTR, PMTAMPER, CALFRACA.  
C 
C  Input on call: 
C 
C    parameter list: MY, ICY, IV, INUM 
C    INUM = 1 leaded 
C           2 unleaded catalyst 
C           3 unleaded no catalyst 
C           4 unleaded catalyst no air 
C           5 unleaded catalyst air  
C 
C 
C    common blocks: 
C 
C  Output on return: PMFRACTS 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  IYPT     I    Discretionary fuel switching index variable  
C  VLNOCT   R    Number of non-catalyst vehicles 
C                using leaded fuel 
C  VUNOCT   R    Number of non-catalyst vehicles 
C                using unleaded fuel 
C  VLYSCT   R    Number of catalyst equipped 
C                vehicles using leaded fuel 
C  VUYSCT   R    Number of catalyst equipped 
C                vehicles using unleaded fuel 
C  FRACA    R    Catalyst fractions: 1: Ox No Air, 2: Ox Air 
C                                    3: 3w No Air, 4: 3w Air 
C 
C  Notes:  
C  Mar-4-1994 @ CSC-bsg Request 340 this was added to overhaul 
C  lead emission factor calculation 
C  May-12-1994 @ CSC-bsg Request 430 this was added to help  
C  calculate various fractions in the model 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DFSBLK.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LEAD.I' 
      INCLUDE 'PART2.I' 
      INCLUDE 'PSIZE.I' 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
      INTEGER INUM 
      INTEGER IV 
      INTEGER MY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IVPART 
      INTEGER IYPT 
      REAL    CATFCT 
      REAL    CATS 
      REAL    VLNOCT 
      REAL    VLYSCT 
      REAL    VUNOCT 
      REAL    VUYSCT 
C 
C 
      PMFRACTS = 0.0 
C 
C  skip Diesels. 
C 
      IF(VVGASMC(IV).EQ.0) RETURN 
C 
C  Map MOBILE6 vehicle class into PART5 space 
C 
      IVPART=VM6M5B(IV) 
C 
C  Get the catalyst fraction of the gasoline vehicle population. 
C  (Use MOBILE6 computed values instead of PART5's FRACA.) 
C 
      CALL PMFRACA(MY,IV,INERR) 
       
      CATFCT = FRACA(1) + FRACA(2) + FRACA(3) + FRACA(4) 
C 
C  determine discretionary fuel switching index variable. 
C 
      CALL PMLEDPTR(IYPT,IV,MY) 
C 
C  determine the tampering and misfueling rates. 
C 
      CALL PMTAMPER(ICY,MY,IV) 
C 
C  determine the fraction of all vehicles which  
C  have NOT had the catalyst removed. 
C 
      CATS = CATFCT * (1.-TAMFRC) 
C 
C  determine the fraction of all vehicles which are  
C  non-catalyst vehicles fueled with leaded gasoline. 
C 
      VLNOCT = (1.-CATS) * (1.-DFS(IYPT,IVPART)) 
C 
C  determine the fraction of all vehicles which are  
C  non-catalyst vehicles fueled with unleaded gasoline. 
C 
      VUNOCT = (1.-CATS) * DFS(IYPT,IVPART) 
C 
C  determine the fraction of all vehicles which are  
C  catalyst equipped vehicles fueled with leaded gasoline. 
C 
      VLYSCT = CATFCT * (RMIS - TAMFRC) 
C 
C  determine the fraction of all vehicles which are  
C  catalyst equipped vehicles fueled with unleaded gasoline. 
C 
      VUYSCT = (1.-VLNOCT-VUNOCT-VLYSCT) 
C 
C  Set fraction depending on INUM requested 
C 
C  inum = 1: add leaded cat to leaded no cat 
C 
      IF(INUM.EQ.1) THEN 
        PMFRACTS=VLNOCT+VLYSCT 
C 
C  inum = 2: use unleaded catalyst fraction 
C 
      ELSEIF(INUM.EQ.2) THEN 
        PMFRACTS=VUYSCT 
C 
C  inum = 3: use unleaded no catalyst fraction 
C 
      ELSEIF(INUM.EQ.3) THEN 
        PMFRACTS=VUNOCT 
C 
C  inum = 4: use sum of ox no air and 3-way no air 
C 
      ELSEIF(INUM.EQ.4) THEN 
        PMFRACTS=FRACA(1)+FRACA(3) 
C 
C  inum = 5: use sum of ox air and 3-way air 
C 
      ELSEIF(INUM.EQ.5) THEN 
        PMFRACTS=FRACA(2)+FRACA(4) 
      ENDIF 
C 
      RETURN 
      END 
