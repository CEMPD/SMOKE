      SUBROUTINE PMLEDPTR(IYPT, IV, MY) 
C 
C  PMLEDPTR selects index pointers for Base Lead Emission Rates(DFS) and 
C  Sales Fractions for Leaded Fuel Vehicles(PBFRAC). 
C 
C Changes: 
C 
C 12 Feb 01 AIR Task 11: Explicit use of PART5 vehicle class IVPART. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by LEADED and UNLDED. 
C 
C  Calls none. 
C 
C  Input on call: 
C 
C    parameter list: IYPT, IV, MY 
C    common blocks: 
C    /DFSBLK/ MYLEAD, MAXLED 
C 
C  Output on return: 
C 
C    parameter list: IYPT 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  IVPART   I    PART5 vehicle class equivalent of MOBILE6 vehicle class  
C 
C  Notes: 5-31-94 @ CSC-bsg Request 434 Integration and regression testing 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DFSBLK.I' 
      INCLUDE 'IVTYPE.I' 
C 
C  Declare parameter list. 
C 
      INTEGER IV 
      INTEGER IYPT 
      INTEGER MY 
C 
C  Declare local variables. 
C 
      INTEGER IVPART 
C 
C  Map MOBILE6 vehicle classes into PART5 space. GasBuses are 
C  treated as HDGV. 
C 
      IVPART=VM6M5B(IV) 
C 
C  Select pointer index(IYPT) for Base Lead Emission Rates(DFS) 
C 
      DO 10 IYPT=1, MAXLED 
        IF(MY .LE. MYLEAD(IYPT,IVPART)) GOTO 90 
   10 CONTINUE 
C 
   90 RETURN 
      END 
