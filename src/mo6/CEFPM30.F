      REAL FUNCTION CEFPM30(IV,MY,ICY,INERR) 
C 
C  Calculates EFPM for gasoline vehicles, by adding the weighted emission 
C  factor for carbon in leaded and unleaded fuel. Values are at PM30. 
C 
C Changes: 
C 
C 16 Oct 01 AIR Task 25: New MOBILE6.2 
C 18 Apr 01 AIR Task 16: New MOBILE61. 
C 13 Feb 01 AIR Task 11: Changed code to use MOBILE6 vehicle classes. 
C 17 Nov 00 AIR Task 06: Removed Dead Code. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by: EFCALXP5 
C 
C  Calls: PMFRACTS, CXEFGF30 
C 
C  Output on Return: total PM emission factor 
C                     
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  CXNOLD    R    Carbon emission factor for not leaded fuel use 
C  CXYSLD    R    Carbon emission factor for leaded fuel use 
C 
C 
C  Notes: 
C  May-10-1994 @ CSC-bsg Request 430 this was created to overhaul 
C  total exhaust particulate emission factors 
C  Mar-24-1994 @ CSC-bsg Request 415 this was created to overhaul 
C  total exhaust particulate emission factors 
C  Jun-17-1994 @ CSC-bsg Request 430 MC use lead instead of gas carbon 
C  Jun-24-1994 @ CSC-bsg Request 430 MC use lead and carbon 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I' 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PMFRACTS 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
      INTEGER IV 
      INTEGER MY 
C 
C  Declare local variables/arrays. 
C 
      REAL    CXNOLD 
      REAL    CXYSLD 
C 
C 
      CALL CXEFGF30(IV,MY,ICY,CXNOLD,CXYSLD,INERR) 
C 
      IF(IV.EQ.VTMC) THEN 
        CEFPM30 = CXNOLD + CXYSLD  
      ELSE 
        CEFPM30 = CXYSLD*PMFRACTS(MY,ICY,IV,1,INERR) + CXNOLD  
      ENDIF 
C 
      RETURN 
      END 
