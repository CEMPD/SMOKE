      INTEGER FUNCTION IPBFUN(ICY) 
C 
C  IPBFUN selects pointer index(IPTR) for Lead content of gas by calendar 
C  year(PB). 
C 
C Changes: 
C 
C 14 Nov 00 AIR Task 06: Changed DO loop maximum to variable NIYRPB. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by  EFCALXP5 
C 
C  Calls none. 
C 
C  Input on call: 
C 
C    parameter list: ICY 
C    common block: /LEDCOM/ IYRPB, NIYRPB 
C 
C  Output on return: 
C    function: IPBFUN 
C 
C  Local variable dictionary: 
C 
C  Name    Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  IPTR     I    Pointer index for Lead content of gas by calendar year 
C 
C  Notes: none. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'LEAD.I' 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IPTR 
C 
C 
      DO 10 IPTR=1,NIYRPB 
        IF(ICY .LE. IYRPB(IPTR)) GOTO 90 
   10 CONTINUE 
C 
   90 IPBFUN=IPTR 
      RETURN 
      END 
