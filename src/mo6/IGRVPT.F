      FUNCTION IGRVPT(MY,IVTAM) 
C 
C  IGRVPT identifies the RVP CF MYG case. 
C 
C  Called by RVPEXH. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, declare MY with INTENT(IN) 
C  18 Nov 98 @DynTel 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IVTAM 
C    common blocks: 
C    /RVPEX2/ MYGRVP 
C 
C  Output on return: 
C 
C    function: IGRVPT 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
      INCLUDE 'RVPEX2.I' 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IG, IVTAM, IGRVPT 
C 
      DO 10 IG=1,4 
      IF(MY.LE.MYGRVP(IG,IVTAM)) GOTO 20 
   10 CONTINUE 
C 
C  There is no "default" index value - MYGRVP is assumed to include for 
C  for each vehicle class an upper bound year = highest allowed as an 
C  evaluation year by MOBILE = 2051. 
C 
C  At the other extreme, there are no stored terms for the first model year 
C  group, since the correction factor for that group is a constant neutral 
C  1.0.  Decrement IG so that IGRVPT correctly indexes the stored terms. 
C 
   20 IGRVPT=IG-1 
C 
      RETURN 
      END 
