C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION IGRVPT(MY,IVTAM)
C
C  IGRVPT identifies the RVP CF MYG case.
C
C  Called by RVPEXH.
C
C  Input on call:
C
C    parameter list: MY,IVTAM
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
C
      COMMON /RVPEX2/ MYGRVP(4,4),RVPCF(25,3,4)
C
      DO 10 IG=1,4
      IF(MY.LE.MYGRVP(IG,IVTAM)) GOTO 20
   10 CONTINUE
C
C  There is no "default" index value - MYGRVP is assumed to include for
C  for each vehicle class an upper bound year = highest allowed as an
C  evaluation year by MOBILE4 = 2020.
C
C  At the other extreme, there are no stored terms for the first model year
C  group, since the correction factor for that group is a constant neutral
C  1.0.  Decrement IG so that IGRVPT correctly indexes the stored terms.
C
   20 IGRVPT=IG-1
C
      RETURN
      END
