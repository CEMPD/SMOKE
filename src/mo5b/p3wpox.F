C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE P3WPOX(IV,P3W,POX)
C
C  Computes the number of catalyst and non-catalyst vehicles,
C
C  Called by BIGRVP, EMIRAT and FUEL
C
C  Calls ITAMPT
C
C  Input on call:
C
C  Parameter list: IV
C    common blocks:
C    /EGSCAL/ TGS
C
C  Output on call:
C
C  Parameter list: P3W, POX
C    common blocks:
C    /LOOKUP/ IVTAM,IQG
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IG1      I    Return value of ITAMPT
C  IQG      I    Pointer in TGS array
C  IVTAM    I    See parameter dictionary
C
C
C  Notes: 
C
C  Sep-22-1994 @ CSC-bsg Request 426  P3WPOX was created to modularize
C  BIGRVP, EMIRAT and FUEL.  BIGRVP is needed to calculate the NOx RFG
C  effects based on percentage of vehicles equipped with 3-way catalysts.
C
C
      INCLUDE 'EGSCAL.I'
      INCLUDE 'LOOKUP.I'
C
      P3W=0.0
      POX=0.0
C
      IF(IV.LE.4) THEN
      IVTAM=IV
      DO 10 IQG=5,6
      IG1=ITAMPT(1)
      P3W=P3W+TGS(IG1,IQG,IVTAM)
   10 CONTINUE
C
      DO 20 IQG=2,3
      IG1=ITAMPT(1)
      POX=POX+TGS(IG1,IQG,IVTAM)
   20 CONTINUE
      POX=POX-P3W
      ENDIF
C
      RETURN
      END
