C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION GSFIDX(ICY,MY,IV,IREJN)
C
C  GSFIDX returns the sales fraction of gas vehicles based on the
C  given MY for the given IV and IREJN, using user-input diesel sales
C  fractions, if appropriate (if DSFLAG=2)
C
C  Called by EFCALX, HCCALX and REGMOD.
C
C  Calls IGSFPT.
C
C  Input on call:
C
C    parameter list: ICY,MY,IV,IREJN
C    common blocks:
C    /FLAGS1/ DSFLAG
C    /GSFCOM/ MAXGSF,GSFRAC
C    /MAXIMA/ MAXYRS
C    /REGISF/ USRGSF,DSFLAG
C
C  Output on return:
C
C    function: GSFIDX
C
C  Notes:
C
C  GSFIDX was added in MOBILE4.1 to support user-entered gas/diesel
C  sales fractions.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'GSFCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'REGISF.I'
C
C  Check DSFLAG and IV and branch if necessary
C
      IF(DSFLAG.NE.2) GO TO 10
      IF(IV.EQ.4.OR.IV.GE.7) GO TO 10
C
C  Handle user-entered gas/diesel sales fractions case
C
      IVTEMP=2
      IF(IV.EQ.1.OR.IV.EQ.5) IVTEMP=1
      JDX=ICY-MY+1
      GSFIDX=USRGSF(JDX,IVTEMP)
      IF(IV.GT.3) GSFIDX=1.0-GSFIDX
      RETURN
C
C  Handle "built-in" gas/diesel sales fraction case
C
   10 GSFIDX=GSFRAC(IGSFPT(MY,IV),IV,IREJN)
      RETURN
C
      END
