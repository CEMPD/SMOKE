      SUBROUTINE ATEMIRAT(MY) 
C 
C  ATEMIRAT computes / assigns emission impact rates for the tampering effects 
C  groups. 
C 
C  Called by ATTAMPER. 
C 
C  Calls ITAMPT and P3WPOX. 
C 
C  Changes: (Last change first) 
C 
C  26 Feb 02: New MOBILE6.3 routine, based on 04 Dec 00 version of EMIRAT 
C 
C  Input on call: 
C 
C    common blocks: 
C    /EGSCAL/ TGS 
C    /EMECAL/ EM3W,EMOX,EMEGR 
C    /LOOKUP/ IVTAM 
C    /MAXIMA/ MAXPOL 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    common blocks: 
C    /EMECAL/ EMI 
C    /LOOKUP/ IQG 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  POX      R    percentage of vehicles equipped with oxidation catalysts 
C  POXN     R    POX normalized to ox cat equipped / all cat equipped 
C  P3W      R    percentage of vehicles equipped with 3-way catalysts 
C  P3WN     R    P3W normalized to 3-way cat equipped / all cat equipped 
C  TK       R    percentage of fleet equipped with catalysts 
C 
 
      IMPLICIT NONE 
C 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'EMECAL.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MAXIMA.I' 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IV, IP, IB, IEK, IG2, IEG, ITAMPT 
      REAL P3W, POX, TK, P3WN, POXN 
C 
C 
C  Compute the percentage of vehicles equipped with 3-way and oxidation 
C  catalysts. 
C 
C     NOTE: P3WPOX is also called from BIGRVP and FUEL, each time with a  
C     slightly different flavor of IV. P3WPOX has been changed to expect  
C     IV, and IVTAM is mapped back to IV here before P3WPOX is called. 
      IV = IVTAM 
      IF (IVTAM.EQ.14) THEN 
        IV = VTGASBUS 
      ENDIF 
      CALL P3WPOX(MY,IV,P3W,POX) 
C 
C  Normalize the percentages in order to use them as weights. 
C 
      TK=P3W+POX 
C 
C  Avoid divide by zero runtime error. 
C 
      IF(TK.EQ.0.0) TK=1.0 
      P3WN=P3W/TK 
      POXN=POX/TK 
C 
C  Compute / assign emissions impact for each bag for HC (air toxic) 
C 
      IP=1 ! Compute HC values used by Air Toxics 
      DO 40 IB=1,3 
C 
C  Air pump disablement (air pump only equipped) emission impact is the 
C  oxidation catalyst emission for the given bag and pollutant.  These 
C  numbers are constant with respect to model year and vehicle type and 
C  therefore are hardcoded in BLOCK DATA into EMI(1,IB,IP).  Note that 
C  EMI(1,IB,IP) = EMOX(1,IB,IP) for all IB,IP.  If EMI for air pump 
C  disablement becomes variable in MOBILE4+, the code will be inserted here. 
C 
C  Next 3 effects groups use the catalyst type weighted emissions. 
C 
      DO 30 IEG=2,4 
      IEK=IEG-1 
      EMI(IEG,IB,IP)=P3WN*EM3W(IEK,IB,IP)+POXN*EMOX(IEK,IB,IP) 
   30 CONTINUE 
C 
C  Get EGR only MY group pointer and look up corresponding emission impact. 
C 
      IG2=ITAMPT(MY,2) 
      EMI(5,IB,IP)=EMEGR(IG2,IB,IP) 
C 
C  EGR disabled / catalyst removed and EGR disabled / misfueled emissions are 
C  a constant function of bag and pollutant.  As with group 1, the values are 
C  hardcoded in BLOCK DATA into EMI(6,IB,IP) and EMI(7,IB,IP), respectively. 
C  Future versions of the model may vary these emissions on some other 
C  parameter(s) such as model year and vehicle type, in which case the 
C  computation / assignment of values would go here. 
C 
   40 CONTINUE 
C 
      RETURN 
      END 
