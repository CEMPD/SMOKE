C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION PIV(MY,IVTEMP)
C
C  PIV sets the Phase-In Values for the Cold CO Standards
C  for LDGV/T.
C
C  Called by BIGTCF
C
C  Calls: None.
C
C  Input on call:
C
C    parameter list: MY,IVTEMP
C    common blocks:
C    /TEMPC5/ MYCOL
C    /TEMPC6/ PIVPCT
C    /TEMPC7/ ICOLD
C
C  Output on return:
C
C    function: PIV
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IPIV     I    Phase-In model year pointer MYCOL=1, MYCOL+1=2, etc.
C
C  Notes:
C
C  The Phase-In Value gives the fraction of a vehicle's model year
C  fleet affected by the standard.  For before MYCOL, the value is
C  zero. For MYCOL, it's 40%, and for MYCOL+1, it's 80%, From
C  MYCOL+2 and up, 100% are covered under the COLD CO standard.
C  Although the CCA sets only one standard, two are still allowed
C  for by the code.
C
C  PIV was added for MOBILE4.1.
C
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(3,2)
      COMMON /TEMPC6/ CSUB(3,2),CMUL(3,2),PIVPCT(3,3,2)
      COMMON /TEMPC7/ ICOLD,ADDCSW,MULCSW
C
C  Initialize the Phase-In Value to zero.
C
      PIV=0.0
C
C  If before MYCOL or not LDGV/T, no Phase-In occurs.
C
      IF(MY.LT.MYCOL(ICOLD).OR.IVTEMP.GE.4) RETURN
C
C  Set the model year pointer.
C
      IPIV=MY-MYCOL(ICOLD)+1
      IF(IPIV.GT.3) IPIV=3
C
C  Assign the Phase-In Value as a fraction.
C
      PIV=PIVPCT(IPIV,IVTEMP,ICOLD)/100.
C
      RETURN
      END
