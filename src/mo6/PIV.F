      REAL FUNCTION PIV(MY,IVX,IPHASE) 
C 
C  PIV sets the Phase-In Values for the Cold CO Standards 
C  for LDGV/T. 
C 
C  Called by MAINTCF 
C 
C  Calls: None. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, declare MY with INTENT(IN) 
C  30 Mar 99  @DynTel-DDJ 2-696 Now called by MAINTCF for hourly emission factors 
C  18 Dec 98 @DynTel 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IVX 
C 
C    common blocks: 
C    /TEMPC5/ MYCOL 
C    /TEMPC6/ PIVPCT 
C    /VVNAMS/ VVLDGAS  (IVTYPE.I) 
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
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'TEMPC5.I' 
      INCLUDE 'TEMPC6.I' 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IMY, IVX, IPHASE 
C 
C  Initialize the Phase-In Value to zero. 
C 
      IF(IPHASE.EQ.2) THEN 
C 
C     MYCOL for phase 2 is year of full implementation, 
C     MYCOL for phase 1 is first implementation year. 
C 
         IMY = MYCOL(IPHASE) - 2 
      ELSE 
         IMY = MYCOL(IPHASE) 
      END IF 
      IF (MY.GE.IMY) THEN 
C 
C  If before MYCOL or not LDGV/T, no Phase-In occurs. 
C 
        IF (MY.EQ.IMY) THEN 
            PIV = PIVPCT(1,IVX,IPHASE)/100. 
        ELSE IF (MY.EQ.(IMY+1)) THEN 
            PIV = PIVPCT(2,IVX,IPHASE)/100. 
        ELSE  IF (MY.EQ.(IMY+2)) THEN 
            PIV = PIVPCT(3,IVX,IPHASE)/100. 
        ELSE 
            PIV = 1.0 
        END IF 
      ELSE 
        PIV = 0.0 
      END IF 
C 
      RETURN 
      END 
