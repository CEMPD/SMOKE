      FUNCTION NMOGFID(MY,IV,IOXY) 
C 
C  NMOGFID calculates the TOG FID correction factor. 
C  Name was changed to NMOGFID from TOGFID 
C 
C  Called by BEF, BEF94P, EFNFTP, FUEL. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations 
C  13 January 1999 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IV,IOXY 
C    common blocks: 
C    /FID/ MAXFID,MYFID,TOGCF 
C    /OXY1/ SHRMKT 
C    /OXY2/ NFUEL 
C    /VVNAMS/ VVGASMC, VVHDDV, VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C    /VVTYPS/ VTLDDV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: NMOGFID 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IFUEL    I    1=Gasoline, 2=Ether, and 3=Ethanol 
C  IGFID    I    my group for FID corrections 
C  IVFID    I    vehicle group for FID corrections( 1=LDGV/MC, 2=LDGT, 
C                3=HDGV, 4=LDDV/T, 5=HDDV) 
C  PCTOXY   R    Local OXY content 
C  TOGEX    R    Effect of Oxy on the TOG CF 
C 
C  Notes: 
C 
C  NMOGFID was added for MOBILE4.1. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FID.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'OXY2.I' 
C 
      INTEGER IVFID, IV, IGFID, IOXY, IFUEL, MY 
C 
      REAL PCTOXY, NMOGFID 
      REAL, DIMENSION(3) :: TOGEX=(/ 0.0, 0.0074, 0.0062 /) 
C 
C  Initialize 
C 
      NMOGFID=0.0 
C 
C  Map the vehicle class into FID groups 
C 
      IVFID=1 
      IF(VVLDGT12(IV).EQ.1.OR.VVLDGT34(IV).EQ.1) IVFID=2 
      IF(VVHDGV(IV).EQ.1) IVFID=3 
      IF(IV.EQ.VTLDDV.OR.VVLDDT(IV).EQ.1) IVFID=4 
      IF(VVHDDV(IV).EQ.1) IVFID=5 
C 
C  For motorcycles, use the pre-1975 LDGV factor 
C 
      IF(IV.EQ.VTMC) THEN  
        IGFID=1 
        GOTO 20 
      ENDIF 
C 
C  Determine the my group 
C 
      DO 10 IGFID=1,MAXFID 
      IF(MY.LE.MYFID(IGFID,IVFID)) GOTO 20 
   10 CONTINUE 
C 
C  Determine the market share weighted TOG FID correction factor. 
C 
   20 DO 30 IFUEL=1,NFUEL 
       PCTOXY = 0.0 
       IF(IFUEL.GT.1.AND.VVGASMC(IV).EQ.1.AND.IOXY.EQ.2) 
     *   PCTOXY = OXYCNT(IFUEL-1)*100. 
       IF(SHRMKT(IFUEL).GT.0.0) 
     *   NMOGFID = NMOGFID + SHRMKT(IFUEL) * 
     *           (TOGCF(IGFID,IVFID,IFUEL) + TOGEX(IFUEL) * PCTOXY) 
   30 CONTINUE 
C 
      RETURN 
      END 
                                               
