      REAL FUNCTION VOCFID(MY,IV) 
C 
C  VOCFID calculates the VOC FID correction factor. 
C 
C  Called by BEF, BEF94P, EFNFTP. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initialization 
C  15 Mar 00 @DynTel-ddj, 1-032 Changed conditional to 
C            (VVLDDSL(IV) .EQ. 1) from 
C            (IV.EQ.VTLDDV.OR.IV.EQ.VTLDDT). 
C  24 Jan 00 @DynTel-MLA 1-011  Replaced VTLDDT with VVLDDT (LDDT renamed). 
C  14 January 1999 @DynTel-ZK Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IV 
C    common blocks: 
C    /FID/ MAXFID,MYFID,VOCCF 
C    /OXY1/ SHRMKT 
C    /OXY2/ NFUEL 
C    /VVNAMS/ VVGASMC, VVHDDV, VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C    /VVTYPS/ VTLDDT, VTLDDV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: VOCFID 
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
C  VOCEX    R    Effect of Oxy on the VOC CF 
C 
C  Notes: 
C 
C  VOCFID was added for MOBILE4.1. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FID.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'OXY2.I' 
C 
C     declare parameter list 
C 
      INTEGER IV 
      INTEGER MY 
C 
C     declare local variables 
C 
      INTEGER IGFID 
      INTEGER IFUEL 
      INTEGER IVFID 
      REAL    PCTOXY 
      REAL, DIMENSION(3) :: VOCEX=(/0.0, 0.0159, 0.0133/) 
C 
C  Initialize 
C 
      VOCFID=0.0 
C 
C  Map the vehicle class into FID groups 
C 
      IVFID=1 
      IF(VVLDGT12(IV).EQ.1.OR.VVLDGT34(IV).EQ.1) THEN  
        IVFID=2 
      ELSEIF(VVHDGV(IV).EQ.1) THEN 
        IVFID=3 
      ELSE IF(VVLDDSL(IV) .EQ. 1) THEN 
        IVFID=4 
      ELSEIF(VVHDDV(IV).EQ.1) THEN 
        IVFID=5 
      ENDIF 
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
C  Determine the market share weighted VOC FID correction factor. 
C 
   20 DO 30 IFUEL=1,NFUEL 
       PCTOXY = 0.0 
       IF(IFUEL.GT.1.AND.VVGASMC(IV).EQ.1) 
     *   PCTOXY = (OXYCNT(IFUEL-1)*100) 
C 
       IF(SHRMKT(IFUEL).GT.0.0) 
     *   VOCFID = VOCFID + SHRMKT(IFUEL) * 
     *           (VOCCF(IGFID,IVFID,IFUEL) + VOCEX(IFUEL) * PCTOXY) 
   30 CONTINUE 
C 
      RETURN 
      END 
