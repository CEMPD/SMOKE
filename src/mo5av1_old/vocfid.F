C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION VOCFID(MY,IV)
C
C  VOCFID calculates the VOC FID correction factor.
C
C  Called by BEF.
C
C  Input on call:
C
C    parameter list: MY,IV
C    common blocks:
C    /FID/ MAXFID,MYFID,VOCCF
C    /OXY1/ SHRMKT
C    /OXY2/ NFUEL
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
C
      COMMON /FID/ MAXFID,MYFID(16,5),TOGCF(16,5,3),VOCCF(16,5,3)
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /OXY2/ NFUEL,MXOXMY,MAXOXY,OXYVAL(25,18,4),OXYBF(25,18,4)
C
      DIMENSION VOCEX(3)
C
      DATA VOCEX / 0.0, 0.0159, 0.0133 /
C
C  Initialize
C
      VOCFID=0.0
C
C  Map the vehicle class into FID groups
C
      IVFID=1
      IF(IV.GE.2.AND.IV.LE.3) IVFID=2
      IF(IV.EQ.4) IVFID=3
      IF(IV.GE.5.AND.IV.LE.6) IVFID=4
      IF(IV.EQ.7) IVFID=5
C
C  For motorcycles, use the pre-1975 LDGV factor
C
      IF(IV.EQ.8) IGFID=1
      IF(IV.EQ.8) GOTO 20
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
       IF(IFUEL.GT.1.AND.(IV.LE.4.OR.IV.EQ.8))
     *   PCTOXY = (OXYCNT(IFUEL-1)*100)
C
       IF(SHRMKT(IFUEL).GT.0.0)
     *   VOCFID = VOCFID + SHRMKT(IFUEL) *
     *           (VOCCF(IGFID,IVFID,IFUEL) + VOCEX(IFUEL) * PCTOXY)
   30 CONTINUE
C
      RETURN
      END
