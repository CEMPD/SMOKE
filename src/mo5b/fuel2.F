C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE FUEL2(KEY,BEFADJ)
C
C  FUEL2 computes the oxy fuel adjustments conversions from NMHC
C  to NMOG and VOC
C
C  Called by BEF and HCCALX
C
C  Input on call:
C
C    parameter list: KEY,BEFADJ
C    common blocks:
C     /FLAGS1/ OXYFLG
C     /FLAGS4/ NMHFLG
C     /OXY1/   SHRMKT,OXYCNT
C     /RFORM1/ RFGON
C
C  Output on return:
C
C    subroutine: FUEL2
C    common blocks:
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ALCEVP   R    conversion factor for alcohol evaporative adjustment
C                (1) is for Diurnal conversion
C                (2) is for Hot Soak conversion
C                (3) is for Running Losses conversion
C  ALCEXH   R    conversion factor for alcohol exhaust adjustment
C                (1) is for conversion from NMHC to NMOG
C                (2) is for conversion from NMHC to VOC
C  ETHEVP   R    conversion factor for ether evaporative adjustment
C                (1) is for Diurnal conversion
C                (2) is for Hot Soak conversion
C                (3) is for Running Losses conversion
C  ETHEXH   R    conversion factor for ether exhaust adjustment
C                (1) is for conversion from NMHC to NMOG
C                (2) is for conversion from NMHC to VOC
C  PCTOXY   R    percent of oxygen
C
C  Notes:
C
C  FUEL2 was added for MOBILE5.6, and copied form FUEL but is used
C  as a subroutine instead of a function.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'RFORM1.I'
C
      DIMENSION ALCEVP(3),ALCEXH(2),ETHEVP(3),ETHEXH(2)
      DIMENSION PCTOXY(2)
      DATA ALCEVP / 3*0.0 /
      DATA ALCEXH / 0.0074, 0.0159 /
      DATA ETHEVP / 3*0.0 /
      DATA ETHEXH / 0.0062, 0.0133 /
C
      IF(OXYFLG.EQ.1.AND.(.NOT.RFGON)) RETURN
      PCTOXY(1)=OXYCNT(1)*100.0
      PCTOXY(2)=OXYCNT(2)*100.0
C
C  Diurnal evaporative conversion
C
      IF(KEY.EQ.1)
     *   BEFADJ=SHRMKT(1)*BEFADJ+
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEVP(1))+
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEVP(1))
C
C  Hot Soak evaporative conversion
C
      IF(KEY.EQ.2)
     *   BEFADJ=SHRMKT(1)*BEFADJ+
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEVP(2))+
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEVP(2))
C
C  Running Losses evaporative conversion
C
      IF(KEY.EQ.3)
     *   BEFADJ=SHRMKT(1)*BEFADJ+
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEVP(3))+
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEVP(3))
C
C  Apply conversion only to (VOC), (TOG) or (NMOG) exhaust
C
      IF(NMHFLG.LT.3) RETURN
      IF(KEY.LT.4) RETURN
C
C  Use the NMOG when NMHFLG is 4 or 5
C
      IF(NMHFLG.GT.3)
     *  BEFADJ=SHRMKT(1)*BEFADJ+
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEXH(1))+
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEXH(1))
C
C  Use the VOC when NMHFLG is 3
C
      IF(NMHFLG.EQ.3)
     *  BEFADJ=SHRMKT(1)*BEFADJ+
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEXH(2))+
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEXH(2))
C
      RETURN
      END
