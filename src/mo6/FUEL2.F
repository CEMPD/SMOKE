      SUBROUTINE FUEL2(KEY,BEFADJ) 
C 
C  FUEL2 computes the oxy fuel adjustments conversions from NMHC 
C  to NMOG and VOC 
C 
C  Called by HCCALX and EVP_HRLY 
C 
C  Changes: (Last change first) 
C 
C  16 May 00 @EPA- BG    bug#88 changed BEFADJ from INTEGER to REAL 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  10 Nov 99 @DynTel-BG 1-006  added index parameter constants, IDIURNAL, IHOTSOAK, 
C                              AND IRUNLOSS 
C  06 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: KEY,BEFADJ 
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
C  BEFADJ   R    adjusted emission 
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
      IMPLICIT NONE  
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'RFORM1.I' 
C       
      INTEGER, INTENT(IN):: KEY 
      REAL, INTENT(INOUT):: BEFADJ 
      INTEGER, PARAMETER :: IHOTSOAK=1,IDIURNAL=2,IRUNLOSS=3 
C 
      REAL, DIMENSION(3) :: ALCEVP=(/0.0,0.0,0.0/) 
      REAL, DIMENSION(2) :: ALCEXH=(/0.0074, 0.0159/) 
      REAL, DIMENSION(3) :: ETHEVP=(/0.0,0.0,0.0/) 
      REAL, DIMENSION(2) :: ETHEXH=(/0.0062, 0.0133/) 
      REAL, DIMENSION(2) :: PCTOXY 
C 
      IF(OXYFLG.EQ.1.AND.(.NOT.RFGON)) RETURN 
      PCTOXY(1)=OXYCNT(1)*100.0 
      PCTOXY(2)=OXYCNT(2)*100.0 
C 
C  Diurnal evaporative conversion 
C 
      IF(KEY.EQ.IDIURNAL) THEN 
        BEFADJ=SHRMKT(1)*BEFADJ+ 
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEVP(1))+ 
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEVP(1)) 
       END IF 
C 
C  Hot Soak evaporative conversion 
C 
      IF(KEY.EQ.IHOTSOAK) 
     *   BEFADJ=SHRMKT(1)*BEFADJ+ 
     *       SHRMKT(2)*BEFADJ*(1+PCTOXY(1)*ETHEVP(2))+ 
     *       SHRMKT(3)*BEFADJ*(1+PCTOXY(2)*ALCEVP(2)) 
C 
C  Running Losses evaporative conversion 
C 
      IF(KEY.EQ.IRUNLOSS) 
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
