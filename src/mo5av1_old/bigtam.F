C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGTAM(IP,IV,IDX)
C
C  BIGTAM takes the additive tampering offset from TAMPER and corrects it for
C  operating mode and temperature.
C
C  Called by BIGCFX.
C
C  Calls none.
C
C  Input on call:
C
C    parameter list: IP,IV,IDX
C    common blocks:
C    /TAMOUT/ TAMBAG
C    /TEMPC3/ CHFRAC,TCF
C
C  Output on return:
C
C    common blocks:
C    /OMTCOM/ OMTTAM, OMTCFF
C
C  Local array subscripts:
C
C  CHFHDG(3)       -  CHFHDG ( IB )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CHFHDG   R    HDGV operating mode fractions (cold, hot stabilized and
C                hot start fractions of HDGV vmt, hardcoded)
C  OMT      R    unweighted product of operating mode vmt fraction times
C                temperature correction factor
C  TOMT     R    HDGV tampering offset adjusted for operating mode
C
C  Notes:
C
C  None.
C
C
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
C
      DIMENSION CHFHDG(3)
C
      DATA CHFHDG/.206,.521,.273/
C
C  HDGV: a special (ignoring user inputted scenario operating mode) opmode
C  correction to the TAMPER offset, the result then corrected by the "bag 2"
C  RVP & T cf (HDGV opmode for TCF is assumed 100% hot stabilized = bag 2).
C
      IF(IV.NE.4) GOTO 20
      TOMT=0.0
      DO 10 IB=1,3
      TOMT=TOMT+CHFHDG(IB)*TAMBAG(IB,IP,IDX,IV)
   10 CONTINUE
      OMTTAM(IDX,IP,IV)=TOMT*TCF(2)
      RETURN
C
C  LDV/T: TAMPER's offset is corrected by bag for both operating mode and
C  temperature.
C
C  Note: multiplicative tcf for correcting tampering offset is neutral (1.0)
C  when LOWCO = 2 ( = cold start LDGV/T 1980/4+ CO -> use additive CO offset).
C
   20 OMTTAM(IDX,IP,IV)=0.0
      DO 30 IB=1,3
      OMT=CHFRAC(IB)*TCF(IB)
      OMTTAM(IDX,IP,IV)=OMTTAM(IDX,IP,IV)+OMT*TAMBAG(IB,IP,IDX,IV)
   30 CONTINUE
C
      RETURN
      END
