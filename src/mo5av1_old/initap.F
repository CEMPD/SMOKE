C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE INITAP
C
C  INITAP is use to repeatedly re-initialize the AER array
C
C  Called by REINIT, and TAMPER
C
C  Input on call:
C
C    Common blocks:
C    /EGSCAL/ AER
C
C  Output on return:
C
C    Common blocks:
C    /EGSCAL/ AER
C
C
C  Notes:
C
C  This subroutine was added to MOBILE5.3 to facilitate the IM240 catalyst
C  check.  It is called to flip back and forth from ATP to IM240 cases.
C
C
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
C
      DO 40 IAY=1,2
        DO 30 IHG=1,2
          DO 20 IC2=1,11
            DO 10 IC1=1,11
              AER(IC1,IC2,IHG,IAY)=0.0
              IF(IC1.EQ.IC2) AER(IC1,IC2,IHG,IAY)=1.0
   10       CONTINUE
   20     CONTINUE
C
   30   CONTINUE
   40 CONTINUE
      RETURN
      END
