C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE SCFSUB(SPD,IV)
C
C  SCFSUB returns the SCF for lower model year vehicles: which are
C  the result of an exponential 5th degree polynomial equation.
C  SCFSUB calculates "high speed" (48 MPH or greater) for LDV's and MC's
C  The SCF results for NOx are not exponentiated.
C
C  Called by BIGIN1 and IDLCAL.
C
C  Input on call:
C  parameter list: SPD,IV
C
C    common blocks:
C    /MAXIMA/ MAXPOL
C    /SPEED1/ SK1C
C    /SPEED4/ MAXSP1
C
C  Output on return:
C
C    common blocks:
C    /SPEED3/ SCUNA1
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ICOEF5   I    pointer to term in 5th degree polynomial scf equation
C  IP       I    indicates pollutant type 1=HC, 2=CO, 3=NOx
C  ISP      I    indicates speed CF MY group (IG)
C  SCF      R    speed correction factor for the given IV, IP and IG
C  TMPSPD   R    tempory average speed for the vehicle type IV to keep SCF
C                from being corrected for speeds greater than 48 mph
C
C
C  Notes:
C
C  Although calculated for all 8 IV, results are only used by BIGIN1 for
C  subgroup IVA (LDGV, LDGT1/2, MC).  Added no implicit variable option to
C  trap spelling errors.  Used standard line number and indention conventions.
C  Placed all local variables in variable dictionary.  SCFSUB was created from
C  BIGIN1 and was implemented to modularize the calculation since it is used
C  in IDLCAL to calculate idle emissions differently then before MOBILE5a.
C  Note that MOBILE5 and MOBILE5a did not calculate idle EF's.
C  July-14-1993 @ ARC-bsg Subtask 277
C
C
CC    IMPLICIT NONE
C
      INCLUDE 'MAXIMA.I'
      INCLUDE 'SPEED1.I'
      INCLUDE 'SPEED3.I'
      INCLUDE 'SPEED4.I'
C
      INTEGER IV,IP,ISP,ICOEF5
      REAL TMPSPD,SCF,SPD
C
C  Build the 5th degree polynomial ("old") equation unadjusted scf table.
C
      TMPSPD=SPD
C
C  This is not possible when called by IDLCAL
C
      IF((IV.LE.3.OR.IV.EQ.8).AND.TMPSPD.GT.48.0) TMPSPD=48.0
C
      DO 60 IP=1,MAXPOL
C
        DO 40 ISP=1,MAXSP1
          SCF=SK1(6,ISP,IP)
          DO 20 ICOEF5=1,5
            SCF=SCF*TMPSPD+SK1(6-ICOEF5,ISP,IP)
  20      CONTINUE
C
          IF(IP.NE.3) THEN
           SCF=EXP(SCF)
          ENDIF
          SCUNA1(ISP,IP,IV)=SCF
  40    CONTINUE
  60  CONTINUE
C
      RETURN
      END
