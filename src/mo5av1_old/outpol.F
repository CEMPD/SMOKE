C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTPOL(IP1,IP2)
C
C  OUTPOL selects loop parameters for pollutant printouts.
C
C  Called by OUTDT1, OUTDT2, OUTDT3 and OUTDT4.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS4/ PRTFLG
C    /MAXIMA/ MAXPOL
C
C  Output on return:
C
C    parameter list: IP1,IP2
C
C  Notes:
C
C  None.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
C
      IP1=PRTFLG
      IP2=PRTFLG
      IF(PRTFLG.LT.4) GOTO 99
      IP1=1
      IP2=MAXPOL
C
   99 RETURN
      END
