C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE NPARAM(DDISTP,ID,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE)
C
C  NPARAM is used to get the correct parameters and leaves the old ones
C  alone as a function of the DDISTP(ID).
C
C  Called by ATPAER
C
C  Input on call:
C    Parameter list: DDISTP,ID,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE
C    common blocks:
C    /ATPAR2/ ATPGM,ATPFQT
C
C  Output on return:
C    Parameter list: NAPGM,NBFQT,NCRATE
C
C
C  Notes:
C
C  This subroutine was added to MOBILE5.3 to facilitate the IM240 test
C  effects on the ATP effectiveness.  It is called to flip back and
C  forth from ATP to IM240 parameters.  The parameters switched are
C  inspection type, frequency and compliance.
C
C
      REAL NCRATE
      INTEGER ATPPGM,ATPFQT,DISTYP
      INTEGER DDISTP(8),APGM,BFQT
C
      COMMON /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP(8)
C
      NAPGM=ATPPGM
      NBFQT=ATPFQT
      NCRATE=CRATP
      IF(DDISTP(ID).EQ.3) THEN
         NAPGM=APGM
         NBFQT=BFQT
         NCRATE=CRATE
      ENDIF
      RETURN
      END
