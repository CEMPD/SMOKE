C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTATP
C
C  OUTATP prints the Anti-Tampering Program (ATP) parameters.
C  and the Pressure and Purge Programs if applicable
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /FLAGS3/ ATPFLG
C    /IOUCOM/ IOUREP
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
C    /STRING/ NOYES
C
C  Local array subscripts:
C
C  PGM(2)    -  PGM ( ATPPGM or PRSPGM or PRGPGM )
C  FQT(2)    -  FQT ( ATPFQT or PRSFQT or PRGFQT )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  FQT      C*8  Fun program frequency names
C  PGM      C*13 Fun program type names
C
C
C  Notes:
C
C  Was modified in MOBILE5.3 to include output for ATPFLG=1-8
C
C
      CHARACTER*4 NOYES,COMMA,PERIOD
      CHARACTER*1 COLON
      INTEGER ATPPGM,ATPFQT,DISTYP,ATPFLG
      INTEGER PRSPGM,PRGPGM,PRSFQT,PRGFQT,PRSFLG,PRGFLG
C
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
      COMMON /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP(8)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS(4)
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
      COMMON /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG(4)
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
C
      CHARACTER*8 FQT(2)
      CHARACTER*13 PGM(2)
C
C                1234567890123   1234567890123
      DATA PGM /'Test Only    ','Test & Repair'/
      DATA FQT /'Annual  ','Biennial'/
C
C  Write Header
C
      WRITE(IOUREP,200)
  200 FORMAT('0Functional Check Program Description:'/
     * '0Check Start   Model Yrs  Vehicle Classes Covered        ',
     * 'Inspection        Comp'/
     * '       (Jan1)  Covered    LDGV   LDGT1  LDGT2  HDGV    ',
     * 'Type      Freq      Rate'/)
C
C  Pressure when PRSFLG is 2.
C
       IF(PRSFLG.NE.1)
     * WRITE(IOUREP,210) LPRSY, LPR1ST, LPRLST,
     * (NOYES(LVPRS(IVTAM)), IVTAM=1, 4),
     * PGM(PRSPGM), FQT(PRSFQT), CRPRS
  210  FORMAT(' Press ',I4, 4X, I4, '-', I4, 2X, 3(A3, 4X), A3, 1X,
     * A13, 1X, A8, 1X, F5.1, '%')
C
C  Purge when PRGFLG is 2.
C
       IF(PRGFLG.NE.1)
     * WRITE(IOUREP,220) LPGSY, LPG1ST, LPGLST,
     * (NOYES(LVPRG(IVTAM)), IVTAM=1, 4),
     * PGM(PRGPGM), FQT(PRGFQT), CRPRG
  220  FORMAT(' Purge ',I4, 4X, I4, '-', I4, 2X, 3(A3, 4X), A3, 1X,
     * A13, 1X, A8, 1X, F5.1, '%')
C
C  ATP when ATPFLG is 2.
C
       IF(ATPFLG.NE.1)
     * WRITE(IOUREP,230) LAPSY, LAP1ST, LAPLST,
     * (NOYES(LVTFLG(IVTAM)), IVTAM=1, 4),
     * PGM(ATPPGM), FQT(ATPFQT), CRATP,
     *                 (NOYES(DISTYP(IDUSER)), IDUSER=1, 8)
  230  FORMAT(' ATP   ',I4, 4X, I4, '-', I4, 2X, 3(A3, 4X), A3, 1X,
     * A13, 1X, A8, 1X, F5.1, '%'/
     * '0Air pump system disablements:       ', A4,
     * ' Catalyst removals:                  ', A4/
     * ' Fuel inlet restrictor disablements: ', A4,
     * ' Tailpipe lead deposit test:         ', A4/
     * ' EGR disablement:                    ', A4,
     * ' Evaporative system disablements:    ', A4/
     * ' PCV system disablements:            ', A4,
     * ' Missing gas caps:                   ', A4)
C
      RETURN
      END
