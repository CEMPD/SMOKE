C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /IMPAR5/ ATPUEF,PPGUEF
C    /IOUCOM/ IOUREP
C    /PPEFFD/ RBHPPF
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
C  ATPEFF   R    ATP program effectiveness to display on output
C  FQT      C*8  Fun program frequency names
C  PGM      C*13 Fun program type names
C  PRGEFF   R    Purge program effectiveness to display on output
C  PRSEFF   R    Pressure program effectiveness to display on output
C
C
C  Notes:
C
C  Was modified in MOBILE5.3 to include output for ATPFLG=1-8
C
C
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Print the pressure/purge  descriptive record
C 29 Feb 96 @DynTel-MLA 2-611 Added output for the 'effectiveness' rates
C     of ATP, Purge, and Pressure Check programs.
C
C  1 Apr 96 @DynTel-MLA 2-611 (part 2) The 'effectiveness' rate is relevant
C     only for test-and-repair programs. Display 1.00 for other prog types.
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'PPEFFD.I'
      INCLUDE 'PRSCH1.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'PRGCH1.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'STRING.I'
C
      REAL          ATPEFF
      CHARACTER*8   FQT(2)
      CHARACTER*13  PGM(3)
      REAL          PRGEFF
      REAL          PRSEFF
C
C                1234567890123   1234567890123
      DATA PGM /'Test Only    ','Test & Repair', 'RBH Program'/
      DATA FQT /'Annual  ','Biennial'/
C
C  Write Header
C
      WRITE(IOUREP,200)
  200 FORMAT('0Functional Check Program Description:'/
     * '0Check Start Model Yrs Vehicle Classes Covered ',
     * '   Inspection        Comp   Eff'/
     * '      (Jan1) Covered   LDGV LDGT1 LDGT2 HDGV   ',
     * 'Type        Freq     Rate   Adj')
C
C  Retest Based Hybrid (RBH) Program.
C
       IF (RBHPPF .EQ. 1) THEN
         IPRS=PRSPGM
         IPRG=PRGPGM
       ELSE
         IPRS=3
         IPRG=3
       ENDIF                  
C
C  Pressure when PRSFLG is 2. NOTE: The test-and-repair 'discount' rate
C  and the test-and-repair 'effectiveness' rate are complementary. The
C  user enters an 'effectiveness' rate for PPGUEF, but we immediately
C  convert it to a 'discount' for internal use. For presentation to the
C  user, convert it back to an 'effectiveness' rate (EFF=1.0-DISC).
C  NOTE: Pressure test effectiveness is relevant only for test-and-repair
C  programs. For all other program types, we will display a value of 1.00.
C
       IF (PRSPGM.EQ.2) THEN
         PRSEFF = 1.0 - PPGUEF
       ELSE
         PRSEFF = 1.0
       ENDIF
C
       IF(PRSFLG.NE.1)
     * WRITE(IOUREP,210) LPRSY, LPR1ST, LPRLST,
     * (NOYES(LVPRS(IVTAM)), IVTAM=1, 4),
     * PGM(IPRS), FQT(PRSFQT), CRPRS, PRSEFF
  210  FORMAT(' Press ',I4, 2X, I4, '-', I4, 1X, A3, 2X, 2(A3, 3X),
     * A3, 1X, A13, 1X, A8, 1X, F5.1, '%',2X,F4.2)
C
C  Purge when PRGFLG is 2.
C  NOTE: Purge test effectiveness is relevant only for test-and-repair
C  programs. For all other program types, we will display a value of 1.00.
C
       IF (PRGPGM.EQ.2) THEN
         PRGEFF = 1.0 - PPGUEF
       ELSE
         PRGEFF = 1.0
       ENDIF
C
       IF(PRGFLG.NE.1)
     * WRITE(IOUREP,220) LPGSY, LPG1ST, LPGLST,
     * (NOYES(LVPRG(IVTAM)), IVTAM=1, 4),
     * PGM(IPRG), FQT(PRGFQT), CRPRG, PRGEFF
  220  FORMAT(' Purge ',I4, 2X, I4, '-', I4, 1X, A3, 2X, 2(A3, 3X),
     * A3, 1X, A13, 1X, A8, 1X, F5.1, '%', 2X, F4.2)
C
C  ATP when ATPFLG is 2.
C  NOTE: ATP program effectiveness is relevant only for test-and-repair
C  programs. For all other program types, we will display a value of 1.00.
C
       IF (ATPPGM.EQ.2) THEN
         ATPEFF = ATPUEF
       ELSE
         ATPEFF = 1.0
       ENDIF
C
       IF(ATPFLG.NE.1)
     * WRITE(IOUREP,230) LAPSY, LAP1ST, LAPLST,
     * (NOYES(LVTFLG(IVTAM)), IVTAM=1, 4),
     * PGM(ATPPGM), FQT(ATPFQT), CRATP, ATPEFF,
     * (NOYES(DISTYP(IDUSER)), IDUSER=1, 8)
  230  FORMAT(' ATP   ',I4, 2X, I4, '-', I4, 1X, A3, 2X, 2(A3, 3X),
     * A3, 1X, A13, 1X, A8, 1X, F5.1, '%', 2X, F4.2/
     * '0Air pump system disablements:       ', A4,
     * ' Catalyst removals:                  ', A4/
     * ' Fuel inlet restrictor disablements: ', A4,
     * ' Tailpipe lead deposit test:         ', A4/
     * ' EGR disablement:                    ', A4,
     * ' Evaporative system disablements:    ', A4/
     * ' PCV system disablements:            ', A4,
     * ' Missing gas caps:                   ', A4)
C
C
C  Print the pressure/purge  descriptive record
C 
      IF (PPDES  .NE. ' ') THEN
        WRITE(IOUREP,'(1X,A)') PPDES
      ENDIF       
C      
      RETURN
      END
