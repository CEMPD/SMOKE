C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETEFR(INERR,*)
C
C  GETEFR reads the I/M Program Effectiveness Record, which contains the
C    user-supplied values for the effectivenes of test-and-repair I/M
C    programs (three values), test-and-repair ATP programs (one value),
C    and test-and-repair Purge/Pressure programs (one value). The record 
C    appears in the One-time Data section if a flag is set on the I/M
C    Program Descriptive Record.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C     /ATPAR4/  EFF, EFFSAV
C     /FLAGS1/  PROMPT
C     /FLAGS2/  IMFLAG
C     /FLAGS3/  ATPFLG
C     /IMPAR2/  ITEST
C     /IMPAR5/  DISCNT, DCNTSV, ATPUEF, PPGUEF
C     /IMPAR6/  IMEFFR, INTYP
C     /IOUCOM/  IOUASK, IOUGEN
C     /PRGCH2/  PRGFLG
C     /PRSCH2/  PRSFLG
C
C  Output on return:
C
C    parameter: INERR
C    common blocks:
C    /ATPAR4/  EFF
C    /IMPAR4/  ATPUEF, PPGUEF, DISCNT
C
C  Local variable / array dictionary:
C
C   Name   Type                      Description
C  ------  ----  ----------------------------------------------------------
C  ATPADJ   R    Adjustment to ATP effectiveness due to presence of IM240.
C  AVGIMF   R    Average of IM effectiveness rates = (HC + CO + NOx) / 3.0
C  EFFREC   C    Character representation of Effectiveness Rec. (The actual
C                record is variable length - some fields are optional.)
C  IATP     I    One dimension of EFF array - see BLOCK DATA 25.
C  ICOL     I    One dimension of EFF array - see BLOCK DATA 25.
C  IROW     I    One dimension of EFF array - see BLOCK DATA 25.
C
C  Notes:
C
C 
      INCLUDE 'ATPAR4.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'PRSCH2.I'
C
      CHARACTER  EFFREC*24
      REAL       ATPADJ, AVGIMF
      INTEGER    IATP, ICOL, IROW
C
C     Initialize the effectiveness rates for test-and-repair ATP programs
C     and purge/pressure programs to the default values. These values
C     wil be used unless the user enters alternate values.
C
      PPGUEF = 0.50
      ATPUEF = 0.0
C
C     If the IMEFFR flag was set, read the user-supplied values.
C
      IF (IMEFFR.EQ.2) THEN
C
C       The record should contain either three or five values. The first 
C       three are the rates for HC, CO, and NOx for test-and-repair I/M 
C       programs (use the same values for both computerized and manual 
C       I/M programs). The next two values are optional. They are the 
C       effectiveness rates for test-and-repair ATP and purge/pressure 
C       programs. (In prompting mode, ask the user to supply these values.)
C
        IF(PROMPT.EQ.2) THEN
          WRITE(IOUASK, 100)
  100     FORMAT(' Enter three values for the effectiveness of',
     &           ' test-and-repair I/M programs in', 
     &          /' reducing HC, CO, and NOx emissions. You may also',
     &           ' enter values for the effectiveness',
     &          /' of test-and-repair ATP programs and the', 
     &           ' effectiveness of test-and-repair', 
     &          /' purge/pressure programs, a total of five numbers.',
     &           ' (The last two are optional.)',
     &          /' The input format is (4(F4.2, 1X), F4.2).')
        ENDIF
C
        EFFREC = ' '
        READ(IOUGEN,120,ERR=98,END=98) EFFREC
C
        READ(UNIT=EFFREC, FMT=121, ERR=98) (DISCNT(IP,2),IP=1,3)
C
C       Use the same values for manual and computerized test-and-repair
C       I/M programs.
C
        DISCNT(1,3) = DISCNT(1,2)
        DISCNT(2,3) = DISCNT(2,2)
        DISCNT(3,3) = DISCNT(3,2)
C
C       Look for effectiveness rates for ATP and Purge/Pressure programs.
C
        IF (EFFREC(16:19).NE.' ') THEN
          READ(UNIT=EFFREC(16:19), FMT=122, ERR=98) ATPUEF
        ENDIF
C
        IF (EFFREC(21:24).NE.' ') THEN
          READ(UNIT=EFFREC(21:24), FMT=122, ERR=98) PPGUEF
        ENDIF
C
  120   FORMAT(A24)
  121   FORMAT(2(F4.2, 1X), F4.2)
  122   FORMAT(F4.2)
C
C       The user-entered value for ATPUEF will be ignored unless an
C       ATP program is explicitly modeled. The value for PPGUEF will
C       be ignored unless a Purge or Pressure program is modeled.
C
        IF (ATPFLG.EQ.1) ATPUEF = 0.0
C
        IF (PRGFLG.EQ.1.AND.PRSFLG.EQ.1) PPGUEF = 0.50
C
      ENDIF
C
C     Reset the effectiveness rates for test-and-repair ATP programs. We 
C     will interpolate betwen the default values for test-only and test-
C     and-repair programs, based on the value supplied by the user. If the 
C     user specified a value of 0.0 for ATPUEF (or if the user supplied no
C     value and ATPUEF defaults to 0.0) the values we place in EFF will
C     be the default values for test-and-repair ATP programs. First, make
C     sure ATPUEF is between 0.0 and 1.0, inclusive.
C
      IF (ATPUEF.LT.0.0.OR.ATPUEF.GT.1.0) THEN
        CALL QUITER(ATPUEF, 0, 167, INERR)
        GOTO 40
      ENDIF
C
C     Reset EFF (to default values if no alternate rates were entered).
C
      DO 30 ICOL = 1,3
C
        DO 30 IATP = 1,2
C
          DO 30 IROW = 1,15
            EFF(2,ICOL,IATP,IROW) = EFFSAV(2,ICOL,IATP,IROW) + ATPUEF*
     & (EFFSAV(1,ICOL,IATP,IROW)-EFFSAV(2,ICOL,IATP,IROW))
   30 CONTINUE
C
C     Return if no atlernate effectiveness rates were entered.
C
      IF (IMEFFR.NE.2) GOTO 90
C
C     IM240 programs are required to cover three ATP inspection types. If 
C     the user is claiming greater than default effectiveness for their
C     test-and-repair IM240 program, the effectiveness of these three ATP
C     inspection types must be enhanced as well. Map the enhanced IM240
C     effectiveness, which varies from 0.50 to 1.00 for three different
C     pollutants, into enhanced effectiveness for ATP programs, which 
C     varies from 0.00 to 1.00. Skip the adjustment if the I/M program
C     effectiveness is no greater than the default value, or if the I/M
C     program effectiveness is less than the claimed ATP effectiveness,
C     or if there are no test-and-repair IM240 programs. 
C
      AVGIMF = (DISCNT(1,2) + DISCNT(2,2) + DISCNT(3,2)) / 3.0
C
      IF (AVGIMF.LE.0.50) GOTO 40
C
      ATPADJ = 2.0 * (AVGIMF - 0.50)
C
      IF (ATPADJ.LE.ATPUEF) GOTO 40
C
      IF (ITEST(1).EQ.4.AND.INTYP(1).GT.1) GOTO 32
C
      IF (IMFLAG.NE.3.AND.IMFLAG.NE.5) GOTO 40
C
      IF (ITEST(2).NE.4.OR.INTYP(2).EQ.1) GOTO 40
C
C     Adjust rows 3, 7, and 9 in EFF.
C
   32 CONTINUE
C
      DO 33 ICOL = 1, 3
C
        DO 33 IATP = 1, 2
          EFF(2, ICOL, IATP, 3) = EFFSAV(2, ICOL, IATP, 3) +
     & ATPADJ * (EFFSAV(1,ICOL,IATP,3)-EFFSAV(2,ICOL,IATP,3))
          EFF(2, ICOL, IATP, 7) = EFFSAV(2, ICOL, IATP, 7) +
     & ATPADJ * (EFFSAV(1,ICOL,IATP,7)-EFFSAV(2,ICOL,IATP,7))
          EFF(2, ICOL, IATP, 9) = EFFSAV(2, ICOL, IATP, 9) +
     & ATPADJ * (EFFSAV(1,ICOL,IATP,9)-EFFSAV(2,ICOL,IATP,9))
   33 CONTINUE
C
C     Make sure all values entered by the user are within range.
C
   40 CONTINUE
C
      DO 50 IP = 1,3
C
        IF (DISCNT(IP,2).LT.0.0.OR.DISCNT(IP,2).GT.1.0)
     & CALL QUITER(DISCNT(IP,2), 0, 157, INERR)
  50  CONTINUE
C
      IF (PPGUEF.LT.0.0.OR.PPGUEF.GT.1.0)
     & CALL QUITER(PPGUEF, 0, 168, INERR)
C
C     Convert PPGUEF from an 'effectiveness' rate, as entered by the
C     user, to a 'discount' rate, as expected by the code.
C
      PPGUEF = 1.0 - PPGUEF
C
   90 RETURN
C            
C   This is the exit for general FORTRAN I/O errors.
C
   98 CALL QUITER(0.0, 0, 113, INERR)
      RETURN 1
C
      END     
