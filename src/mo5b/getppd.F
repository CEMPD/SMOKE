C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETPPD(INERR,*)
C
C  GETPPG verifies Purge and Pressure  program types. If one is Hybrid and 
C  the other is not ERROR. 
C  Reads the PPEFF.D file. The first descriptive record write into report. 
C
C  Called by ONESEC.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /IOUCOM/ IOUPFF
C    /PPEFFD/ PPEFF(25,2),PPEFF1(25,2),PPDES
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /PPEFFD/ PPEFF(25,2),PPDES,RBHPPF 
C
C  Notes:
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Subroutine was daveloped in accordance with this task request.  
C 
      INCLUDE 'IOUCOM.I'
      INCLUDE 'PPEFFD.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'PRSCH2.I'
C
C  Verify Purge and Pressure  program types. If one is Hybrid and 
C  the other is not ERROR. 
C
      IF ( (PRSFLG .EQ. 2) .AND. (PRGFLG .EQ. 2) .AND.
     *   (PRGPGM .NE. PRSPGM) ) THEN
           CALL QUITER(0.,0,162,INERR)
           GOTO 99
      ENDIF	  
      
C
C  Read data file PPEFF.D and load values into real array PPEFF
C
      OPEN (IOUPFF, FILE='PPEFF.D',ERR=98)
C      
      READ (IOUPFF,'(A)',ERR=98) PPDES 
C
      DO 10 J1 = 1,2
      READ (IOUPFF,100,ERR=98) (PPEFF(I1,J1), I1 = 1,25)
 100    FORMAT (25F5.3)
  10  CONTINUE
      CLOSE(IOUPFF)
C      
C  Set the Pressure and Purge flags to Test Only and RBHPPF flag to 2         
C
      IF  (PRSPGM .EQ. 3) PRSPGM=1
      IF  (PRGPGM .EQ. 3) PRGPGM=1
      RBHPPF=2
C
  90  RETURN 
C
  98  CALL QUITER(0.,0,161,INERR)
  99  RETURN1
      END
