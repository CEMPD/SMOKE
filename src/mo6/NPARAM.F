      SUBROUTINE NPARAM(DDISTP,ID,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE) 
C 
C  NPARAM is used to get the correct parameters and leaves the old ones 
C  alone as a function of the DDISTP(ID). 
C 
C  Called by ATPAER 
C 
C  Changes: (Last change first) 
C 
C  25 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C    argument list: DDISTP,ID,NAPGM,APGM,NBFQT,BFQT,NCRATE,CRATE 
C    common blocks: 
C    /ATPAR2/ ATPGM,ATPFQT,CRATP,DISTYP 
C 
C  Output on return: 
C    argument list: NAPGM,NBFQT,NCRATE 
C 
C 
C  Notes: 
C 
C  This subroutine was added to MOBILE5.3 to facilitate the IM240 test 
C  effects on the ATP effectiveness.  It is called to flip back and 
C  forth from ATP to IM240 parameters.  The parameters switched are 
C  inspection type, frequency and compliance. 
C 
C   5 Mar 96 @DynTel-MLA 2-611 Fix an old bug uncovered while testing. The 
C      IM240 parameters should not override the ATP parameters for 
C      inspection types that are covered by both programs if the ATP is a 
C      test-only program that would deliver greater benefits. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATPAR2.I' 
C 
      REAL NCRATE, CRATE 
      INTEGER DDISTP(8), APGM, BFQT, NAPGM, NBFQT, ID 
C 
      NAPGM=ATPPGM 
      NBFQT=ATPFQT 
      NCRATE=CRATP 
      IF(DDISTP(ID).EQ.3) THEN 
C 
C       If DDISTP(ID) is 3, inspection type ID is covered by an IM240 
C       program. NAPGM, NBFQT, and NCRATE have been set by default above  
C       to the ATP program paramters. This block of code switches back to  
C       the I/M program parameters APGM, BFQT, and CRATE, unless this  
C       inspection type is covered by a test-only ATP that would give  
C       greater benefits. 
C 
        IF (DISTYP(ID).EQ.2.AND.ATPPGM.EQ.1) THEN 
C 
C         There is a test-only ATP. The I/M program can still win if 
C         the I/M is also test-only and the I/M parameters are more 
C         beneficial than the ATP paramters. 
C 
          NAPGM = ATPPGM 
          IF (APGM.EQ.1.AND.BFQT.LT.ATPFQT) NBFQT = BFQT 
          IF (APGM.EQ.1.AND.CRATE.GT.CRATP) NCRATE = CRATE 
C 
        ELSE 
C 
C         No ATP, or ATP is not test-only. Use I/M program parameters. 
C 
          NAPGM = APGM 
          NBFQT = BFQT 
          NCRATE = CRATE 
C 
        ENDIF 
C 
      ENDIF 
      RETURN 
      END 
