      SUBROUTINE ATBIGCFX(IDX,MY,IV,IP,IH) 
C 
C  ATBIGCFX computes correction factors applied to composite emission 
C  factors, e.g., speed, optional (humidity, trailer towing, a/c), 
C  temperature, by-bag RVP, operating mode, methane and tampering. 
C 
C  Called by ATHRLOOP. 
C 
C  Calls BIGIN1 and BIGTCF. 
C 
C  Changes : (Last change first) 
C 
C  27 Feb 02 New MOBILE6.3 routine based on 27 Jul 01 version of BIGCFX. 
C      
C 
C  Input on call: 
C 
C    Constants:  
C    argument list: IDX,JDX,MY,IV,IP,IH 
C    common blocks: 
C    /IVTYPE/ VM6M5B 
C    /MYRCAL/ TF 
C 
C  Output on return: 
C 
C    common blocks: 
C    /TEMPC4/ LOWCO 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE  
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'TEMPC4.I' 
      INCLUDE 'IVTYPE.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IH 
C 
C     Declare local variables 
C 
      INTEGER IVTMP5 
C 
C  Initialize variables to be used in upcoming loops. 
C 
      CALL BIGIN1(IH) 
C 
C  Skip calculations if travel fraction is 0.0, zeroing emission factor 
C 
      IF(TF(IDX,IV).EQ.0.0) GOTO 999 
C 
C  Reset CO offset switch and calculate temperature correction factor (TCF). 
C  There is no TCF in MOBILE4 for diesel IV (IVTNP5=0).  When the ambient 
C  temperature is the FTP standard 75 degrees (MDLOHI=0), TCF represents the 
C  RVP CF for RVP MYG cases 2 and 3 for LDGV/T & HDGV, otherwise it remains the 
C  initialized default neutral value (1.0) for all IV (1-8). 
C 
      LOWCO=1 
      IVTMP5=VM6M5B(IV) 
      IF(IVTMP5.GE.1) CALL BIGTCF(IP,MY,IV,IVTMP5) 
C 
  999 RETURN 
      END 
