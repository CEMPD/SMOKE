      SUBROUTINE ICCOPT(MY) 
C 
C  ICCOPT selects the pointer to the additive and multiplicative cold CO 
C  standard corrections vectors, plus sets the switches which determine when 
C  these corrections are applied. 
C 
C  Called by BIGTCF. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, declare MY with INTENT(IN) 
C  17 Nov 98 @DYnTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    argument list: MY 
C    common blocks: 
C    /FLAGS4/ COLDFG 
C    /TEMPC5/ NCC,MYCOL 
C 
C  Output on return: 
C 
C    common block: 
C    /TEMPC7/ ICOLD,ADDCSW,MULCSW 
C 
C  Notes: 
C 
C  ICCOSW was added for the cold CO standard(s) version and MOBILE4.1. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'TEMPC5.I' 
      INCLUDE 'TEMPC7.I' 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IG 
C 
      ICOLD=0 
      IF(MY.LT.MYCOL(1)) GOTO 20 
      IF(COLDFG.EQ.2) ADDCSW=2 
      IF(COLDFG.EQ.3) MULCSW=2 
C 
      DO 10 IG=1,NCC 
      ICOLD=NCC-IG+1 
      IF(MY.GE.MYCOL(ICOLD)) GOTO 20 
   10 CONTINUE 
C 
C  Either MY.LT.MYCOL(1) or get a match. 
C 
   20 RETURN 
      END 
