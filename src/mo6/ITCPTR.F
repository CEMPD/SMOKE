      FUNCTION ITCPTR(MY,IP,IVTEMP) 
C 
C  ITCPTR selects the pointer to the temperature correction coefficients. 
C 
C  Called by BIGCFX. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, created dummy model year variable, MYQTR, to be 
C            passed to QUITER 
C  24 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IP,IVTEMP 
C    common block: /TEMPC2/ MAXTCF,MYGTCF 
C 
C  Output on return: 
C 
C    function: ITCPTR 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'TEMPC2.I' 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IG, MYQTR, IP, IVTEMP, ITCPTR 
C 
      DO 10 IG=1,MAXTCF 
      IF(MY.LE.MYGTCF(IG,IP,IVTEMP)) GOTO 20 
   10 CONTINUE 
C 
C  Set default. 
C 
      IG=1 
      MYQTR=MY 
      CALL QUITER(0.,15,61,MYQTR) 
   20 ITCPTR=IG 
C 
      RETURN 
      END 
