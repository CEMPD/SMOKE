      REAL FUNCTION SCFSP(IP,IV,MY,SP,INERR)
C
C SFCSP calculate the speed correction factor for
C specific speed using MOBILE5 methods
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C 17 Aug 99 @DynTel-ddj 2-698 Changed to a function and removed IDX
C               and ISR from parameter list.
C 12  May  1999 @Dyntel-HQ 2-698 simplify the code for calculation of SCF
C  1 March 1999 @DynTel-HQ 2-698 calculate SFC for specific speed
C                          using MOBILE5 method.
C    Called by SCFSLOW
C
C    Calls BIGSC2, BIGSC3.
C
C    Constants: 
C
C    argument list: ICY
C    common blocks:
C    /SPEED4/ LB1STS,IVA
C
C
      IMPLICIT NONE
      INCLUDE 'SPEED4.I'
C
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(INOUT) :: INERR
      REAL BIGSC2, BIGSC3
      REAL SP
C
      IVA=IV
C
C  For Speed cf, either use old MOBILE3 algorithm,
C  Pre 1980 LDGV and LDGT
C
       IF(MY.LT.LB1STS(IV)) THEN
          SCFSP = BIGSC2(IP,IV,MY,SP,INERR)
C
C  Or use new equation form, based on speed and pollutant type.
C  Low, Mid and High Speed for all three IP's
C  Post 1980 LDGV and LDGT
C
        ELSE
          SCFSP = BIGSC3(IP,MY,SP,INERR)
        ENDIF
C
100   RETURN
      END
