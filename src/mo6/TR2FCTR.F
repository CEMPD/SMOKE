      REAL FUNCTION TR2FCTR(MY,ACTUAL_MY,IVGAS,IPASSFAIL) 
C 
C    This function determines weighting for phase in of tier2 vehicles. These 
C    multiplicative factors are from tier2 emissions reductions and fractions 
C    of vehicles introduced into the fleet of vehicles in a given model year. 
C    The fractions are stored in the two arrays T2_EVAP and T2_RED in 
C    BD_TIER2.FOR 
C 
C    Called by DIRNLDAY, HOTSOK, and 
C 
C    Changes: (Last change first) 
C 
C    27 Feb 01 @EPA-BAg bug314 revised tier2 effects on heavy duty vehicles 
C 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: MY,ICY,IVGAS 
C    functions : 
C    modules   : 
C 
C    common blocks: 
C    /EVAPHD/ HDWGT 
C    /FLAGS5/ HDR2007 
C    /TIER2/  T2_EVAP,T2_RED, T2_HDWGT 
C 
C 
C  Output on return: 
C 
C    function     : TR2FCTR 
C    argument list: 
C    modules      : 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name          Type              Description 
C  ------         ----  ------------------------------------------------------- 
C  NAGE            I    Age relative to 2004 (starts at1); vehicle age for 
C                       tier2 fractions 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAPHD.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'TIER2.I' 
C 
      INTEGER, INTENT(IN) :: IPASSFAIL 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: IVGAS 
C 
      INTEGER NAGE 
      INTEGER, SAVE :: IVGAS_MC = 14 
      INTEGER, SAVE :: INO =0 
      INTEGER, SAVE :: IPASSBOTH = 1 
C 
      NAGE = ACTUAL_MY - 2004 + 1 
      TR2FCTR=1.0 
C 
      IF( IVGAS==IVGAS_MC .OR. 
     *   (IVGAS>5.AND.HDR2007==INO) .OR. 
     *   (IVGAS<6.AND.IPASSFAIL/=IPASSBOTH)) GOTO 99 
C 
      MODEL_YEARS : SELECT CASE (MY) 
        CASE (2004:2015) MODEL_YEARS 
          TR2FCTR=T2_EVAP(NAGE,IVGAS)*(1.0-T2_RED(IVGAS)) 
     &              + (1.0 - T2_EVAP(NAGE,IVGAS)) 
C 
        CASE (2016:) MODEL_YEARS 
          TR2FCTR=1.0-T2_RED(IVGAS) 
C 
        CASE DEFAULT MODEL_YEARS 
          TR2FCTR=1.0 
C 
      END SELECT MODEL_YEARS 
C 
  99  IF(MY<2004.OR.(IVGAS>=6.AND.MY<2008)) TR2FCTR=1.0 !be sure tier2 effects don't exist 
C                                                        prior to 2004 for light duty gas 
C                                                        and prior to 2008 for heavy duty gas 
C 
      RETURN 
      END FUNCTION TR2FCTR 
