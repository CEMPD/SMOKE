      SUBROUTINE TEX1(IP,IV,IDX,ISR,MY) 
C 
C    Computes the tampering offset. 
C 
C    13 Sep 99 New Function 
C 
C  Called by EFCALX. 
C 
C  Calls: 
C 
C  Changes: (Last change first) 
C 
C  07 Jul 00 AIR Task 05: Removed unused KKSTD from parameter list. 
C  01 Jun 00 @EPA-djb Bug123, Changed TAMFLG=0 case for HDGV to 0.0 for 
C                     all model years. 
C 
C  Input on call: 
C    argument list: IP,IV,IDX,ISR,BASEF 
C    common blocks: 
C    /FLAGS1/ TAMFLG 
C    /TAMOUT/ TAMBAG,TAMBAS 
C 
C 
C  Output on return: 
C    common blocks: 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   BASEF   R    basic emission returned from BEFUN 
C   IVT     I    vehicle type mapping to tampering arrays 
C   TEXB    R    tampering excess for the case with no anti-tampering 
C   TEXC    R    tampering excess with anti-tampering effects 
C   TEXD    R    Effect of anti-tampering on the base tampering excess 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'TAMOUT.I' 
      INCLUDE 'IVTYPE.I' 
C 
C    Declare called functions 
C 
C    Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: MY 
C 
C    Declare local variables 
C 
      INTEGER  IVT 
      REAL     TEXB 
      REAL     TEXC 
      REAL     TEXD 
C 
C    Only gasoline fueled vehicles which are not motorcycles 
C    get tampering excess. 
C 
      TEX(IV,IP,ISR,IDX) = 0.0 
C 
C    Map the vehicle types into one of three categories. 
C 
      IF(VVGAS(IV).EQ.0)   RETURN 
      IVT = IV 
      IF(IV.EQ.25) IVT = 14 
C 
C    Start emissions are calculated from the difference between Bag 1 
C    and Bag 3 emissions and converted to grams per start. 
C 
C    Running emissions are calculated from the LA4 VMT weighted 
C    combination of Bag 2 and Bag 3 (in grams per mile). 
C 
C    First calculate the base (uncontrolled) tampering excess. 
C 
      IF(ISR.EQ.2) THEN 
       TEXB=(TAMBAS(1,IP,IDX,IVT)-TAMBAS(3,IP,IDX,IVT))*3.49 
      ELSE 
       TEXB=0.521*TAMBAS(2,IP,IDX,IVT)+0.479*TAMBAS(3,IP,IDX,IVT) 
      ENDIF 
      IF(TEXB.LT.0.0) TEXB=0.0 
C 
C    Next calculate the controlled tampering excess. 
C 
      IF(ISR.EQ.2) THEN 
       TEXC=(TAMBAG(1,IP,IDX,IVT)-TAMBAG(3,IP,IDX,IVT))*3.49 
      ELSE 
       TEXC=0.521*TAMBAG(2,IP,IDX,IVT)+0.479*TAMBAG(3,IP,IDX,IVT) 
      ENDIF 
      IF(TEXC.LT.0.0) TEXC=0.0 
C 
C 
C    The difference between the controlled and uncontrolled case is the 
C    effect of anti-tampering programs on the basic exhaust emission 
C    rate. 
C 
      TEXD=TEXB-TEXC 
      IF(TEXD.LT.0.0) TEXD=0.0 
C 
C    Pre-1981 model years will have the remaining tampering excess, 
C    after any anti-tampering program effects, *added* to the 
C    basic exhaust emission rate.  1981 and newer model years will 
C    have the effect of anti-tampering programs *subtracted* from the 
C    basic exhaust emission rate (since tampering is already included). 
C    Vehicles certified to emission standards more stringent than Tier 1 
C    (KKSTD>1) are assumed to have no tampering excess at all. 
C    It is assumed that the effect of anti-tampering programs cannot 
C    equal more than a fixed fraction of the basic exhaust emission rate. 
C 
      IF(IV.LE.5) THEN 
        IF(MY.LT.1981) THEN 
          TEX(IV,IP,ISR,IDX)=TEXC 
        ELSEIF(MY.GE.1981) THEN 
          TEX(IV,IP,ISR,IDX)=TEXD*(-1.0) 
        ENDIF 
      ELSEIF(IV.GE.6) THEN 
          TEX(IV,IP,ISR,IDX)= 0.206*TAMBAG(1,IP,IDX,IVT) + 
     *                        0.521*TAMBAG(2,IP,IDX,IVT) + 
     *                        0.273*TAMBAG(3,IP,IDX,IVT) 
      ENDIF 
C 
      IF(TAMFLG.EQ.0)  THEN 
        IF(IV.LE.5) THEN 
          IF(MY.LT.1981) THEN 
             TEX(IV,IP,ISR,IDX) = 0.0 
          ELSEIF(MY.GE.1981) THEN 
             TEX(IV,IP,ISR,IDX) = TEXB*(-1.0) 
          ENDIF 
        ELSEIF(IV.GE.6) THEN 
          TEX(IV,IP,ISR,IDX)= 0.0 
        ENDIF 
      ENDIF 
C 
C  All 1996 and newer vehicles are assumed to have no tampering effects. 
C 
      IF(MY.GE.1996) THEN 
        TEX(IV,IP,ISR,IDX) = 0.0 
      ENDIF 
C 
      RETURN 
      END 
