      SUBROUTINE CHK94UEF(IV,IP,MY,ISR,INH,ISTD,KMILES,USERINP, 
     &                    USEREF) 
C 
C  CHK94UEF returns user-supplied emission factors for the 94+ light 
C  duty gasoline vehicles and trucks. 
C 
C  Called by EF_POST81. 
C 
C  Calls  ...  
C 
C  Changes: (Last change first) 
C 
C    15 Sep 99 @DynTel-ddj 2-000  Added ISTD index to USEREFS array. 
C 
C  Input on call or from calls: 
C 
C    argument list: INH,IP,ISR,ISTD,IV,KMILES,MY 
C    common blocks: 
C 
C  Output on return: 
C 
C    argument list: USERINP,USEREF. 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C  TEMP     R    Pointer to defined type EMISSREC (see module EFUSER94). 
  
      USE EFUSER94 
      IMPLICIT NONE 
 
      ! Declare parameters 
      INTEGER       :: INH 
      INTEGER       :: IP 
      INTEGER       :: ISR 
      INTEGER       :: ISTD 
      INTEGER       :: IV 
      REAL          :: KMILES 
      INTEGER       :: MY 
      REAL          :: USEREF 
      LOGICAL       :: USERINP 
 
      ! Declare local variables 
      TYPE (EMISSREC), POINTER  ::  TEMP 
 
      USERINP = .FALSE. 
      USEREF = -100.0 
C 
      IF (.NOT. ALLOCATED(USEREFS)) GO TO 99 
C 
      TEMP => USEREFS(IV,IP,INH,ISR,ISTD)%ROOT 
 
      DO WHILE(ASSOCIATED(TEMP)) 
        IF (MY >= TEMP%MYFIRST .AND. MY <= TEMP%MYLAST .AND. 
     &    ISTD == TEMP%STDCAT ) THEN 
          USERINP = .TRUE. 
          USEREF = TEMP%ZML + KMILES * TEMP%DR1 
          EXIT 
        ENDIF 
        TEMP => TEMP%NEXT 
      END DO  ! WHILE (ASSOCIATED) 
 
99    RETURN 
 
      END SUBROUTINE CHK94UEF 
 
