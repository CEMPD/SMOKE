      SUBROUTINE DALSCN() 
C 
C  Return any storage that was dynamically allocated for scenario inputs. 
C 
C  Called by MOBILE. 
C 
C  Calls ... 
C 
C  Input on call: 
C 
C  Output on return: 
C 
C    Dynamically allocated storage is returned. 
C 
C  Changes : (Last change first) 
C 
C   15 Sep 99 @DynTel-ddj  Added ISTD to USEREFS array. 
C 
C  Module variables: 
C 
C    USEREFS - Array of derived type PTRNODE (see module EFUSER94). 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   INH     I    Loop control variable. 
C   IP      I    Loop control variable. 
C   ISR     I    Loop control variable. 
C   IV      I    Loop control variable. 
C   TEMP   Der.  Pointer to dynamically allocated derived type EMISSREC. 
C 
C  Notes: 
C 
C 
      USE EFUSER94 
      IMPLICIT NONE 
 
      INTEGER                   ::  INH 
      INTEGER                   ::  IP  
      INTEGER                   ::  ISR 
      INTEGER                   ::  IV 
      INTEGER                   ::  ISTD 
      TYPE (EMISSREC), POINTER  ::  TEMP 
 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
 
      IF (ALLOCATED(USEREFS)) THEN 
 
        DO ISR = 1,2 
          DO INH = 1,2 
            DO IP=1,3 
              DO IV = 1,5 
                DO ISTD = 1, 10 
                  DO WHILE(ASSOCIATED 
     &                    (USEREFS(IV,IP,INH,ISR,ISTD)%ROOT)) 
                    TEMP => USEREFS(IV,IP,INH,ISR,ISTD)%ROOT 
                    USEREFS(IV,IP,INH,ISR,ISTD)%ROOT => 
     &                    USEREFS(IV,IP,INH,ISR,ISTD)%ROOT%NEXT 
                    DEALLOCATE(TEMP) 
                  END DO  ! WHILE 
                END DO  ! ISTD DO 
              END DO  ! IV 
            END DO  ! IP 
          END DO  ! INH 
        END DO  ! ISR 
 
        DEALLOCATE(USEREFS) 
 
      END IF 
 
   99 RETURN 
 
      END SUBROUTINE DALSCN 
