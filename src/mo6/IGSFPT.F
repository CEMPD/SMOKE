      FUNCTION IGSFPT(MY) 
C 
C  IGSFPT selects the pointer to the gas/diesel sales fractions. 
C 
C  Called by DSFIDX. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  25 Jan 2000 @EPA-djb reflect changes in array. 
C  18 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IV 
C    common blocks: 
C    /DSFCOM/ MAXGSF 
C    /REGION/ IREJN 
C 
C  Output on return: 
C 
C    function: IGSFPT 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
      INCLUDE 'DSFCOM.I' 
      INTEGER MY, IGSFPT 
C 
      IF(MY.LE.1972) THEN 
        IGSFPT = 1 
        RETURN 
      ELSEIF(MY.GE.1996) THEN 
        IGSFPT = 25 
        RETURN 
      ELSE 
        IGSFPT = MY - 1971 
      ENDIF 
C 
      RETURN 
      END 
