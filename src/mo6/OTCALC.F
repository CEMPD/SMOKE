      REAL FUNCTION OTCALC(ID,IAYBTR,IHGEVP,IAYEVP) 
C 
C  OTCALC is a REAL-valued function that returns a component of the Overall 
C  Tampering Rate (OTR) equation. 
C 
C  Called by EVPEME 
C 
C  Changes: (Last change first) 
C 
C 
C  Input on Call: 
C 
C    argument list: ID,IAYBTR,IHGEVP,IAYEVP 
C    common blocks: 
C    /EVPDAT/ EVPAER 
C    TAMEQ4/ BTR 
C 
C  Output on return: 
C 
C    function: OTCALC 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IAYBTR   I    IAY belonging to the BTR portion of the equation 
C  IHGEVP   I    IHG belonging to the EVPAER portion of the equation 
C  IAYEVP   I    IAY belonging to the EVPAER portion of the equation 
C 
C  Notes: 
C 
C  None. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'EVPDAT.I' 
      INCLUDE 'TAMEQ4.I' 
C 
C     Declare parameter list 
C 
      INTEGER IAYBTR 
      INTEGER IAYEVP 
      INTEGER ID 
      INTEGER IHGEVP 
C 
      OTCALC=BTR(ID,IAYBTR)*EVPAER(IHGEVP,IAYEVP) 
C 
      RETURN 
C 
      END 
