C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION OTCALC(ID,IAYBTR,IHGEVP,IAYEVP)
C
C  OTCALC is a REAL-valued function that returns a component of the Overall
C  Tampering Rate (OTR) equation.
C
C  Called by EVPEME
C
C  Input on Call:
C
C    parameter list: ID,IAYBTR,IHGEVP,IAYEVP
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
      INCLUDE 'EVPDAT.I'
      INCLUDE 'TAMEQ4.I'
C
      OTCALC=BTR(ID,IAYBTR)*EVPAER(IHGEVP,IAYEVP)
C
      RETURN
C
      END
