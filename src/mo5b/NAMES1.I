C Version "@(#)$Id$ $Source$ $Date$ 

C  /NAMES1/:
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C  NAMFLG   C*8   Control flag names for prompts and error messages
C  NAMMMR   C*8   Prompts for replacement AMAR & JULMYR, by JDX & IV
C  NAMRVP   C*8   Name of RVP falling outside allowable range
C  NAMTEM   C*4   'mini' & 'maxi' prefixes for fuel tank temp error messages
C
C  Notes:
C
      CHARACTER*8 NAMFLG,NAMMMR,NAMRVP      
      CHARACTER*4 NAMTEM      
      COMMON /NAMES1/ NAMFLG(26),NAMMMR(2,3,8),NAMRVP(4),NAMTEM(2)
