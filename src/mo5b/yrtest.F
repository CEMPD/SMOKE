C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE YRTEST(IY,NERR,IY1,IY2,INERR)
C
C  YRTEST converts 2-digit input year to 4-digit year.
C
C  Called by GETATP, GETERQ, GETIMP, GETNAT, PARSEC and RDCITY.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: IY,NERR,IY1,IY2,INERR
C    common blocks:
C    /YEARS4/ IYSTART
C
C  Output on return:
C
C    parameter list: IY,INERR
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  II       I    IY converted to 4 digit year
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'YEARS4.I'
C
      II=IY+1900
      IF(II.LT.IYSTART) II=II+100
      IF(IY.LT.0.OR.II.LT.IY1.OR.II.GT.IY2)
     *   CALL QUITER(0.,IY,NERR,INERR)
      IY=II
      RETURN
      END
