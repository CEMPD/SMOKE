C Version "@(#)$Id$ $Source$ $Date$ 

C
C  PARAME sets the upper bound of calendar year.  The upper bound of calendar
C  year should be within 100 years from lower bound of calendar year(NSYR).
C  NSYR and NEYR are used to initialize IYSTART and IYEND.  IYSTART and IYEND 
C  are used to pass values to procedures, write statement, etc.
C
C
C  Parameter Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -------------------------------------------------
C   NEYR     I    defines upper bound of calendar year
C
C  Notes:
C
C  24-May-1996 @ DynTel-yc Req. 2-621 PARAME was created to make it easier
C  to change the lower and upper bounds of calendar year.
C  3-Sep-1996 @DynTel-yc Req. 0-448 PARAME was divided into two parts to
C              avoid warning message.
C
      INTEGER NEYR
      PARAMETER (NEYR=2051)
