C Version "@(#)$Id$ $Source$ $Date$ 

C
C  PARAMS sets the lower bound of calendar year.  The lower bound of calendar
C  year should be within 100 years from upper bound of calendar year(NEYR).
C  NSYR and NEYR are used to initialize IYSTART and IYEND.  IYSTART and IYEND 
C  are used to pass values to procedures, write statement, etc.
C
C
C  Parameter Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -------------------------------------------------
C   NSYR     I    defines lower bound of calendar year
C
C  Notes:
C
C  3-Sep-1996 @DynTel-yc Req. 0-448 PARAMS was splitted from PARAME to avoid
C             warning message about unreferenced symbol for NSYR. 
C
      INTEGER NSYR
      PARAMETER (NSYR=1952)
