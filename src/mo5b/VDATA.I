C Version "@(#)$Id$ $Source$ $Date$ 

C
C  VDATA contains the variable that holds the value of the name of the version
C  of the model.
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C  VERSON    C    holds the value of the name of the version of the model.
C
C
C  Notes:
C
C  12-Dec-1994 @ CSC-rwk Req. 2-000 VDATA was created to make it easier to
C  change the name of the version of the model.  This was done for quality
C  control.
C
      CHARACTER*80 VERSON
      COMMON /VDATA/ VERSON
