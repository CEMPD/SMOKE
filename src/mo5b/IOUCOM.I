C Version "@(#)$Id$ $Source$ $Date$ 

C
C  IOUCOM contains the variables that hold the values of the unit number
C   for INPUT and OUTPUT files.
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C
C   IOUIMD   I    unit number that connected with I/M credit file
C   IOUGEN   I    unit number that connected with MOBILE5 inputs
C   IOUREP   I    unit number that connected with report  file
C   IOUERR   I    unit number that connected with error   file
C   IOUASK   I    unit number that connected with console
C   IOUIMD   I    unit number that connected with PPEFF.D file
C   IOUTTC   I    unit number that connected with TTC.D file
C
C  Notes:
C
      INTEGER IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK,IOUPFF,IOUTTC 
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK,IOUPFF,IOUTTC
