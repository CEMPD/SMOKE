C Version "@(#)$Id$ $Source$ $Date$ 

C
C  IMPAR7 contains several I/M program variables.
C
C  Common Block Dictionary:
C
C   Name    Type                  Description
C  ------   ----  -----------------------------------------------------
C  IREAD     I    Shows whether prev. I/M credits were high/low altitude
C  SVALTC    I    Last I/M credits read: 0=None, 1=Standard, 2=Alternate
C  TNAME     C    Name of alternate Tech1&2 I/M credit file.
C
C
C  Notes:
C
C  5 Jun 96 @DynTel-MLA 2-620 Added SVALTC, a flag used by GETSC1. If it is
C    zero, no I/M credits have been read. If it is 1, the I/M credit arrays
C    hold the default I/M credits. If it is 2, the arrays hold alternate 
C    I/M credits.
C
      CHARACTER*40 TNAME
      INTEGER IREAD, SVALTC
      COMMON /IMPAR7/ TNAME, IREAD, SVALTC                                                                   
