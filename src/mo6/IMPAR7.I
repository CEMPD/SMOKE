C 
C  IMPAR7 contains several I/M program variables. 
C 
C  Changes: (Last change first) 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  IREAD     I    Shows whether prev. I/M credits were high/low altitude 
C  NAMIMC    C    Name of alternate Tech4+ I/M credit file. 
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
C  1 Apr 98 @DynTel-MLA 2-659  Added NAMIMC for the transitonal labeled  
C           input version of Mobile6. 
C 
C 13 Apr 98 @DynTel-BAG 2-661 Added NAMIMC and increased the length of  
C           TNAME to 80 characters 
C  The data is initialized in BD_IM2 
C 
      CHARACTER*80 TNAME 
      CHARACTER*80 NAMIMC  
      INTEGER IREAD, SVALTC 
      COMMON /IMPAR7/ TNAME, NAMIMC, IREAD, SVALTC                                                                    
