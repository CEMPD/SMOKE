C 
C  RSDINF contains the variables that holds the value for the 
C         Remote Sensing Device (RSD) program  
C 
C  Changes: (Last change first) 
C 
C  15 Sep 00 AIR Task 08: Increase arrays from 5 to 7 to allow for Evap I/M. 
C  17 Feb 99 @DynTel-MLA 2-663  Replaced RSDFLG with RSDCLEAN and RSDDIRTY. 
C 
C  18 Nov 98 @DynTel-MLA 2-663  Changed RSDFLG from a single integer to 
C            an array to store separate values for each I/M program. 
C         
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  RSDFLG    I    Flag that indicate RSD  Program   
C 
C  Notes: 
C 
C  The data is initialized in GETIMC.FOR 
C 
      INTEGER RSDCLEAN,RSDDIRTY 
      COMMON /RSDINF/ RSDCLEAN(7),RSDDIRTY(7) 
