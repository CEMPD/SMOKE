C 
C  TTCINF contains the variables that holds the values for the 
C         Technical Training and Certification Benefits(TTC)  
C 
C  Changes: (Last change first) 
C 
C  15 Sep 00 AIR Task 08: Increased I/M related arrays from 5 to 7. 
C  19 Nov 98 @DynTel-MLA 2-663  Removed the IHYB variable since Retest 
C            Based Hybrid programs are not supported in Mobile6. Changed 
C            the dimensions of IASM, ITTCFG, TTCDAT, and TTCEFF to allow 
C            TTC information to be stored separately for each I/M program. 
C 
C  The data is initialized in BD38 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  TTCDAT    R    holds the values that have to be for the current run. 
C  TTCDEF    R    holds the default values of the TTC credits. 
C  ITTCFG    I    TTC flag is set by user and indicate: 
C                          =1 "NO" for TTC program, 
C                          =2 "YES" for TTC program, 
C                          =3 One Time Data Section contains "TTC record" 
C  Subroutine GETIMP can sets the additional values of ITTCFG flag to indicate: 
C                          =4 To use the user's percent of the 
C                             TTC effectiveness 
C                          =5 To read the TTC.D file and to use the user's 
C                             percent of the TTC effectiveness              
C  TTCEFF    R    Percent of the TTC effects on I/M benefits 
C**IHYB      I    Flags that indicate Hybrid Program for I/M tests  
C  IASM      I    Flags that indicate ASM procedure  for I/M tests  
C  TTCSML    C    Smart line from TTC.D file  
C 
C 
C  Notes: 
      INTEGER ITTCFG,IASM 
      REAL    TTCDAT,TTCEFF,TTCDEF 
      CHARACTER*80 TTCSML 
      COMMON /TTCINF/ TTCDAT(24,3,18,7),ITTCFG(7),TTCEFF(7),IASM(7), 
     *                TTCDEF(24,3,18),TTCSML(7)   
