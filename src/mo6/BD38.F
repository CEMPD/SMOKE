      BLOCK DATA BD38 
C 
C  BLOCK DATA Subprogram 38: /EVAL/ 
C 
C  BLOCK DATA 38 initializes the storage arrays used for by-month calculations. 
C 
C  Changes: (Last change first) 
C 
C 
C  Common block array subscripts: 
C 
C  MONTH(12)         -  MONTH  ( MEVAL ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------- 
C  /EVAL/: 
C 
C  MEVAL    I    month of evaluation for by-month output, 1=Jan, 2=Feb, etc. 
C  MONTH    C*4  4 character abbreviations of the months 
C  TTCDAT   R    Holds the values that have to be for the current run 
C  TTCDEF   R    Holds the default values of the TTC credits. 
C  ITTCFG   I    TTC (I&M technician training) flag is set by user and indicate: 
C                         =1 "NO" for TTC Program, 
C                         =2 "YES" for TTC Program, 
C                         =3 One Time Data Section contains "TTC reco 
C  Subroutine GETIMP can sets the additional values of ITTCFG flag to 
C                         =4 To use the user's percent of the  
C                            TTC effectiveness 
C                         =5 To read the TTC.D file and to use the us 
C                            percent of the TTC effectiveness 
C  TTCEFF   R    Percent of the TTC effects on I/M benefits 
C  IHYB     I    Flags that indicate Hybrid Program for I/M tests 
C  IASM     I    Flags that indicate ASM procedure for I/M tests 
C  TTCSML   C    Smart line from TTC.D file 
C 
C                                     
C  Notes: 
C 
C  BD 38 was added for MOBILE5v1. 
C 
      IMPLICIT NONE 
      INCLUDE 'EVAL.I' 
C 
C  /EVAL/: 
C 
      DATA MEVAL/1/ 
      DATA MONTH/'Jan.','Feb.','Mar.','Apr.','May.','June', 
     *           'July','Aug.','Sep.','Oct.','Nov.','Dec.'/ 
C 
      END 
