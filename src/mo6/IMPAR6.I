C 
C  IMPAR6 contains several I/M program variables. 
C 
C  Changes: (Last change first) 
C 
C 25 Aug 00 AIR Task 08  Expanded arrays IFREQ and INTYP to 
C           accomodate 7 IM programs. 
C 17 Feb 99 @DynTel-MLA 2-663 Expanded arrays IFREQ and INTYP to hold 
C           information about each I/M program. 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  IFILET    I    Flag - user-supplied or default I/M cred. file names 
C  IFREQ     I    I/M program inspection frequency - 1=annual, 2=bienn. 
C  IMEFFR    I    Program Effectiveness flag from I/M prog. desc. rec. 
C  INTYP     I    I/M program type - 1=test-and-repair, 2=test-only 
C 
C 
C  Notes: 
C 
C 8 Feb 96 @DynTel-MLA 2-611 Added IMEFFR, a flag to show whether the 
C    user has entered a program effectiveness record to supply alternate  
C    effectiveness rates for test-and-repair I/M, ATP, and functional  
C    pressure and purge tests. 
C 
C The data is initialized in BD_IM2 
C 
      INTEGER IFREQ,INTYP,IFILET,IMEFFR 
      COMMON /IMPAR6/ IFREQ(7),INTYP(7),IFILET,IMEFFR 
