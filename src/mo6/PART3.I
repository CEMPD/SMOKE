C 
C  PART3 contains the PM results. 
C 
C  Changes: (Last change first) 
C 
C  09 May 01: AIR Task 11/5: Expanded arrays to handle PM idle. 
C  23 Feb 01: AIR Task 11: New INCLUDE. 
C 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C /PART3/ 
C   DPPM    R     PM Database Emissions 
C   EFPM    R     PM Emissions weighted for each veh. class 
C   VPM     R     PM Emissions weighted over all vehicle classes 
C 
C The data is initialized in BD43. 
C 
      REAL     BERADJ 
      REAL     DBPM 
      REAL     EFPM 
      REAL     VPM 
C 
      COMMON /PART3/ EFPM(11,28),VPM(11),DBPM(10,4,24,25,28), 
     *               BERADJ(4,24,25,28) 
