C 
C  RESUL4 contains the variables that accumulate separate values for start 
C  and running emissions and an array that is used to hold the sum of the 
C  evaporative emissions. 
C 
C  Changes: (Last change first) 
C 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C 
C /RESUL4/ 
C 
C   IRUNEF   I    IRUNEF            Token for running emissions in EFIVT 
C   ISTARTEF I    ISTARTEF          Token for start emissions in EFIVT 
C   ITEVAP   I    ITEVAP            Token for total evap. emiss. in EFIVT 
C   RUNEF    R    RUNEF(IP, IV)     Running emissions for each IP, IV 
C   STARTEF  R    STARTEF(IP, IV)   Start emissions for each IP, IV 
C   TOTEVAP  R    TOTEVAP(IV)       Sum of HC evap. emissions for each IV 
C   VRUN     R    VRUN(IP)          Emissions weighted over veh. classes 
C   VSTART   R    VSTART(IP)        Emissions weighted over veh. classes 
C 
C The data is initialized in BD27. 
C 
      INTEGER  IRUNEF 
      INTEGER  ISTARTEF 
      INTEGER  ITEVAP 
      REAL     RUNEF(3,29) 
      REAL     STARTEF(3,29) 
      REAL     TOTEVAP(29) 
      REAL     VRUN(3) 
      REAL     VSTART(3) 
      REAL     VTEVAP 
C 
      COMMON /RESUL4/ RUNEF,STARTEF,TOTEVAP,IRUNEF,ISTARTEF,ITEVAP, 
     &                VRUN,VSTART,VTEVAP 
