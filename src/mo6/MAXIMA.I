C 
C The data is initialized in BD16. 
C 
C Changes : 
C 
C Apr 13 99 @DynTel-bag 2-000 Removed MAXIH_EVAP; merged variables for 
C           exhaust and evaporative hourly loop. 
C Oct 30 98 @DynTel-bag 2-677 Added MAXIH_EVAP evaporative hourly loop. 
C Sep 16 98 @DynTel-MLA 2-659 Added MAXIH and MAXEMIT for hourly loop. 
C Sep 98 @DynTel-RJD 0-632 Added new MCYRS variable that specifies the 
C            maximum model year index for motorcycles included in basic 
C            emission factor calculations. 
C 
      INTEGER MAXPOL,MAXREG,MAXYRS,MCYRS,MAXIH,MAXEMIT 
      COMMON /MAXIMA/ MAXPOL,MAXREG,MAXYRS,MCYRS,MAXIH,MAXEMIT 
 
