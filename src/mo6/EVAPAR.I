C 
C  Changes: (Last change first) 
C 
C  27 Feb 01 @EPA-bag bug314  removed hourly diurnal parameters; they are now in 
C            DIRNLDAT.FOR 
C  30 Oct 98 @DynTel-bag 2-677  Added parameters for hourly diurnal calculations. 
C 
C The data is initialized in BD11. 
C  
      INTEGER IFDS,ISTDC,ISTDT,ISTDP 
      REAL FINJ 
      COMMON /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ(2) 
