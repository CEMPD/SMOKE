C 
C  Changes: (Last change first) 
C 
C  17 Dec 01 @EPA-BAG        Remove sawtooth 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C 25 Aug 00 AIR Task 08: Increased number of I/M programs from 5 to 7 
C           to allow for EVAP I/M. Added EVFBEN. 
C 17 Feb 99 @DynTel-RJD 2-663  Removed WAVPCT AND added WAVPRE81 and 
C           WAV81PLUS. 
C 
C 28 Oct 98 @DynTel-MLA 2-663  Broke WAIVER into 2 parts, WAVPRE81 for 
C           pre-1981 vehicles and WAV81PLUS for 81+ vehicles. 
C 
C 28 Oct 98 @DynTel-MLA 2-663  Added a dimension to ICYIM. It now saves 
C           end years for I/M programs as well as start years. The maximum 
C           number of I/M programs was changed to 5 in all arrays. Added 
C           variables WAVPCT, EXEMPTAGE, and IMGRACE for Mobile6. 
C 
C The data is initialized in BD_IM1. 
C 
      INTEGER EXEMPTAGE,ICYIM,IMGRACE,ISTRIN,MODYR 
      REAL CRIM,WAVPRE81,WAV81PLUS,EVFBEN,OLDFRAC,NEWFRAC,BIENRATIO 
      COMMON /IMPAR1/ ICYIM(2,7),ISTRIN(7),MODYR(2,7),WAVPRE81(7), 
     &        WAV81PLUS(7),EXEMPTAGE(7),IMGRACE(7),CRIM(7),EVFBEN(4,7), 
     &        OLDFRAC,NEWFRAC,BIENRATIO(25,3) 
