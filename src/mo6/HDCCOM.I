C 
C  Changes: (Last change first) 
C 
C  26 Sep 00 @EPA-ELG   Added HDCFAC2 array. 
C  29 Sep 98 @DynTel-HQ 2-678  Add CONVFLG for conversion factor  
C 
C The data is initialized in BD29. 
C 
      REAL HDCFAC,HDCFAC2 
      INTEGER LOWYR,MAXCIG,HDCFAC_FLG 
      COMMON /HDCCOM/ HDCFAC(56,19),HDCFAC2(56,19),LOWYR,MAXCIG, 
     *                HDCFAC_FLG 
