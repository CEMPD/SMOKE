C 
C The data is initialized in BD30. 
C  
C  Changes: (Last change first) 
C 
C 06 Oct 00 @EPA-BAG  bug155 Removed INDUSTRY AVERAGE FUEL parameters IAFCNT 
C           and GLIAFCNT 
C 26 Jun 00 @EPA-elg  Changed IPHASE to JPHASE. 
C 20 Apr 00 @EPA-ddj E-7  Added FUEL_P. 
C 25 Feb 00 @DynTel-ddj 1-028, Removed obsolete ASTM code. 
C  4 Apr 00 AIR Task 02: Removed IASTM, ASTMCL, and GLASTM from list. 
C  6 May 98 @DynTel-MLA 2-659  Added GBLASTM for Mobile6.  
C 
C 
      LOGICAL RFGON 
C 
      INTEGER GLRFG 
      INTEGER GLWIN 
      INTEGER FUEL_P 
      INTEGER JPHASE 
      INTEGER RFGFLG 
      INTEGER RFG_REGION 
      INTEGER SEAFLG 
      INTEGER WINFLG 
C 
      COMMON /RFORM1/ JPHASE, RFGFLG, WINFLG, SEAFLG, 
     &                GLRFG,  GLWIN,  RFGON,  RFG_REGION, 
     &                FUEL_P 
