C Data is initialized in BD_IM1 and OBD_HI_FRCN 
C 
C 06 Oct 00 @EPA-BAG BUG155 Remove OBDPC 
C 15 Sep 00 AIR Task 08: Changed MIL to WMIL for easier variable tracking 
C 14 Jun 00 @DynTel-eg bug146, Added AR94B 
C 29 Apr 99 @DynTel-ddj 2-694, Added array R_S_CONV and removed 
C           non-used variables/arrays. 
C 09 Mar 99 @DynTel-RJD 2-694 New data for Mobile 6 
C 
      REAL HIGH3, HIGH4, WMIL, OBDF,  OBDIM, OBDREP, R_S_CONV, 
     &     AR94B 
      COMMON /LDGOBD/ HIGH3(25,3,5), HIGH4(25,3,5), R_S_CONV(11,3,2), 
     &                WMIL, OBDF(3), OBDIM, OBDREP, AR94B 
