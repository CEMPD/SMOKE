C 
C     The data is initialized in BD17. 
C  
C  Changes: (Last change first) 
C 
C 
C 08 May 00 AIR Task 03: Added TEMEVP and TDIFF (from MOBILE5b) 
C  8 Apr 99 @DynTel-bag 2-000 Removed ambient temperature, AMBT. 
C 11 Mar 99 @DynTel-bag 2-664  Added TEMP_DIST for output of type of 
C           user supplied temperatures. Deleted TEMEVP array. 
C 30 Oct 98 @DynTel-bag 2-677  Added DEL_TEM(24) and GTDEL_TEM(24) 
C  6 Aug 98 @DynTel-bag 2-665  Added TEMHRLY(24) and GLTHRLY(24)  
C           for Mobile6.  
C  6 May 98 @DynTel-MLA 2-659  Added GBLTMIN and GBLTMAX for Mobile6.  
C 
      LOGICAL TEMP_DIST 
      REAL TEMMIN, TEMMAX, TEMEXH, TEMEVP, TEMHRLY, DEL_TEM 
      REAL GLTMIN,GLTMAX,GLTHRLY,GLDEL_TEM 
      COMMON /TEMPS/ TEMMIN,TEMMAX,TEMEXH,TEMEVP,TEMHRLY(24), 
     &               DEL_TEM(24),GLTMIN,GLTMAX,GLTHRLY(24), 
     &               GLDEL_TEM(24),TEMP_DIST 
