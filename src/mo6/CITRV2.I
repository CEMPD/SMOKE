C 
C  Changes: (Last change first) 
C 
C  12 Jun 00 @EPA-djb Removed RVPX and added RVPNOX and OXYRVP. 
C  17 May 00 AIR Task 03: Replace RVPHS array with single variable. 
C  03 Dec 98 @DynTel-BAG 2-664 Added hourly arrays for hot soak and exhaust RVP's. 
C            Removed RVPLOS. 
C 
C The data is initialized in BD30. 
C 
      REAL RVPHS,RVP090,RVP115,RVPNOX,OXYRVP 
      COMMON /CITRV2/ RVPHS,RVP090,RVP115,RVPNOX,OXYRVP 
