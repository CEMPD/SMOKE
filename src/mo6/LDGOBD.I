C Data is initialized in BD_IM1 and OBD_HI_FRCN
C
C 23 Sep 02 AIR Task 35B: Added OBDIMFLG to prevent using OBD EXH RESPONSE
C           whenever an I/M program is specified.
C 24 Apr 24 AIR Task 32: Changed OBDF to OBDXF, OBDXM, OBDVF, OBDVA
C           for OBD EXHAUST and EVAP RESPONSE, added OBDXFLG, OBDVFLG
C           OBDXLIN, OBDVLIN
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C 06 Oct 00 @EPA-BAG BUG155 Remove OBDPC
C 15 Sep 00 AIR Task 08: Changed MIL to WMIL for easier variable tracking
C 14 Jun 00 @DynTel-eg bug146, Added AR94B
C 29 Apr 99 @DynTel-ddj 2-694, Added array R_S_CONV and removed
C           non-used variables/arrays.
C 09 Mar 99 @DynTel-RJD 2-694 New data for Mobile 6
C
      REAL OBDHINIM, OBDHI_IM, WMIL, OBDXF, OBDXM, OBDVF, OBDIM,
     &     OBDREP, R_S_CONV, AR94B
      INTEGER OBDVA
      LOGICAL OBDXFLG, OBDVFLG, OBDXLIN, OBDVLIN, OBDIMFLG
      COMMON /LDGOBD/ OBDHINIM(25,3,5), OBDHI_IM(25,3,5),
     &     R_S_CONV(11,3,2), WMIL, OBDXF(3), OBDXM(2), OBDVF(3),
     &     OBDVA(2), OBDIM, OBDREP, AR94B, OBDXFLG, OBDVFLG,
     &     OBDXLIN, OBDVLIN, OBDIMFLG
