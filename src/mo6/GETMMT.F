      SUBROUTINE GETMMT(INPREC,NSCEN,INERR,RC,TEMP_TYP)
C
C  GETMMT reads in the minimum and maximum temperatures and either reads or
C  determines hourly temperatures for a scenario.
C
C  Called by DSPTCH.
C
C  Calls QUITER.
C
C  Changes : 
C
C  25 Jul 02 AIR Task Bug 393: fixed inconsistent temperature limits and
C            corrected hourly temperature error message.
C  24 Jul 01 @EPA-bag remove code for rounding off temperatures
C  31 Oct 00 @EPA-bag bug281  Input hourly temperatures with a Free Format
C                        READ statement.
C  25 May 00 @EPA-elg    Changed input from Fixed Format to Free Format.
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations      
C  16 Mar 99 @DynTel-BAG 2-664  Added logical variable, TEMP_DIST, for output of type
C            of user supplied temperatures.
C  30 Oct 98 @DynTel-BAG 2-677  Added hourly temperature determinations; uses
C           variable TEMP_TYP to select between temperature profile 
C           determinations. 
C       
C  Input on call:
C
C    argument list: TEMP_TYP, INPREC, NSCEN, INERR
C
C    common blocks:
C
C    /RVPNAT/  RVPLIM
C
C  Output on return:
C
C    argument list: INERR, RC
C
C    common block:
C    /TEMPS/   GLTMIN, GLTMAX, GLTHRLY, GLDEL_TEM, TEMMIN, TEMMAX,
C              TEMHRLY, DEL_TEM, TEMP_DIST
C    /RVPNAT/  RVPDIU
C
C  Local variable / array dictionary:
C
C  Name                 Type              Description
C  ------               ----  ----------------------------------------------------
C  STNDRD_DEL_TEMP(IH)   R   Hourly temperature differences determined from the 
C                             Code of Federal Regulations' daily temperature profile
C  IH                    I   Hour index
C  SCAL_FCTR             R   Scaling factor to determine hourly temperatures based
C                             on TEMMIN, TEMMAX, and STNDRD_DEL_TEMP
C
C  Notes:
C
C  GETMMT is new for Mobile6. This code used to reside in RDCITY.
C
      IMPLICIT NONE
      INCLUDE 'RVPNAT.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'IOUCOM.I'
C
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER, INTENT(OUT) :: RC
      INTEGER, INTENT(IN) :: NSCEN
      CHARACTER INPREC*(*)
      INTEGER   TEMP_TYP, IH
      REAL      SCAL_FCTR
      REAL, DIMENSION(24) :: STNDRD_DEL_TMP=(/
C           Hourly temperature increases/decreases determined from the
C           Code of Federal Regulations' standardized temperature
C           profile (40CFR86)
     &       0.5, 3.0, 4.8, 4.9, 4.2, 3.7, 2.0, 0.7,
     &       0.2,-0.5,-1.4,-2.4,-3.1,-3.1,-2.7,-1.9,
     &      -1.9,-1.8,-1.4,-1.1,-0.8,-0.6,-0.7,-0.6 /)
C
      RC = 0
C
C Use the variable TEMP_TYP to determine the method of filling
C the hourly temperature array TEMHRLY(HR); this will be done 
C irrespective of whether hourly data is requested
C
      SELECT CASE(TEMP_TYP)
C
C Hourly temeperatures are determined from the min. and max. daily temperature 
C and the standardized daily temperature profile
C
      CASE(0)
C
C  Read and verify the min. and max. daily temperatures.
C
        TEMP_DIST = .FALSE.
C
        READ(INPREC(21:),*,IOSTAT=RC,ERR=96)
     &        TEMMIN, TEMMAX
C
C  Check temperature ranges.  These bounds are for catching data entry errors
C  only.  What constitutes reasonable min & max temperature values depends on
C  local conditions, especially AMBT.
C
        IF(TEMMIN.LT.0.0.OR.TEMMIN.GT.100.0)
     *      CALL QUITER(TEMMIN,1,91,INERR)
C
        IF(TEMMAX.LT.0.0.OR.TEMMAX.GT.120.0)
     *      CALL QUITER(TEMMAX,2,91,INERR)
C
        IF(TEMMAX.LT.TEMMIN) CALL QUITER(TEMMAX,0,92,INERR)
C
C Determine the hourly temperatures from TEMMIN, TEMMAX, and STNDRD_DEL_TMP
C and round the temperatures to the nearest 0.1F.
C
        TEMHRLY(1) = TEMMIN
        SCAL_FCTR = (TEMMAX - TEMMIN)/24.
        DO IH = 2,24
          DEL_TEM(IH) = SCAL_FCTR*STNDRD_DEL_TMP(IH-1)
          TEMHRLY(IH) = TEMHRLY(IH - 1) + DEL_TEM(IH)
        END DO
C
C User input of hourly temperatures
C
      CASE(1)
C
C  Read and verify the hourly temperatures.
C
        TEMP_DIST = .TRUE.
        READ(INPREC(21:),*,ERR=95,END=96) (TEMHRLY(IH), IH = 1,12)
        READ(IOUGEN,*,ERR=95,END=96) (TEMHRLY(IH), IH = 13,24)
C
C Find the maximum and minimum temperatures. Also, check temperature
C ranges. These bounds are for catching data entry errors only.
C
        TEMMIN = TEMHRLY(1)
        TEMMAX = TEMHRLY(1)
        DEL_TEM(1) = 0.
C
         IF(TEMHRLY(1).LT.0.0.OR.TEMHRLY(1).GT.120.0)
     *      CALL QUITER(TEMHRLY(1),1,40,INERR)
C
        DO IH = 2,24
C
         IF(TEMHRLY(IH).LT.0.0.OR.TEMHRLY(IH).GT.120.0)
     *      CALL QUITER(TEMHRLY(IH),IH,40,INERR)
C
          IF( TEMHRLY(IH) .GT. TEMMAX ) TEMMAX = TEMHRLY(IH)
          IF( TEMHRLY(IH) .LT. TEMMIN ) TEMMIN = TEMHRLY(IH)
          DEL_TEM(IH) = TEMHRLY(IH) - TEMHRLY(IH-1) 
        END DO
C
      END SELECT
C
C  Compute diurnal evap RVP limits.  The upper bound is dependent on 
C  the max and min temperatures.
C
      IF(TEMMIN.LE.95.) RVPDIU(2)=26.823-0.004781*TEMMIN-0.13407*TEMMAX
      IF(TEMMIN.GT.95.) RVPDIU(2)=26.000-0.004781*TEMMIN-0.13407*TEMMAX
      IF(RVPDIU(2).GT.RVPLIM(2)) RVPDIU(2)=RVPLIM(2)
C
      IF (NSCEN.EQ.0) THEN
        GLTMIN = TEMMIN
        GLTMAX = TEMMAX
        DO IH = 1,24
          GLDEL_TEM(IH) = DEL_TEM(IH)
          GLTHRLY(IH) = TEMHRLY(IH)
        END DO
      ENDIF
C
      GOTO 99
C
C  Error or EOF on any attempted read.
C
   95 RC = 1
      INERR = INERR + 1
      GOTO 99
   96 RC = -1
      INERR = INERR + 1
      GOTO 99
C
   99 RETURN
      END
