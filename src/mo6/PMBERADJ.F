      SUBROUTINE PMBERADJ(ICY,INERR) 
C 
C  This routine first computes the total gasoline PM BER by summing the  
C  component parts (carbon and SO4). (NOTE: The SO4 here is at the 
C  associated BER sulfur levle, *NOT* the inuse fuel level. NOTE: The 
C  GASPM value is at PM30, so that the Total PM is at the regulation PM30 
C  level.) If the total BER is above the regulation BER, then the carbon  
C  component adjustment factor is computed. 
C 
C  Changes : 
C 
C  14 Jan 02 AIR Task 25: Properly directed error output to IOUERR and 
C                         corrected error message formatting. 
C  04 Jan 02 AIT Task 25: Updated negative adjustment error handling. 
C  20 Oct 01 AIR Task 25: New Module 
C 
C     Called by:  PMCALX 
C 
C  Input on call: 
C 
C    Parameter list: ICY,INERR 
C 
C    Common blocks:  See Includes below 
C 
C  Output on return:  
C    /PART3.I/ EFPM,VPM,DBPM 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I'   ! IOUERR 
      INCLUDE 'IVTYPE.I'   ! MAXVEH,VVGASMC 
      INCLUDE 'MAXIMA.I'   ! MAXYRS, MAXIH 
      INCLUDE 'MYRCAL.I'   ! TF 
      INCLUDE 'PART1.I'    ! PRTCHK, D_PPM, BERSULF, BERMAX 
      INCLUDE 'PART3.I'    ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'   ! FVMT, HVMT 
      INCLUDE 'VMXCOM.I'   ! VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER, INTENT(IN)    :: ICY     ! Calendar Year 
      INTEGER, INTENT(INOUT) :: INERR   ! Error 
C 
C  Declare local variables/arrays 
C 
      INTEGER IDX        ! IDX=(MAXYRS+1)-JDX 
      INTEGER IH         ! Hour index 
      INTEGER IROAD      ! Road index 
      INTEGER IV         ! Vehicle Class 
      INTEGER JDX        ! Age 
      INTEGER MY         ! Model year 
      REAL    PMBER      ! Basic emmission rates, g/mi 
C 
C 
C  Initialize 
C 
      BERADJ=1. ! Neutral adjustment factor 
       
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if not a gasoline vehicle 
C 
      IF(VVGASMC(IV).EQ.0) CYCLE 
C 
C  Skip if VMTMIX is 0. 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Loop for all ages 
C 
      DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
      IDX=(MAXYRS+1)-JDX 
C 
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
      IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
      MY=ICY-JDX+1 
C 
C  Loop for each hour 
C 
      DO IH = 1,MAXIH 
C 
C  Loop for each road type 
C 
       DO IROAD = 1,4 
C 
C  Build the Total PM here. It is the sum of the carbon and SO4, 
C  with the SO4 computed at the associated sulfur level. 
C 
         PMBER = DBPM(4,IROAD,IH,JDX,IV) ! Carbon at PM30 
     *         + DBPM(1,IROAD,IH,JDX,IV) ! SO4, at associated sulfur level 
C 
C  Compute the gasoline BER limit adjustment. Limits are computed by the 
C  PMBERPAR subroutine. If BERMAX is negative, no limit is applied. 
C  Warn if adjustment becomes negative. If so, reset to zero.  
C 
         IF(BERMAX(JDX,IV).GE.0. .AND. PMBER.GT.BERMAX(JDX,IV)) THEN 
C 
         BERADJ(IROAD,IH,JDX,IV)=(BERMAX(JDX,IV)           ! Adj=(Limit-S04)/GASPM 
     *                          -DBPM(1,IROAD,IH,JDX,IV)) 
     *                          /DBPM(4,IROAD,IH,JDX,IV) 
C 
         IF(BERADJ(IROAD,IH,JDX,IV).LT.0.) THEN 
         WRITE(IOUERR,10) BERADJ(IROAD,IH,JDX,IV),IROAD,IH,MY,IV 
   10    FORMAT(' Error: The gasoline Carbon PM adjustment is',F8.4, 
     *          ' for IROAD=',I2,', IH=',I2,' MY=',I4,' IV=',I2/ 
     *          ' This is likely caused by a mismatch between the'/ 
     *          ' gasoline sulfur level, the carbon ZML input via'/ 
     *          ' the file PMGZML.CSV, and the regulatory maximum'/ 
     *          ' PM emission rate.') 
         INERR=INERR+1 
         RETURN 
         END IF 
C 
         END IF 
C   
         END DO  !IROAD 
C 
        END DO  !IH 
C 
       END DO  !JDX 
C 
      END DO  !IV 
C 
      RETURN 
      END 
