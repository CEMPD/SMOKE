      SUBROUTINE PMDRCPSOF(ICY) 
C 
C  Computes the diesel RCP and SOF from the 
C  difference of EMFAC_PM and SO4. (Note that 
C  EMFAC_PM = RCP + SOF + SO4.) 
C  RCP and SOF are dependent on particle size but 
C  are assumed to be independent of fuel sulfur content. 
C  SOF is called OCARBON and RCP is called ECARBON 
C  in the descriptive output. 
C 
C  Changes : 
C 
C  17 Apr 01 AIR Task 11: Changed code to handle expanded 
C            FVMT dimensions and prevent negative RCP+SOF 
C  02 Mar 01 AIR Task 11: New Routine 
C 
C  Called by PMCALX 
C 
C  Calls: None 
C 
C  Common Blocks: See Includes below. 
C 
C  Local variable / array dictionary: 
C 
C   Name     Type             Description 
C  ------    ----  ----------------------------------------------------- 
C  IDX         I   IDX=(MAXYRS+1)-JDX 
C  IH          I   Hour 
C  IROAD       I   Road Type 
C  IV          I   Vehicle Class 
C  IVPART      I   MOBILE6 vehicle class maped into PART6 vehicle class 
C  JDX         I   Age 
C  RCPSOF      R   RCP + SOF (EMFAC_PM - SO4) 
C  SUMRCPIH    R   Temporary summation of RCP 
C  SUMRCPJDX   R   Temporary summation of RCP 
C  SUMSOFIH    R   Temporary summation of SOF 
C  SUMSOFJDX   R   Temporary summation of SOF 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FRACT.I'        ! SOFPNT 
      INCLUDE 'IOUCOM.I'       ! IOUREP,IOUERR 
      INCLUDE 'IVTYPE.I'       ! MAXVEH, VDM6PRT, VVDSL 
      INCLUDE 'MAXIMA.I'       ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'       ! TF 
      INCLUDE 'PART3.I'        ! VPM, DBPM, EFPM 
      INCLUDE 'PSIZE.I'        ! PSDIE 
      INCLUDE 'SPEED9.I'       ! FVMT, HVMT 
      INCLUDE 'VMXCOM.I'       ! VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
C 
C  Declare local variables/arrays 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER IVPART 
      INTEGER JDX 
      INTEGER MY 
      REAL    RCPSOF 
      REAL    SUMRCPIH 
      REAL    SUMRCPJDX 
      REAL    SUMSOFIH 
      REAL    SUMSOFJDX 
C 
C 
C  Re-iInitialize to hold final values. 
C 
      VPM(2)=0.  !SOF 
      VPM(3)=0.  !RCP 
C 
C  RCP and SOF can only be calculated for diesel vehicles. 
C 
      DO IV=1,MAXVEH 
C 
       IF(VVDSL(IV).EQ.0) CYCLE 
C 
C  Re-initialize to hold final values. 
C 
       EFPM(2,IV)=0.  !SOF 
       EFPM(3,IV)=0.  !RCP 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Map the MOBILE6 vehicle classes into PART5 diesel space. 
C 
       IVPART = VDM6PRT(IV) 
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
C  Initialize 
C 
        SUMRCPJDX=0. 
        SUMSOFJDX=0. 
C 
C  Loop for hour 
C 
        DO IH=1,MAXIH 
C 
         SUMRCPIH=0. 
         SUMSOFIH=0. 
C 
C  Loop for road type 
C 
         DO IROAD=1,4 
C 
C  Compute SOF+RCP. The SO4 being subtracted here is from PMDSO4SO2 and 
C  is based on the BER associated sulfur level, not the in use level. 
C 
          RCPSOF = DBPM(10,IROAD,IH,JDX,IV) - DBPM(1,IROAD,IH,JDX,IV) 
C 
C  Produce an error message when negative RCPSOF values occur.  
C  Rest such values to zero. 
C 
        IF(RCPSOF.LT.0.) THEN 
CDFK        WRITE(IOUOUT,10) RCPSOF,ICY,MY,IV 
CDFK        WRITE(IOUERR,10) RCPSOF,ICY,MY,IV 
CDFK  10    FORMAT('*** Error:  a negative Diesel RCPSOF value of ',F8.3,/ 
CDFK     &         '    occurred for CY=',I4,' MY=,',I4,' and IV=',I2) 
        RCPSOF=0.0 
        END IF 
C 
C  Adjust for particle size 
C 
          RCPSOF = RCPSOF * PSDIE 
C 
C  Compute SOF and RCP separately. Compute composites. 
C 
          DBPM(2,IROAD,IH,JDX,IV) = RCPSOF * SOFPNT(IVPART)          !SOF 
          DBPM(3,IROAD,IH,JDX,IV) = RCPSOF - DBPM(2,IROAD,IH,JDX,IV) !RCP 
C 
          SUMSOFIH = SUMSOFIH+DBPM(2,IROAD,IH,JDX,IV)*FVMT(IROAD,IH,IV) !SOF 
          SUMRCPIH = SUMRCPIH+DBPM(3,IROAD,IH,JDX,IV)*FVMT(IROAD,IH,IV) !RCP  
C 
         END DO  !IROAD 
C 
         SUMSOFJDX = SUMSOFJDX + SUMSOFIH * HVMT(IH) !SOF 
         SUMRCPJDX = SUMRCPJDX + SUMRCPIH * HVMT(IH) !RCP 
C 
        END DO  !IH 
C 
        EFPM(2,IV)= EFPM(2,IV) + SUMSOFJDX * TF(IDX,IV) !SOF 
        EFPM(3,IV)= EFPM(3,IV) + SUMRCPJDX * TF(IDX,IV) !RCP 
C 
       END DO !JDX 
C 
       VPM(2) = VPM(2) + EFPM(2,IV) * VMTMIX(IV) !SOF 
       VPM(3) = VPM(3) + EFPM(3,IV) * VMTMIX(IV) !RCP 
C 
      END DO !IV 
C 
      RETURN 
      END 
