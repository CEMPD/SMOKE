      FUNCTION HDCPIC(MY,IV) 
C 
C  HDCPIC picks the HDV gpbhp to gpm units conversion factor. 
C 
C  Called by HDCGPM. 
C 
C  Changes: (Last change first) 
C 
C  26 Sep 00 @EPA-EG      Added HDCFAC Flag Structure. 
C  28 Apr 00 @EPA-elg Added ligt duty gas and light duty diesel vehicle 
C             assignment for HDCPIC 
C   05 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C   7 Oct 98  @DynTel-HQ 2-678 Comment statements and VVHDDV array index 
C  28 Sep 98  @DynTel-HQ 2-678 Conversion Factor 
C 
C  Input on call: 
C 
C    argument list: MY,IV 
C    common blocks: 
C    /HDCCOM/ HDCFAC,LOWYR,MAXCIG 
C    /VVNAMS/ VVHDDV, VVHDGV  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: HDCPIC 
C 
C  Notes: 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'HDCCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INTEGER IG, MY, IV, IVHD 
      REAL HDCPIC 
C 
C  There is a different conversion factor for each model year from 
C  LOWYR up through MAXCIG's year.  Below LOWYR, use LOWYR's. 
C  Above MAXCIG's year, use MAXCIG's year's. 
C 
      IG=MY-LOWYR+1 
      IF(IG.LT.1) IG=1 
      IF(IG.GT.MAXCIG) IG=MAXCIG 
C 
C 
      IF(IV.LE.5 .OR. IV.EQ.14 .OR. IV.EQ.15 .OR. 
     *   IV.EQ.24 .OR. IV.EQ.28) THEN 
C 
         HDCPIC = 1.00 
C 
      ELSE 
C 
C   
C     IV = 6-13 for Heavy Duty GAS Vehicle Type 
C     IV = 16-23 for Heavy Duty Diesel Vehicle Type 
C     IV = 25 for gas bus 
C     IV = 26 for commercial diesel bus 
C     IV = 27 for school diesel bus 
C     IV = 28 for transit diesel bus 
C 
        IF(VVHDGV(IV).EQ.1) IVHD=IV-5 
        IF(VVHDDV(IV).EQ.1) IVHD=IV-7 
        IF(IV.EQ.VTGASBUS)  IVHD=17 
        IF(IV.EQ.VTURBAN)   IVHD=18 
        IF(IV.EQ.VTSCHOOL)  IVHD=19 
C 
C 
C     IVHD=1-8 for VVHDGV 
C     IVHD=9-16 for VVHDDV 
C     IVHD=17 for GAS BUS 
C     IVHD=18-19 for transit/commercial and school diesel bus 
C 
        IF(HDCFAC_FLG.EQ.0) THEN 
          HDCPIC=HDCFAC(IG,IVHD) 
        ELSEIF(HDCFAC_FLG.EQ.1) THEN 
          HDCPIC=HDCFAC2(IG,IVHD) 
        ENDIF 
C 
      ENDIF 
C 
      RETURN 
      END 
