      SUBROUTINE  OBD_HI_FRCN(ICY,INERR)
C
C  OBD_HI_FRCN returns the proportion of high emitters amongst the 1994+
C  light duty gasoline vehicles and trucks equipped with OBD.(See equations
C  8 through 11 on page 15 of M6.EXH.009 or see eqns. 15-18 on pgs. 14&15
C  in M6.EXH.007.)
C
C  Called by EFCALX
C
C  Calls EF_HI_FRCN
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C  15 Sep 00 AIR Task 08: Replaced MIL with WMIL for easier variable tracking.
C            Removed improper initialization of OBDIM. OBDIM is now properly
C            intialized in BD_IM1 and REINIT.
C  15 Aug 00 @EPA-djb Bug 224, Changed logic for calculating high emitter cap.
C            removed NORML_FRCN.
C  26 Jul 00 @EPA-jk  Bug154 for default age (JDX) case replaced NORML_FRCN with
C                     1.0-OBDHINIM(JDX-1,IP,IV) to determine GROWTH fraction
C  14 Jun 00 @EPA-EG  Bug146 Revised assignments to OBD and
C                            OBDIM
C  08 Jun 00 @EPA-djb Bug113 Added ISTD to EF_HI_FRCN call.
C  23 May 00 @EPA-EG     Changed CASE(1) determination of FRCN;
C                        removed division by NORML_FRCN
C   7 Jan 00 @DynTel-ddj 1-024, Capped OBDHINIM and OBDHI_IM at 1.0.
C            Selection of OBD moved prior to the JDX loop.
C   8 Sep 99 @DynTel-BAG 2-660 Renamed HIEMIT to EF_HI_FRCN
C  30 Apr 99 @DynTel-ddj 2-694, Removed dead code.
C  02 Mar 99 @Dyntel-RJD 2-694 This function is new for Mobile 6.
C
C  Input on call or from calls:
C
C    argument list: ICY
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /LDGOBD/ WMIL,OBDXF,OBDXM,OBDXLIN,OBDIM
C    /MAXIMA/ MAXYRS
C
C  Output on return:
C
C    common blocks:
C    /LDGOBD/ OBDHINIM,OBDHI_IM
C
C  Local variable / array dictionary:
C
C   Name      Type              Description
C  ------     ----  -----------------------------------------------------
C  HI_YEAR_DIF  R     Change in fraction of high emitters over a period
C                     of 1 year
C  GROWTH       R     HI_YEAR_DIF divided by fraction of normals;
C                     the EPA documentation refers to this quantity as the
C                     high emissions vehicle growth rate
C  OBD          R     Fraction of motorists who respond to malfunction
C                     indicator light
C  VMTn         R     Vehicle miles travelled in units of 10000 miles at
C                     two (n=1,2) dates one year apart.
C  YEAR1        R     Fraction of high emitters in current model year.
C  YEAR2        R     Fraction of high emitters in the following model year.
C  A            R     Slope in linear OBD equation
C  B            R     Intercept in linear OBD equation
C
C  Notes:
C
C
      IMPLICIT NONE
      INCLUDE 'CUMCOM.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'LDGIMC.I'
      INCLUDE 'LDGOBD.I'
      INCLUDE 'MAXIMA.I'
C
C    Declare external functions
C
      REAL EF_HI_FRCN
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: ICY
      INTEGER, INTENT(INOUT) :: INERR
C
C     Declare local variables
C
      INTEGER IP
      INTEGER IV
      INTEGER JDX
      INTEGER MY
      INTEGER JMAX
      REAL    FRCN
      REAL    HI_YEAR_DIFF
      REAL    GROWTH
      REAL    OBD
      REAL    VMT1
      REAL    VMT2
      REAL    YEAR1
      REAL    YEAR2
      REAL    A
      REAL    B
C
      IF (ICY.LT.2018) THEN
        JMAX = ICY - 1993
      ELSE
        JMAX = 25
      END IF
C
      DO IV = VTLDGV, VTLDGT4
         DO IP = 1, 3
           DO JDX = 1, MAXYRS
              OBDHINIM(JDX,IP,IV) = 0.0
              OBDHI_IM(JDX,IP,IV) = 0.0
           END DO
C
C Calculate the fraction of high emitters that have been identified
C under OBD or OBD I/M programs on the emissions evaluation date.
C
           DO JDX = 1, JMAX
             MY = ICY - JDX + 1
C
C            mileage at age JDX
             VMT2 = CUMMIL(JDX,IV)/10000.0  !CUMMIL is in units of miles
C
C Assign the OBD Exhaust Response rate.
C If LINEAR option is selected (OBDXLIN=.TRUE.), then interpolate between
C 1st and 2nd responses (using Y=AX+B)
C
             IF (VMT2.LE.OBDXM(1)) THEN        ! 1st Inflection and Response
                OBD = OBDXF(1)             
             ELSE IF (VMT2.LE.OBDXM(2)) THEN   ! 2nd Inflection and Response
                OBD = OBDXF(2)
                IF(OBDXLIN) THEN               ! Linear fit (Y=AX+B) from 1st to 2nd Response
                  A=(OBDXF(2)-OBDXF(1))
     &             /(OBDXM(2)-OBDXM(1))
                  B=(OBDXM(2)*OBDXF(1)
     &              -OBDXF(2)*OBDXM(1))
     &             /(OBDXM(2)-OBDXM(1))
                  OBD=A*VMT2+B
                END IF          
             ELSE
                OBD = OBDXF(3)                 ! 3rd Response
             END IF
C
C
             SELECT CASE (JDX)
C
             CASE(1)  !special case with JDX=1 and HIGH(JDX-1=0,IP,IV)=0; see p.15 M6.EXH.007
C
               FRCN = EF_HI_FRCN(MY,IP,IV,VMT2,1,2,INERR)
               IF(FRCN .GT. 1.0) FRCN = 1.0
C
              OBDHINIM(JDX,IP,IV) = (1.0 - OBD*WMIL) *FRCN  !OBD without I&M program(see eq.15-18 pgs.14&15 M6.EXH.007)
              OBDHI_IM(JDX,IP,IV) = (1.0 - OBDIM*WMIL) *FRCN!OBD with I&M program(see eq.15-18 pgs.14&15 M6.EXH.007)
C
             CASE DEFAULT ! all vehicle ages except JDX=1
C
C              mileage at age JDX-1
               VMT1 = CUMMIL(JDX-1,IV)/10000.0
C              high emitter growth rate (see eqs.15-18 pgs.14&15 M6.EXH.007)
               YEAR1 = EF_HI_FRCN(MY,IP,IV,VMT2,1,2,INERR)
               IF(YEAR1 .GT. 1.0) YEAR1 = 1.0
               YEAR2 = EF_HI_FRCN(MY+1,IP,IV,VMT1,1,2,INERR)
               IF(YEAR2 .GT. 1.0) YEAR2 = 1.0
               HI_YEAR_DIFF = YEAR1 - YEAR2
               GROWTH = 0.0
               IF(OBDHINIM(JDX-1,IP,IV) .LT. 1.0) THEN
                 GROWTH = HI_YEAR_DIFF / (1.0-OBDHINIM(JDX-1,IP,IV))
               ENDIF
C
              OBDHINIM(JDX,IP,IV)=OBDHINIM(JDX-1,IP,IV)     !OBD without I&M program(see eq.11a p.15 M6.EXH.007)
     &                           +(1.0 - OBD*WMIL)
     &                           *(GROWTH*(1.0 - OBDHINIM(JDX-1,IP,IV)))
C
              OBDHI_IM(JDX,IP,IV) = OBDHI_IM(JDX-1,IP,IV)   !OBD with I&M program(see eq.11b p.15 M6.EXH.007)
     &                           + (1.0 - OBDIM*WMIL)
     &                           *(GROWTH*(1.0 - OBDHI_IM(JDX-1,IP,IV)))
C
             END SELECT
               OBDHINIM(JDX,IP,IV) = MIN(1.0, OBDHINIM(JDX,IP,IV) )
               OBDHI_IM(JDX,IP,IV) = MIN(1.0, OBDHI_IM(JDX,IP,IV) )
           END DO    !JDX
        END DO       !IP
      END DO         !IV
C
      RETURN
      END
