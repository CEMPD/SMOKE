        SUBROUTINE GETPAR( RSOLAR, PRES, ZEN, PARDB, PARDIF )

C***********************************************************************
C  subroutine body starts at line  
C
C  DESCRIPTION:
C  
C        Based on code from Bart Brashers (10/2000), which was based on
C        code from Weiss and Norman (1985).  
C
C
C  PRECONDITIONS REQUIRED:
C     Solar radiation (W/m2) and pressure (mb)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    3/01 : Prototype by JMV
C 
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

      IMPLICIT NONE

C........ Inputs

        REAL , INTENT  (IN) :: RSOLAR   ! modeled or observed total radiation (W/m2)
        REAL , INTENT  (IN) :: PRES     ! atmospheric pressure (mb)
        REAL , INTENT  (IN) :: ZEN      ! solar zenith angle (radians)
 
C........ Outputs

        REAL, INTENT (OUT) :: PARDB     ! direct beam PAR( umol/m2-s)
        REAL, INTENT (OUT) :: PARDIF    ! diffuse PAR ( umol/m2-s)

C...........   PARAMETERS and their descriptions:

        REAL, PARAMETER :: WATT2UMOL = 4.6  ! convert W/m^2 to umol/m^2-s (4.6)

C      
        REAL RATIO              ! transmission fraction for total radiation
        REAL OT                   ! optical thickness
        REAL RDVIS                ! possible direct visible beam (W/m^2)
        REAL RFVIS              ! possible visible diffuse (W/m^2)
        REAL WA                 ! water absorption in near-IR (W/m^2)
        REAL RDIR               ! direct beam in near-IR (W/m^2)
        REAL RFIR               ! diffuse near-IR (W/m^2)
        REAL RVT                ! total possible visible radiation (W/m^2)
        REAL RIRT               ! total possible near-IR radiation (W/m^2)
        REAL FVIS               ! fraction of visible to total 
        REAL FVB                ! fraction of visible that is direct beam
        REAL FVD                ! fraction of visible that is diffuse

        CHARACTER(16) :: PROGNAME = 'GETPAR'   !  program name

        CHARACTER(256)  MESG

C***************************************
C   begin body of subroutine  

C............ Assume that PAR = 0 if zenith angle is greater than 87 degrees
C............ or if solar radiation is zero

        IF (ZEN .GE. 1.51844 .OR. RSOLAR .LE. 0.) THEN
             PARDB  = 0.
             PARDIF = 0.
             RETURN
        ENDIF
  
C............ Compute clear sky (aka potential) radiation terms

        OT    = PRES / 1013.25 / COS(ZEN)              !Atmospheric Optical thickness
        RDVIS = 600. * EXP(-.185 * OT) * COS(ZEN)      !Direct visible beam, eqn (1)
        RFVIS = 0.42 * (600 - RDVIS) * COS(ZEN)        !Visible Diffuse, eqn (3)
        WA    = 1320 * .077 * (2. * OT)**0.3           !water absorption in near-IR, eqn (6)
        RDIR  = (720. * EXP(-0.06 * OT)-WA) * COS(ZEN) !Direct beam near-IR, eqn (4)
        RFIR  = 0.65 * (720. - WA - RDIR) * COS(ZEN)   !Diffuse near-IR, eqn (5)

        RVT   = RDVIS + RFVIS                    !Total visible radiation, eqn (9)
        RIRT  = RDIR + RFIR                      !Total near-IR radiation, eqn (10) 
        FVIS  = RVT/(RIRT + RVT)                 !Fraction of visible to total radiation, eqn 7
        RATIO = RSOLAR /(RIRT + RVT)             !Ratio of "actual" to clear sky solar radiation

C............ Compute fraction of visible that is direct beam

        IF (RATIO .GE. 0.89) THEN
           FVB = RDVIS/RVT * 0.941124
        ELSE IF (RATIO .LE. 0.21) THEN
           FVB = RDVIS/RVT * 9.55E-3
        ELSE
           FVB = RDVIS/RVT * (1.-((0.9 - RATIO)/0.7)**0.666667)
        ENDIF
        FVD = 1. - FVB

C............ Compute PAR (direct beam and diffuse) in umol/m2-sec

        PARDB  = RSOLAR * FVIS * FVB * WATT2UMOL
        PARDIF = RSOLAR * FVIS * FVD * WATT2UMOL      


        RETURN 

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx


C...........   Internal buffering formats............ 94xxx

        END SUBROUTINE GETPAR
