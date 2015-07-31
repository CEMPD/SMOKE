
      REAL FUNCTION PLUMRIS( HS, TS, VS, DS )

C***********************************************************************
C  subroutine body starts at line  82
C
C  DESCRIPTION:  Computes effective plume height using Briggs algorithm.  See:
C   1) Briggs, Gary A., 1971: Some Recent Analyses of Plume Rise Observation
C      pp 1029 - 1032 in PROCEEDINGS OF THE SECOND INTERNATIONAL CLEAN AIR
C      CONGRESS, edited by H. M. Englun and W. T. Beery. Academic Press, 
C      New York.
C   2) Briggs, Gary A., 1972: Discussion on Chimney Plumes in Neutral
C      and Stable Surroundings. ATMOS. ENVIRON. 6, 507 - 510. (Jul 72).
C
C  REVISION  HISTORY:
C       Copied 8/99 from plumris.F v4.2 in SMOKE prototype
C       Adapted 10/95 by Carlie J. Coats, Jr., from UAM/EPS BEH072().
C       Discontinuity in EPS DISTF calculation at F=55 resolved so
C       that DISTF --> 595.0  as  F --> 55.0 from either side.
C       Uses standard conditions temperature T=293 deg K, 
C                                pressure    P=960 mb, 
C                                wind speed  U=  2 m/s,
C                       Pasquill stability KST=  2
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
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
C************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CONST3.EXT'      ! physical and mathematical constants


C...........   ARGUMENTS and their descriptions:

        REAL      HS    !  physical stack height (m)
        REAL      TS    !  stack gas temperature (deg k)
        REAL      VS    !  stack gas exit velocity (m/sec)
        REAL      DS    !  inside stack diameter (m)

C...........   PARAMETERS and their descriptions:

        REAL        T           !  default ambient air temperature (deg k)
        REAL        P           !  default ambient air pressure (mb)
        REAL        U           !  default wind speed (m/sec)
        REAL        D3          !  one-third

        PARAMETER ( T    = 293.0 ,
     &              P    = 960.0 ,
     &              U    =   2.0 ,
     &              D3   =   1.0 / 3.0 )

C...........   LOCAL VARIABLES:

        REAL      F     ! buoyancy flux (m**4/sec**3)
C       REAL      DELHF ! final plume rise (m)
C       REAL      DISTF ! distance of final plume rise from source (m)


C***********************************************************************
C   begin body of function PLUMRIS

        IF ( TS .LE. T ) THEN
            PLUMRIS = MAX( HS, 3.0 )
            RETURN
        END IF

        IF ( DS .LE. 0.0 ) DS = 0.2
        IF ( HS .LE. 0.0 ) HS = 3.0
        IF ( VS .LE. 0.0 ) VS = 0.5
        
        F = 0.25 * GRAV * VS * DS * DS * ( TS - T ) / TS

C        IF( F .LT. 55.0 ) THEN
C            DISTF = 3.5 * 13.8906395 * F**0.625
C        ELSE
C            DISTF = 3.5 * 34.22187854 * F**0.4
C        END IF
C        DELHF = ( 1.6 / U ) * ( F * DISTF * DISTF )**D3
C        PLUMRIS = HS + DELHF

        IF( F .LT. 55.0 ) THEN
            PLUMRIS = HS + 21.31311057 * F**0.75 / U
        ELSE
            PLUMRIS = HS + 38.87776061 * F**0.6 / U
        END IF

        RETURN
        END

