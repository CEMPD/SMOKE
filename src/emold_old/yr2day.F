
        REAL FUNCTION YR2DAY ( YEAR )

C********************************************************************
C       function body starts at line  58
C
C  FUNCTION:
C
C      Returns the year to day conversion factor for a given year
C       
C
C  REVISION HISTORY:
C
C       Create by M Houyoux: 4/97
C
C  ARGUMENT LIST DESCRIPTION:
C
C    Input arguments:
C
C        YEAR   - 4 digit year
C
C    Output arguments:  none
C
C  RETURNS   user response after checking its range; or default.
C
C********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

        IMPLICIT NONE

C.......   ARGUMENTS:

        INTEGER         YEAR  ! 4 digit year

C......................................................................
C       begin YR2DAY

        IF ( MOD( YEAR,4 ) .NE. 0 ) THEN           !  nonleap years

            YR2DAY = 1.0 / 365.0

        ELSE IF ( MOD( YEAR,100 ) .NE. 0 ) THEN    !  noncentury leap years

            YR2DAY = 1.0 / 366.0

        ELSE IF ( MOD( YEAR,400 ) .NE. 0 ) THEN    !  century nonleap years

            YR2DAY = 1.0 / 365.0

        ELSE                                       !  leap centuries

            YR2DAY = 1.0 / 366.0

        END IF

        RETURN

C................   end body of YR2DAY  .......................................

        END

