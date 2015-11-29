
        SUBROUTINE APPLREAC( NSRC, NREAC, KEY1, KEY2, PRJFLAG, LMKTPON, 
     &                       EMIS, EMIST, IDX, REPEM, PRJFC, MKTPEN,
     &                       RMAT, RCTINFO )

C***********************************************************************
C  subroutine body starts at line 59
C
C  DESCRIPTION:
C      This subroutine applies the reactivity matrix for the species of
C      interest (indicated by the KEY, which has already been checked to
C      ensure species is in the reactivity matrix).  The replacement
C      emissions are applied to the sources in the array of inventory 
C      emissions.  The base emissions from the reactivity matrix are used
C      to generate a temporal factor from the ratio of base emissions to
C      the inventory emissions, which could be temporalized.  If the inventory
C      emissions are already projected, the projection factor from the 
C      reactivity matrix is applied.  
C
C      It also populates a market penetration array for all sources.  The
C      default value is 0 for sources that don't get reactivity controls.  If
C      the flag is set to ignore market penetration, then the penetration for
C      reactivity sources is 1.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC             ! number inventory sources
        INTEGER     , INTENT (IN) :: NREAC            ! no. in react matrix
        INTEGER     , INTENT (IN) :: KEY1             ! position in EMIS
        INTEGER     , INTENT (IN) :: KEY2             ! position in RMAT
        LOGICAL     , INTENT (IN) :: PRJFLAG          ! true: emissions are proj 
        LOGICAL     , INTENT (IN) :: LMKTPON          ! true: use mkt pen 
        REAL        , INTENT (IN) :: EMIS   ( NSRC,* )! inventory emissions
        REAL        , INTENT (IN) :: EMIST  ( NSRC,* )! inv or hourly emissions
        INTEGER     , INTENT (IN) :: IDX    ( NREAC ) ! src idx of react. matrix
        REAL        , INTENT (IN) :: REPEM  ( NREAC ) ! replacement emissions
        REAL        , INTENT (IN) :: PRJFC  ( NREAC ) ! stored orig. speciation
        REAL        , INTENT (IN) :: MKTPEN ( NREAC ) ! market penetration
        REAL        , INTENT (IN) :: RMAT   ( NREAC,* ) ! reactivity facs
        REAL        , INTENT(OUT) :: RCTINFO( NSRC,2 )! per-source react. info

C.........  Other local variables
        INTEGER         J, S        !  counters and indices

        REAL            EM          !  tmp emission value
        REAL            PFAC        !  tmp projection factor

        CHARACTER(16) :: PROGNAME = 'APPLREAC' ! program name

C***********************************************************************
C   begin body of subroutine APPLREAC

C.........  Initialize reactivity arrays for all sources
        RCTINFO( :,1 ) = EMIS( :,KEY1 )
        RCTINFO( :,2 ) = 0.

C.........  Exit from subroutine if reactivity matrix does not apply for this
C           species
        IF( KEY2 .LE. 0 ) RETURN

C.........  Note that inventory and hourly (if present) have already been 
C           checked to insure the same base and projection years.  
C.........  Compute the reactivity-controlled emissions by source and set the
C           market penetration.
        PFAC = 1.
        DO J = 1, NREAC

C.............  Retrieve the source index from the reactivity matrix
            S = IDX( J )

C.............  If the inventory emissions are already projected, apply the
C               reactivity projection factor
            IF( PRJFLAG ) PFAC = PRJFC( J )

C.............  Compute output emissions as the product of the replacement 
C               emissions, times the projection factor, times the temporal 
C               adjustment (or 1. if EMIS=EMIST). Screen for divide by zero.
            EM = 0.
            IF( EMIS( S,KEY1 ) .GT. 0 ) THEN

                EM = REPEM( J ) * PFAC * 
     &               EMIST( S,KEY1 ) / EMIS( S,KEY1 )

            END IF

            RCTINFO( S,1 ) = EM * RMAT( J,KEY2 )

C.............  If using market penetration, then set to value in matrix
            IF( LMKTPON ) THEN

                RCTINFO( S,2 ) = MKTPEN( J )

C.............  If ignoring market penetration, then set to unity
            ELSE

                RCTINFO( S,2 ) = 1.

            END IF

        END DO

        RETURN

        END SUBROUTINE APPLREAC
