
        SUBROUTINE GENEFMET( NSRC, MXPSI, NVLDTMM, NIACT, 
     &                       METIDX, TIPSI )
   
C***********************************************************************
C  subroutine GENEFMET body starts at line < >
C
C  DESCRIPTION:
C
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module contains the cross-reference tables
        USE MODXREF

C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS
        INTEGER   FIND1
        EXTERNAL  FIND1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NSRC            ! no. sources
        INTEGER, INTENT (IN) :: MXPSI           ! max. no. PSIs for all actvtys
        INTEGER, INTENT (IN) :: NVLDTMM         ! no. valid tempr combos
        INTEGER, INTENT (IN) :: NIACT           ! no. activites
        INTEGER, INTENT (IN) :: METIDX( NSRC,4 )! min/max combo indicator
        INTEGER, INTENT(OUT) :: TIPSI( MXPSI,NVLDTMM,NIACT ) ! PSI/min-max flag

C...........  Local variables
        INTEGER      H, I, J, K, S, V   ! indices and counters

        INTEGER      LPSI    ! parameter scheme index from previous iteration
        INTEGER      PSI     ! tmp parameter scheme index

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

C***********************************************************************
C   begin body of subroutine GENEFMET

C.........  Initialize PSI/tmpr-combo flag to zero
        IF( FIRSTIME ) THEN
            TIPSI = 0   ! array
            FIRSTIME = .FALSE.
        END IF

c.........  Loop through sources
        DO S = 1, NSRC

C.............  Loop through activities
            DO V = 1, NIACT

C.................  Get index to unsorted PSIs cross-reference
                K = EFSIDX( S,V )

C.................  Loop through PSIs and flag PSI/min-max-tmpr
                LPSI = -9
                DO H = 1, 24
 
                    PSI = IPSIA( K,H )

                    IF( PSI .NE. LPSI ) THEN

C.........................  Find PSI in sorted list for position
                        I = FIND1( PSI, NPSI( V ), PSILIST( 1,V ) )

C.........................  For each min/max combo, set flag
                        CALL SET_TIPSI( METIDX( S,1 ) )
                        CALL SET_TIPSI( METIDX( S,2 ) )
                        CALL SET_TIPSI( METIDX( S,3 ) )
                        CALL SET_TIPSI( METIDX( S,4 ) )

                        LPSI = PSI

                    END IF

                END DO ! end loop on hours

            END DO     ! end loop on activities

        END DO         ! end loop on sources

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram sets TIPSI if the index J is
C               non-zero.  All other variables are accessed from main program

            SUBROUTINE SET_TIPSI( J )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: J    ! 2nd index of TIPSI

C----------------------------------------------------------------------

            IF( J .GT. 0 ) THEN
                TIPSI( I, J, V ) = 1
            END IF

            RETURN

            END SUBROUTINE SET_TIPSI

        END SUBROUTINE GENEFMET



