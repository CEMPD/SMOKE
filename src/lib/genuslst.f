
        SUBROUTINE GENUSLST

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine processes generates the source-category-specific lists
C      of inventory characteristics in the MODLISTS module. 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDC
        EXTERNAL   CRLF, FINDC

C...........   Sorting index
        INTEGER, ALLOCATABLE :: INDX( : )

C...........   Other local variables
        INTEGER          J, J1, J2, S
        INTEGER          IOS                 ! allocate i/o status
        INTEGER          PSIC                ! previous iteration SIC 
        INTEGER          SIC                 ! tmp SIC 

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine

        CHARACTER*300          MESG          ! message buffer
        CHARACTER(LEN=SCCLEN3) PSCC          ! previous iteration SCC
        CHARACTER(LEN=SCCLEN3) PSCCL         ! previous iteration left SCC
        CHARACTER(LEN=SCCLEN3) SCCL          ! tmp left SCC
        CHARACTER(LEN=SCCLEN3) TSCC          ! tmp SCC

        CHARACTER*16  :: PROGNAME = 'GENUSLST' ! program name

C***********************************************************************
C   begin body of subroutine GENUSLST

C.........  Include all processing in a firstime loop, because this information
C           only needs to be created once per program run, and all of the
C           outputs are stored in MODLISTS

        IF( FIRSTIME ) THEN

C.............  Allocate memory for sorting index
            ALLOCATE( INDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )
            
C.............  Check if CSCC is allocated.  
C.............  If it is, generate unique list of SCCs and left SCCs
            IF( ALLOCATED( CSCC ) ) THEN

C.................  Initialize SCC sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO

C.................  Sort all SCCs in the point sources inventory in increasing
C                   order
                CALL SORTIC( NSRC, INDX, CSCC )

C.................  Count number of unique SCCs
                PSCC = '-9'
                J1 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    TSCC = CSCC( J )

                    IF( TSCC .NE. PSCC ) J1 = J1 + 1

                    PSCC = TSCC

                END DO 
                NINVSCC = J1

C.................  Allocate memory for SCC lists
                ALLOCATE( INVSCC( NINVSCC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSCC', PROGNAME )
                ALLOCATE( INVSCL( NINVSCC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSCL', PROGNAME )

C.................  Create unique SCCs lists
                PSCC = '-9'
                PSCCL = '-9'
                J1 = 0
                J2 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    TSCC  = CSCC( J )
                    SCCL  = CSCC( J )( 1:LSCCEND )

                    IF( TSCC .NE. PSCC ) THEN
                        J1 = J1 + 1
                        INVSCC( J1 ) = TSCC
                        PSCC = TSCC
                    END IF

                    IF( SCCL .NE. PSCCL ) THEN
                        J2 = J2 + 1
                        INVSCL( J2 ) = SCCL
                        PSCCL = SCCL
                    END IF

                END DO 
                NINVSCL = J2

            END IF   ! End SCC processing

C.............  Check if ISIC is allocated.  
C.............  If it is, generate unique list of SICs, and generate list of
C               where which SCCs go with which SICs
            IF( ALLOCATED( ISIC ) ) THEN

C.................  Initialize SIC sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO

C.................  Sort all SICs in the point sources inventory in increasing
C                   order
                CALL SORTI1( NSRC, INDX, ISIC )

C.................  Count number of unique SICs
                PSIC = -9
                J1 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    SIC = ISIC( J )

                    IF( SIC .NE. PSIC ) J1 = J1 + 1

                    PSIC = SIC

                END DO 
                NINVSIC = J1

C.................  Allocate memory for SIC lists
                ALLOCATE( INVSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSIC', PROGNAME )
                ALLOCATE( IBEGSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IBEGSIC', PROGNAME )
                ALLOCATE( IENDSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IENDSIC', PROGNAME )

C.................  Initialize the position of the earliest and latest SCC
C                   in the unique SCC list for each SIC.
                IBEGSIC = NINVSCC   ! array
                IENDSIC = 0         ! array

C.................  Populate unique list of SICs and set the position of the
C                   earliest and latest SCC in the unique SCCs list.  This
C                   assumes that the SCCs are grouped together by SIC.

C NOTE: Need to check to see if this assumption is true.  Furthermore, need to
C       check that there is only one SIC per SCC.

                PSIC = -9
                J1 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    SIC = ISIC( J )

C.....................  Find this SCC in SCC list
                    J2 = FINDC( CSCC( J ), NINVSCC, INVSCC )

                    IF( SIC .NE. PSIC ) THEN

                        J1 = J1 + 1
                        INVSIC( J1 ) = SIC
                        PSIC = SIC

                        IF( J2 .LT. IBEGSIC( J1 ) ) IBEGSIC( J1 ) = J2
                        IF( J2 .GT. IENDSIC( J1 ) ) IENDSIC( J1 ) = J2

                    END IF

                END DO 

            END IF   ! End SIC processing

            FIRSTIME = .FALSE.

C.............  Deallocate local allocatable arrays
            DEALLOCATE( INDX )

        END IF  ! End of firstime

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GENUSLST
