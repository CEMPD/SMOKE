
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
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDC
        EXTERNAL        CRLF, FINDC

C...........   Sorting index
        INTEGER, ALLOCATABLE :: INDX( : )

C...........   Concatonated SIC and SCC for all sources
        CHARACTER(LEN=SICLEN3+SCCLEN3), ALLOCATABLE :: CSICSCC( : )

C...........   Other local variables
        INTEGER          J, J1, J2, L1, L2, S
        INTEGER          IOS                 ! allocate i/o status
        INTEGER          FIP                 ! current cntry/st/co code
        INTEGER          PFIP                ! previous iteration cntry/st/co 
        INTEGER          PSIC                ! previous iteration SIC 
        INTEGER          SIC                 ! tmp SIC 

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine

        CHARACTER*8            FMTSIC        ! format buffer for SIC
        CHARACTER*300          MESG          ! message buffer
        CHARACTER(LEN=SCCLEN3) PSCC          ! previous iteration SCC
        CHARACTER(LEN=SCCLEN3) PSCCL         ! previous iteration left SCC
        CHARACTER(LEN=SCCLEN3) SCCL          ! tmp left SCC
        CHARACTER(LEN=SCCLEN3) TSCC          ! tmp SCC
        CHARACTER(LEN=SICLEN3) CSIC          ! tmp char SIC
        CHARACTER(LEN=SICLEN3) PCSIC         ! previous char SIC

        CHARACTER*16  :: PROGNAME = 'GENUSLST' ! program name

C***********************************************************************
C   begin body of subroutine GENUSLST

C.........  Include all processing in a firstime loop, because this information
C           only needs to be created once per program run, and all of the
C           outputs are stored in MODLISTS

        IF( FIRSTIME ) THEN

            WRITE( FMTSIC, 94300 ) '(I', SICLEN3, '.', SICLEN3, ')'

C.............  Allocate memory for sorting index
            ALLOCATE( INDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )

C.............  Check if IFIP is allocated.  
C.............  If it is, generate unique list of country/state/county codes
            IF( ALLOCATED( IFIP ) ) THEN

C.................  Count number of unique codes
                PFIP = IMISS3
                J1 = 0
                DO S = 1, NSRC
 
                    FIP = IFIP( S )
                    IF( FIP .NE. PFIP ) J1 = J1 + 1
                    PFIP = FIP

                END DO
                NINVIFIP = J1

C.................  Allocate memory for country/state/county lists
                ALLOCATE( INVIFIP( NINVIFIP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVIFIP', PROGNAME )

C.................  Create unique country/state/county codes list
                PFIP = IMISS3
                J1 = 0
                DO S = 1, NSRC
 
                    FIP = IFIP( S )
                    IF( FIP .NE. PFIP ) THEN
                        J1 = J1 + 1
                        INVIFIP( J1 ) = FIP
                        PFIP = FIP
                    END IF

                END DO

            END IF   ! End of IFIP allocated or not

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
                J1 = 0        ! to skipp TSCC = 0
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

            END IF   ! End SCC processing
            NINVSCL =  J2

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

C.................  Create unique SIC list
                PSIC = -9
                J1 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    SIC = ISIC( J )

                    IF( SIC .NE. PSIC ) THEN
                        J1 = J1 + 1
                        INVSIC( J1 ) = SIC
                        PSIC = SIC
                    END IF

                END DO 

C.................  Create list of SCCs for each SIC, based on the contents of
C                   the inventory...

C.................  Allocate memory for combo string
                ALLOCATE( CSICSCC( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSICSCC', PROGNAME )

C.................  Initialize SIC sorting index and create combo SIC/SCC string
                DO S = 1, NSRC

                    SIC  = ISIC( S )
                    TSCC = CSCC( S )

                    WRITE( CSIC, FMTSIC ) SIC

                    INDX( S ) = S
                    CSICSCC( S ) = CSIC // TSCC
                   
                END DO

C.................  Sort SIC/SCC combination
                CALL SORTIC( NSRC, INDX, CSICSCC )

C.................  Count the total number of SCCs for all SICs
                J1 = 0
                PCSIC = '-9'
                PSCC  = '-9'
                L1    = SICLEN3 + 1
                L2    = SICLEN3 + SCCLEN3
                DO S = 1, NSRC

                    J = INDX( S )
 
                    CSIC = CSICSCC( J )(  1:SICLEN3 )
                    TSCC = CSICSCC( J )( L1:L2      )

                    IF( CSIC .NE. PCSIC .OR. 
     &                  TSCC .NE. PSCC       ) THEN
                        J1 = J1 + 1
                        PCSIC = CSIC
                        PSCC  = TSCC
                    END IF

                END DO
                NSCCPSIC = J1

C.................  Allocate memory for SCCs per SICs array and indices
                ALLOCATE( SCCPSIC( NSCCPSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SCCPSIC', PROGNAME )
                ALLOCATE( IBEGSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IBEGSIC', PROGNAME )
                ALLOCATE( IENDSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IENDSIC', PROGNAME )

C.................  Store the total SCCs for all SICs and the needed indexes
                J1 = 0
                J2 = 0
                PCSIC = '-9'
                PSCC  = '-9'
                DO S = 1, NSRC

                    J = INDX( S )
 
                    CSIC = CSICSCC( J )(  1:SICLEN3 )
                    TSCC = CSICSCC( J )( L1:L2      )

                    IF( CSIC .NE. PCSIC .OR. 
     &                  TSCC .NE. PSCC       ) THEN
                        J1 = J1 + 1
                        IF( J1 .LE. NSCCPSIC ) SCCPSIC( J1 ) = TSCC
                        PSCC  = TSCC
                    END IF

                    IF( CSIC .NE. PCSIC ) THEN
                        J2 = J2 + 1                        
                        IF( J2 .LE. NINVSIC ) IBEGSIC( J2 ) = J1
                        PCSIC = CSIC
                    END IF

                    IF( J2 .LE. NINVSIC ) IENDSIC( J2 ) = J1

                END DO

C.................  Confirm that SCC-per-SIC count was correct
                IF( J1 .NE. NSCCPSIC ) THEN
                    MESG = 'INTERNAL ERROR: SCC-per-SIC count was '//
     &                     'insufficient for processing.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C...,,,,..........  Confirm that SIC count is same as J2 in above loop
                IF( J2 .NE. NINVSIC ) THEN

                    MESG = 'INTERNAL ERROR: SIC count in SCC-per-SIC '//
     &                     'loop is not equal to NINVSIC.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

            END IF   ! End SIC processing

            FIRSTIME = .FALSE.

C.............  Deallocate local allocatable arrays
            IF( ALLOCATED( CSICSCC ) ) DEALLOCATE( CSICSCC )

            DEALLOCATE( INDX )

        END IF  ! End of firstime

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE GENUSLST
