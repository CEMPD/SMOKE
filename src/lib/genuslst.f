
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

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: CSOURC, IFIP, CSCC, ISIC, CMACT, 
     &                      CORIS, CBLRID, CPDESC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, NINVSCC, NINVSCL, NINVSIC, 
     &                      NINVSIC2, NINVMACT, NINVORIS, 
     &                      INVIFIP, INVSCC, INVSCL, INVSIC,
     &                      INVSIC2, INVMACT, INVORIS,
     &                      INVORFP, IORSMTCH, INVODSC, ORISBLR,
     &                      OBSRCBG, OBSRCNT, ORISPNT, OPSRCBG,
     &                      OPSRCNT, NORISBLR, NORISPNT, ORISFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, LSCCEND

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         FINDC
        LOGICAL         SETSCCTYPE
        
        EXTERNAL        CRLF, FINDC, SETSCCTYPE

C...........   Sorting index
        INTEGER, ALLOCATABLE :: INDX( : )

C...........   Local allocateable arrays for ORIS lists
        INTEGER, ALLOCATABLE :: FOIDXA  ( : )  ! sorting index for oris
        INTEGER, ALLOCATABLE :: OBIDXA  ( : )  ! sorting index for oris//blr
        INTEGER, ALLOCATABLE :: OPIDXA  ( : )  ! sorting index for oris//point
        INTEGER, ALLOCATABLE :: INVORFPA( : )  ! FIPS code for ORIS IDs
        INTEGER, ALLOCATABLE :: OBSRCBGA( : )  ! unsrtd 1st src per ORIS/boiler
        INTEGER, ALLOCATABLE :: OBSRCNTA( : )  ! unsrtd src count per ORIS/boiler
        INTEGER, ALLOCATABLE :: OPSRCBGA( : )  ! unsrtd 1st src per ORIS/point
        INTEGER, ALLOCATABLE :: OPSRCNTA( : )  ! unsrtd src count per ORIS/point

        CHARACTER(ORSLEN3), ALLOCATABLE :: INVORISA( : )  ! ORIS
        CHARACTER(OBRLEN3), ALLOCATABLE :: ORISBLRA( : )  ! ORIS // boiler
        CHARACTER(OPTLEN3), ALLOCATABLE :: ORISPNTA( : )  ! ORIS // point
        CHARACTER(DSCLEN3), ALLOCATABLE :: INVODSCA( : ) ! plant description from inventory

C...........   Other local variables
        INTEGER          I, J, J1, J2, L1, L2, N, S
        INTEGER          IOS                 ! allocate i/o status
        INTEGER          FIP                 ! current cntry/st/co code
        INTEGER          PFIP                ! previous iteration cntry/st/co 
        INTEGER          PSIC                ! previous iteration SIC 
        INTEGER          PSIC2               ! previous iteration 2-digit SIC 
        INTEGER          SIC                 ! tmp SIC 
        INTEGER          SIC2                ! tmp 2-digit SIC 

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
        LOGICAL, SAVE :: FIRSTORS = .TRUE.   ! true: first run of ORIS arrays
        LOGICAL          SCCFLAG             ! true: SCC type is different from previous

        CHARACTER(8)       FMTSIC        ! format buffer for SIC
        CHARACTER(300)     MESG          ! message buffer
        CHARACTER(SCCLEN3) PSCC          ! previous iteration SCC
        CHARACTER(SCCLEN3) PSCCL         ! previous iteration left SCC
        CHARACTER(SCCLEN3) SCCL          ! tmp left SCC
        CHARACTER(SCCLEN3) TSCC          ! tmp SCC
        CHARACTER(SICLEN3) CSIC          ! tmp char SIC
        CHARACTER(SICLEN3) PCSIC         ! previous char SIC
        CHARACTER(MACLEN3) TMACT         ! tmp char MACT code
        CHARACTER(MACLEN3) PMACT         ! previous char MACT code
        CHARACTER(BLRLEN3) BLID          ! tmp boiler ID
        CHARACTER(BLRLEN3) PBLID         ! previous boiler ID
        CHARACTER(ORSLEN3) CORS          ! tmp DOE plant ID
        CHARACTER(ORSLEN3) PCORS         ! previous DOE plant ID
        CHARACTER(OBRLEN3) PCORSBLR      ! previous DOE plant ID // boiler
        CHARACTER(OPTLEN3) PCORSPNT      ! previous DOE plant ID // point
        CHARACTER(CHRLEN3) PNT           ! point (IDA char1)
        CHARACTER(CHRLEN3) PPNT          ! previous point (IDA char1)
        CHARACTER(DSCLEN3) PDSC          ! tmp plant description

        CHARACTER(16) :: PROGNAME = 'GENUSLST' ! program name

C***********************************************************************
C   begin body of subroutine GENUSLST

C.........  Include all processing in a firstime loop, because this information
C           only needs to be created once per program run, and all of the
C           outputs are stored in MODLISTS

        IF( FIRSTIME ) THEN

            MESG = 'Generating unique lists from inventory data...'
            CALL M3MSG2( MESG )

            WRITE( FMTSIC, 94300 ) '(I', SICLEN3, '.', SICLEN3, ')'

C.............  Allocate memory for sorting index
            ALLOCATE( INDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )

C.............  Check if IFIP is allocated.  
C.............  If it is, generate unique list of country/state/county codes
            IF( ASSOCIATED( IFIP ) ) THEN

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
            IF( ASSOCIATED( CSCC ) ) THEN

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
                    
C.....................  Set type of SCC                
                    SCCFLAG = SETSCCTYPE( TSCC )
                    SCCL  = TSCC( 1:LSCCEND )

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
            IF( ASSOCIATED( ISIC ) ) THEN

C.................  Initialize SIC sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO

C.................  Sort all SICs in the inventory in increasing order
                CALL SORTI1( NSRC, INDX, ISIC )

C.................  Count number of unique SICs and number of unique
C                   2-digit SICs
                PSIC  = -9
                PSIC2 = -9
                J1 = 0
                J2 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    SIC  = ISIC( J )
                    SIC2 = SIC/100

                    IF( SIC  .NE. PSIC  ) J1 = J1 + 1
                    IF( SIC2 .NE. PSIC2 ) J2 = J2 + 1 

                    PSIC  = SIC
                    PSIC2 = SIC2

                END DO 
                NINVSIC  = J1
                NINVSIC2 = J2

C.................  Allocate memory for SIC lists
                ALLOCATE( INVSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSIC', PROGNAME )
                ALLOCATE( INVSIC2( NINVSIC2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSIC2', PROGNAME )

C.................  Create unique SIC list
                PSIC  = -9
                PSIC2 = -9
                J1 = 0
                J2 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    SIC = ISIC( J )
                    SIC2 = SIC/100

                    IF( SIC .NE. PSIC ) THEN
                        J1 = J1 + 1
                        INVSIC( J1 ) = SIC
                        PSIC = SIC
                    END IF

                     IF( SIC2 .NE. PSIC2 ) THEN
                        J2 = J2 + 1
                        INVSIC2( J2 ) = SIC2
                        PSIC2 = SIC2
                    END IF

               END DO 

            END IF   ! End SIC processing

C.............  Check if CMACT is allocated.  
C.............  If it is, generate unique list of MACT codes
            IF( ASSOCIATED( CMACT ) ) THEN

C.................  Initialize MACT sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO                

C.................  Sort all MACTs in the inventory in increasing order
                CALL SORTIC( NSRC, INDX, CMACT )

C.................  Count number of unique MACTs
                PMACT = '-9'
                J1 = 0
                DO S = 1, NSRC
                
                    J = INDX( S )
                    
                    TMACT = CMACT( J )
                    
                    IF( TMACT /= PMACT ) J1 = J1 + 1
                    
                    PMACT = TMACT
                    
                END DO
                NINVMACT = J1

C.................  Allocate memory for MACT lists
                ALLOCATE( INVMACT( NINVMACT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVMACT', PROGNAME )

C.................  Create unique MACTs list
                PMACT = '-9'
                J1 = 0
                DO S = 1, NSRC
                
                    J = INDX( S )
                    
                    TMACT = CMACT( J )
                    
                    IF( TMACT /= PMACT ) THEN
                        J1 = J1 + 1
                        INVMACT( J1 ) = TMACT
                        PMACT = TMACT
                    END IF
                    
                END DO
                
            END IF    ! End MACT processing
                
            DEALLOCATE( INDX )

            FIRSTIME = .FALSE.

        END IF  ! End of firstime

C.........  Create list of FIPs, ORIS IDs and boiler IDs from the inventory
C           and how many sources for each...
        IF ( ORISFLAG             .AND.
     &       FIRSTORS             .AND.
     &       ALLOCATED ( CORIS  ) .AND.
     &       ALLOCATED ( CBLRID )       ) THEN

            MESG = 'Generating ORIS lists...'
            CALL M3MSG2( MESG )

C.............  First, count the number of unique records
            NINVORIS = 0
            NORISBLR = 0
            NORISPNT = 0
            PFIP     = -9
            PCORS    = ' '
            PBLID    = ' '
            PPNT     = ' '
            DO S = 1, NSRC

                FIP  = IFIP  ( S )
                CORS = CORIS ( S )
                BLID = CBLRID( S )
                PNT  = CSOURC( S )( PTBEGL3(3):PTENDL3(3) )
 
C.................  Skip missing ORIS IDs
                IF ( CORS .EQ. ' ' ) CYCLE

C.................  Count unique ORIS IDs
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS      ) THEN
                    NINVORIS = NINVORIS + 1
                END IF

C.................  Count unique ORIS // boiler combos
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS .OR.  
     &               BLID .NE. PBLID      ) THEN
                    NORISBLR = NORISBLR + 1

                END IF

C.................  Count unique ORIS // point combos
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS .OR.  
     &               PNT  .NE. PPNT       ) THEN
                    NORISPNT = NORISPNT + 1

                END IF

                PFIP  = FIP
                PCORS = CORS
                PBLID = BLID
                PPNT  = PNT

            END DO

C.............  Allocate memory for sorted and unsorted FIPS/ORIS list
            ALLOCATE( FOIDXA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FOIDXA', PROGNAME )
            ALLOCATE( INVORISA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORISA', PROGNAME )
            ALLOCATE( INVORIS( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORIS', PROGNAME )
            ALLOCATE( INVORFPA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORFPA', PROGNAME )
            ALLOCATE( INVORFP( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORFP', PROGNAME )
            ALLOCATE( IORSMTCH( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IORSMTCH', PROGNAME )
            ALLOCATE( INVODSCA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVODSCA', PROGNAME )
            ALLOCATE( INVODSC( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVODSC', PROGNAME )
            INVORIS = ' '
            INVORFP = 0
            IORSMTCH = .FALSE.
            INVODSCA = ' '

C.............  Allocate memory for sorted and unsorted ORIS/boiler
            ALLOCATE( OBIDXA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBIDXA', PROGNAME )
            ALLOCATE( ORISBLRA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ORISBLRA', PROGNAME )
            ALLOCATE( ORISBLR( NORISBLR ), STAT=IOS )   ! ORIS/boiler combo
            CALL CHECKMEM( IOS, 'ORISBLR', PROGNAME )
            ALLOCATE( OBSRCBGA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCBGA', PROGNAME )
            ALLOCATE( OBSRCBG( NORISBLR ), STAT=IOS )   ! starting src no.
            CALL CHECKMEM( IOS, 'OBSRCBG', PROGNAME )
            ALLOCATE( OBSRCNTA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCNTA', PROGNAME )
            ALLOCATE( OBSRCNT( NORISBLR ), STAT=IOS )   ! count of sources
            CALL CHECKMEM( IOS, 'OBSRCNT', PROGNAME )
            ORISBLR = ' '
            OBSRCBG = 0
            OBSRCNT = 0

C.............  Allocate memory for sorted and unsorted ORIS/point
            ALLOCATE( OPIDXA( NORISPNT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OPIDXA', PROGNAME )
            ALLOCATE( ORISPNTA( NORISPNT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ORISPNTA', PROGNAME )
            ALLOCATE( ORISPNT( NORISPNT ), STAT=IOS )   ! ORIS/point combo
            CALL CHECKMEM( IOS, 'ORISPNT', PROGNAME )
            ALLOCATE( OPSRCBGA( NORISPNT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OPSRCBGA', PROGNAME )
            ALLOCATE( OPSRCBG( NORISPNT ), STAT=IOS )   ! starting src no.
            CALL CHECKMEM( IOS, 'OPSRCBG', PROGNAME )
            ALLOCATE( OPSRCNTA( NORISPNT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OPSRCNTA', PROGNAME )
            ALLOCATE( OPSRCNT( NORISPNT ), STAT=IOS )   ! count of sources
            CALL CHECKMEM( IOS, 'OPSRCNT', PROGNAME )
            ORISPNT = ' '
            OPSRCBG = 0
            OPSRCNT = 0

C.............  Store arrays
            NINVORIS    = 0
            NORISBLR = 0
            NORISPNT = 0
            PFIP     = -9
            PCORS    = ' '
            PBLID    = ' '
            PPNT     = ' '
            DO S = 1, NSRC

                FIP  = IFIP  ( S )
                CORS = CORIS ( S )
                BLID = CBLRID( S )
                PNT  = CSOURC( S )( PTBEGL3(3):PTENDL3(3) )
                PDSC = CPDESC( S )

C.................  Skip missing ORIS IDs
                IF ( CORS .EQ. ' ' ) CYCLE

C.................  Unsorted oris/FIPS array
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS      ) THEN
                    NINVORIS = NINVORIS + 1
                    INVORISA( NINVORIS ) = CORS
                    INVORFPA( NINVORIS ) = FIP
                    INVODSCA( NINVORIS ) = PDSC
                    FOIDXA  ( NINVORIS ) = NINVORIS
                END IF

C.................  Unsorted oris/boiler array
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS .OR.  
     &               BLID .NE. PBLID      ) THEN
                    
                    NORISBLR = NORISBLR + 1
                    OBIDXA  ( NORISBLR ) = NORISBLR
                    ORISBLRA( NORISBLR ) = CORS // BLID
                    OBSRCBGA( NORISBLR ) = S
                    OBSRCNTA( NORISBLR ) = 1

                ELSE IF ( FIP  .EQ. PFIP  .AND.
     &                    CORS .EQ. PCORS .AND.
     &                    BLID .EQ. PBLID       ) THEN

                    OBSRCNTA( NORISBLR ) = OBSRCNTA( NORISBLR ) + 1

                END IF

C.................  Unsorted oris/point array
                IF ( FIP  .NE. PFIP  .OR.
     &               CORS .NE. PCORS .OR.  
     &               PNT  .NE. PPNT       ) THEN
                    
                    NORISPNT = NORISPNT + 1
                    OPIDXA  ( NORISPNT ) = NORISPNT
                    ORISPNTA( NORISPNT ) = CORS // PNT
                    OPSRCBGA( NORISPNT ) = S
                    OPSRCNTA( NORISPNT ) = 1

                ELSE IF ( FIP  .EQ. PFIP  .AND.
     &                    CORS .EQ. PCORS .AND.
     &                    PNT  .EQ. PPNT       ) THEN

                    OPSRCNTA( NORISPNT ) = OPSRCNTA( NORISPNT ) + 1

                END IF

                PFIP  = FIP
                PCORS = CORS
                PBLID = BLID
                PPNT  = PNT

            END DO

C.............  Sort arrays
            CALL SORTIC( NINVORIS, FOIDXA, INVORISA )
            CALL SORTIC( NORISBLR, OBIDXA, ORISBLRA )
            CALL SORTIC( NORISPNT, OPIDXA, ORISPNTA )

C.............  Store sorted arrays
c note: add check to ensure that their is a 1-1 ORIS/FIP assignment (the same
c N: ORIS doesn't appear in two FIPS)
            N = 0
            PCORS = ' '
            DO I = 1, NINVORIS
                J = FOIDXA( I )

                IF ( INVORISA( J ) .NE. PCORS ) THEN
                    N = N + 1
                    INVORIS( N ) = INVORISA( J )
                    INVORFP( N ) = INVORFPA( J )
                    INVODSC( N ) = INVODSCA( J )
                END IF

                PCORS = INVORISA( J )

            END DO
            NINVORIS = N

c note: with this code, if there is a oris/boiler reduction because of
c    n: duplicates that are out of order in the inventory, (multiple plants
c    n: have the same ORISID//boiler that aren't together in the inventory),
c    n: then the count of sources per oris/boiler is wrong.
            N = 0
            PCORSBLR = ' '
            DO I = 1, NORISBLR
                J = OBIDXA( I )

                IF ( ORISBLRA( J ) .NE. PCORSBLR ) THEN
                    N = N + 1
                    ORISBLR( N ) = ORISBLRA( J )
                    OBSRCBG( N ) = OBSRCBGA( J )
                    OBSRCNT( N ) = OBSRCNTA( J )
                END IF

                PCORSBLR = ORISBLRA( J )

            END DO
            NORISBLR = N

            N = 0
            PCORSPNT = ' '
            DO I = 1, NORISPNT
                J = OPIDXA( I )

                IF ( ORISPNTA( J ) .NE. PCORSPNT ) THEN
                    N = N + 1
                    ORISPNT( N ) = ORISPNTA( J )
                    OPSRCBG( N ) = OPSRCBGA( J )
                    OPSRCNT( N ) = OPSRCNTA( J )
                END IF

                PCORSPNT = ORISPNTA( J )

            END DO
            NORISPNT = N

C.............  Deallocate unneeded unsorted arrays
            DEALLOCATE( FOIDXA, INVORISA, INVORFPA, OBIDXA, ORISBLRA, 
     &                  OBSRCBGA, OBSRCNTA )
                
            FIRSTORS = .FALSE.
        
        END IF   ! End boiler processing

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE GENUSLST
