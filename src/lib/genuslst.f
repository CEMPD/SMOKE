
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
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC, CIFIP, CSCC, CISIC, CINTGR, CMACT, 
     &                      CORIS, CBLRID, CPDESC, CNAICS, CVTYPE

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, NINVSCC, NINVSCL, NINVSIC, 
     &                      NINVSIC2, NINVMACT, NINVORIS, 
     &                      INVCFIP, INVSCC, INVSCL, INVSIC,
     &                      INVSIC2, INVMACT, INVORIS,
     &                      INVORFP, IORSMTCH, INVODSC, ORISBLR,
     &                      OBSRCBG, OBSRCNT, NORISBLR, NOBLRSRC,
     &                      OBSRCNM, ORISFLAG, NINVNAICS, INVNAICS,
     &                      NINVVTYP, INVVTYP

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, LSCCEND

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       INTEGER         FINDC
C       INTEGER         INDEX1
        LOGICAL         SETSCCTYPE, CHKEXPSCC, CHKEXPSIC
        
C        EXTERNAL        CRLF, FINDC, INDEX1, SETSCCTYPE
        EXTERNAL     SETSCCTYPE

C...........   Sorting index
        INTEGER, ALLOCATABLE :: INDX( : )

C...........   Local allocateable arrays for ORIS lists
        INTEGER, ALLOCATABLE :: FOIDXA  ( : )  ! sorting index for oris
        INTEGER, ALLOCATABLE :: OBIDXA  ( : )  ! sorting index for oris//blr
        INTEGER, ALLOCATABLE :: OBSRCNTA( : )  ! unsrtd src count per ORIS/boiler

        CHARACTER(ORSLEN3), ALLOCATABLE :: INVORISA( : )  ! ORIS
        CHARACTER(FIPLEN3), ALLOCATABLE :: INVORFPA( : )  ! FIPS code for ORIS IDs
        CHARACTER(OBRLEN3), ALLOCATABLE :: ORISBLRA( : )  ! ORIS // boiler
        CHARACTER(DSCLEN3), ALLOCATABLE :: INVODSCA( : ) ! plant description from inventory

C...........   Other local variables
        INTEGER          I, J, J1, L1, L2, N, NS, S
        INTEGER       :: J2 = 0
        INTEGER          IOS                 ! allocate i/o status

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
        LOGICAL, SAVE :: FIRSTORS = .TRUE.   ! true: first run of ORIS arrays
        LOGICAL          SCCFLAG             ! true: SCC type is different from previous

        CHARACTER(300)     MESG          ! message buffer
        CHARACTER(FIPLEN3) CFIP          ! current cntry/st/co code
        CHARACTER(FIPLEN3) PCFIP        ! previous iteration cntry/st/co
        CHARACTER(VTPLEN3) PVTYP         ! previous vehicle type
        CHARACTER(VTPLEN3) TVTYP         ! tmp vehicle type
        CHARACTER(SCCLEN3) PSCC          ! previous iteration SCC
        CHARACTER(SCCLEN3) PSCCL         ! previous iteration left SCC
        CHARACTER(SCCLEN3) SCCL          ! tmp left SCC
        CHARACTER(SCCLEN3) TSCC          ! tmp SCC
        CHARACTER(SICLEN3) CSIC          ! tmp char SIC
        CHARACTER(SICLEN3) CSIC2         ! tmp 2-digit SIC
        CHARACTER(SICLEN3) PCSIC         ! previous char SIC
        CHARACTER(SICLEN3) PCSIC2        ! previous 2-char SIC
        CHARACTER(MACLEN3) TMACT         ! tmp char MACT code
        CHARACTER(MACLEN3) PMACT         ! previous char MACT code
        CHARACTER(NAILEN3) TNAICS        ! tmp char NAICS code
        CHARACTER(NAILEN3) PNAICS        ! previous char NAICS code
        CHARACTER(BLRLEN3) BLID          ! tmp boiler ID
        CHARACTER(BLRLEN3) PBLID         ! previous boiler ID
        CHARACTER(ORSLEN3) CORS          ! tmp DOE plant ID
        CHARACTER(ORSLEN3) PCORS         ! previous DOE plant ID
        CHARACTER(OBRLEN3) PCORSBLR      ! previous DOE plant ID // boiler
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

C.............  Allocate memory for sorting index
            ALLOCATE( INDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )

C.............  Check if CIFIP is allocated.  
C.............  If it is, generate unique list of country/state/county codes
            IF( ASSOCIATED( CIFIP ) ) THEN

C.................  Count number of unique codes
                PCFIP = ' '
                J1 = 0
                DO S = 1, NSRC
 
                    CFIP = CIFIP( S )
                    IF( CFIP .NE. PCFIP ) J1 = J1 + 1
                    PCFIP = CFIP

                END DO
                NINVIFIP = J1

C.................  Allocate memory for country/state/county lists
                ALLOCATE( INVCFIP( NINVIFIP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVCFIP', PROGNAME )
                INVCFIP = ''

C.................  Create unique country/state/county codes list
                PCFIP = ' '
                J1 = 0
                DO S = 1, NSRC
 
                    CFIP = CIFIP( S )
                    IF( CFIP .NE. PCFIP ) THEN
                        J1 = J1 + 1
                        INVCFIP( J1 ) = CFIP
                        PCFIP = CFIP
                    END IF

                END DO

            END IF   ! End of CIFIP allocated or not

C.............  Check if CVTYPE is allocated
C.............  If it is, generate unique list of vehicle types
            IF( ALLOCATED( CVTYPE ) ) THEN
            
C.................  Initialize vehicle type sorting index
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO
                
C.................  Sort vehicle types
                CALL SORTIC( NSRC, INDX, CVTYPE )
                
C.................  Count number of unique vehicle types
                PVTYP = '-9'
                J1 = 0
                DO S = 1, NSRC
                    J = INDX( S )
                    TVTYP = CVTYPE( J )
                    
                    IF( TVTYP /= PVTYP ) J1 = J1 + 1
                    
                    PVTYP = TVTYP
                END DO
                NINVVTYP = J1

C.................  Allocate memory for vehicle type list
                ALLOCATE( INVVTYP( NINVVTYP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVVTYP', PROGNAME )
                
C.................  Create list of unique vehicle types
                PVTYP = '-9'
                J1 = 0
                DO S = 1, NSRC
                    J = INDX( S )
                    TVTYP = CVTYPE( J )
                    
                    IF( TVTYP /= PVTYP ) THEN
                        J1 = J1 + 1
                        INVVTYP( J1 ) = TVTYP
                        PVTYP = TVTYP
                    END IF
                END DO
            END IF  ! end vehicle type processing

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
                J1 = 0        ! to skip TSCC = 0
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

C.....................  Don't include expanded SCCs in left SCC list
                    IF( CHKEXPSCC( TSCC ) ) CYCLE

                    IF( SCCL .NE. PSCCL ) THEN
                        J2 = J2 + 1
                        INVSCL( J2 ) = SCCL
                        PSCCL = SCCL
                    END IF

                END DO 

            END IF   ! End SCC processing
            NINVSCL =  J2

C.............  Check if CISIC is allocated.  
C.............  If it is, generate unique list of SICs
            IF( ASSOCIATED( CISIC ) ) THEN

C.................  Initialize SIC sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO

C.................  Sort all SICs in the inventory in increasing order
                CALL SORTIC( NSRC, INDX, CISIC )

C.................  Count number of unique SICs
                PCSIC = '-9'
                J1 = 0
                DO S = 1, NSRC

                    J = INDX( S )
                    
                    CSIC = CISIC( J )
                    
                    IF( CSIC .NE. PCSIC ) J1 = J1 + 1

                    PCSIC = CSIC

                END DO 
                NINVSIC = J1

C.................  Allocate memory for SIC lists
                ALLOCATE( INVSIC( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSIC', PROGNAME )
                ALLOCATE( INVSIC2( NINVSIC ), STAT=IOS )
                CALL CHECKMEM( IOS,' INVSIC2', PROGNAME )

C.................  Create unique SIC lists
                PCSIC = '-9'
                PCSIC2 = '-9'
                J1 = 0
                J2 = 0
                DO S = 1, NSRC

                    J = INDX( S )

                    CSIC = CISIC( J )
                    CSIC2 = CSIC( 1:SICLEN3-2 )

                    IF( CSIC .NE. PCSIC ) THEN
                        J1 = J1 + 1
                        INVSIC( J1 ) = CSIC
                        PCSIC = CSIC
                    END IF

C.....................  Don't include expanded SICs in 2-digit SIC list
                    IF( CHKEXPSIC( CSIC ) ) CYCLE

                    IF( CSIC2 .NE. PCSIC2 ) THEN
                        J2 = J2 + 1
                        INVSIC2( J2 ) = CSIC2
                        PCSIC2 = CSIC2
                    END IF

               END DO 

            END IF   ! End SIC processing
            NINVSIC2 = J2

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

C.............  Check if CNAICS is allocated.  
C.............  If it is, generate unique list of NAICS codes
            IF( ASSOCIATED( CNAICS ) ) THEN

C.................  Initialize NAICS sorting index     
                DO S = 1, NSRC
                    INDX( S ) = S
                END DO                

C.................  Sort all NAICS in the inventory in increasing order
                CALL SORTIC( NSRC, INDX, CNAICS )

C.................  Count number of unique NAICS
                PNAICS = '-9'
                J1 = 0
                DO S = 1, NSRC
                
                    J = INDX( S )
                    
                    TNAICS = CNAICS( J )
                    
                    IF( TNAICS /= PNAICS ) J1 = J1 + 1
                    
                    PNAICS = TNAICS
                    
                END DO
                NINVNAICS = J1

C.................  Allocate memory for NAICS lists
                ALLOCATE( INVNAICS( NINVNAICS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVNAICS', PROGNAME )

C.................  Create unique NAICS list
                PNAICS = '-9'
                J1 = 0
                DO S = 1, NSRC
                
                    J = INDX( S )
                    
                    TNAICS = CNAICS( J )
                    
                    IF( TNAICS /= PNAICS ) THEN
                        J1 = J1 + 1
                        INVNAICS( J1 ) = TNAICS
                        PNAICS = TNAICS
                    END IF
                    
                END DO
                
            END IF    ! End NAICS processing
                
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
            NOBLRSRC = 0
            PCORS    = ' '
            PBLID    = ' '
            DO S = 1, NSRC

                CORS = CORIS ( S )
                BLID = CBLRID( S )
 
C.................  Skip missing ORIS IDs
                IF ( CORS == ' ' ) CYCLE

C.................  Count unique ORIS IDs; ORIS IDs are not guaranteed to be in
C                   blocks, but this will count the maximum number
                IF ( CORS /= PCORS ) THEN
                    NINVORIS = NINVORIS + 1
                END IF

C.................  Skip missing boiler IDs
                IF ( BLID == ' ' ) CYCLE
                
                NOBLRSRC = NOBLRSRC + 1

C.................  Count unique ORIS // boiler combos
                IF ( CORS /= PCORS .OR.  
     &               BLID /= PBLID      ) THEN
                    NORISBLR = NORISBLR + 1
                END IF

                PCORS = CORS
                PBLID = BLID

            END DO

C.............  Allocate memory for unsorted ORIS lists
            ALLOCATE( FOIDXA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FOIDXA', PROGNAME )
            ALLOCATE( INVORISA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORISA', PROGNAME )
            ALLOCATE( INVORFPA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORFPA', PROGNAME )
            ALLOCATE( INVODSCA( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVODSCA', PROGNAME )
            INVODSCA = ' '

C.............  Allocate memory for unsorted ORIS/boiler lists
            ALLOCATE( OBIDXA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBIDXA', PROGNAME )
            ALLOCATE( ORISBLRA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ORISBLRA', PROGNAME )
            ALLOCATE( OBSRCNTA( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCNTA', PROGNAME )

C.............  Store arrays
            NINVORIS = 0
            NORISBLR = 0
            PCORS    = ' '
            PBLID    = ' '
            DO S = 1, NSRC

                CFIP = CIFIP ( S )
                CORS = CORIS ( S )
                BLID = CBLRID( S )
                PDSC = CPDESC( S )

C.................  Skip missing ORIS IDs
                IF ( CORS == ' ' ) CYCLE

C.................  Unsorted ORIS arrays
                IF ( CORS /= PCORS ) THEN
                    J = INDEX1( CORS, NINVORIS, INVORISA )
                    
                    IF( J <= 0 ) THEN
                        NINVORIS = NINVORIS + 1
                        FOIDXA  ( NINVORIS ) = NINVORIS
                        INVORISA( NINVORIS ) = CORS
                        INVORFPA( NINVORIS ) = CFIP
                        INVODSCA( NINVORIS ) = PDSC
                    ELSE
                        IF( INVORFPA( J ) /= CFIP ) THEN
                       	    MESG = 'WARNING: Different FIPS codes ' //
     &                            'found for ORIS ID ' // CORS
     &                            // '.  Will use ' // CFIP //
     &                            ' for reporting.'
                            CALL M3MESG( MESG )
                        END IF
                        
                        IF( INVODSCA( J ) /= PDSC ) THEN
                            MESG = 'WARNING: Different plant ' //
     &                        'descriptions found for ORIS ID ' // CORS
     &                        // '.  Will use ' // INVODSCA( J ) //
     &                        ' for reporting.'
                            CALL M3MESG( MESG )
                        END IF
                    END IF
                END IF

C.................  Skip missing boiler IDs
                IF ( BLID == ' ' ) CYCLE

C.................  Unsorted oris/boiler array
                J = INDEX1( CORS // BLID, NORISBLR, ORISBLRA )
                    
                IF( J <= 0 ) THEN
                    NORISBLR = NORISBLR + 1
                    OBIDXA  ( NORISBLR ) = NORISBLR
                    ORISBLRA( NORISBLR ) = CORS // BLID
                    OBSRCNTA( NORISBLR ) = 1
                ELSE
                    OBSRCNTA( J ) = OBSRCNTA( J ) + 1
                END IF

                PCORS = CORS
                PBLID = BLID

            END DO
            
            IF( EFLAG ) THEN
                MESG = 'Problem counting ORIS IDs'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Sort arrays
            CALL SORTIC( NINVORIS, FOIDXA, INVORISA )
            CALL SORTIC( NORISBLR, OBIDXA, ORISBLRA )

C.............  Allocate memory for sorted ORIS lists
            ALLOCATE( INVORIS( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORIS', PROGNAME )
            ALLOCATE( INVORFP( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVORFP', PROGNAME )
            ALLOCATE( INVODSC( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVODSC', PROGNAME )
            ALLOCATE( IORSMTCH( NINVORIS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IORSMTCH', PROGNAME )
            IORSMTCH = .FALSE.

C.............  Allocate memory for sorted ORIS/boiler lists
            ALLOCATE( ORISBLR( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ORISBLR', PROGNAME )
            ALLOCATE( OBSRCNT( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCNT', PROGNAME )
            ALLOCATE( OBSRCBG( NORISBLR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCBG', PROGNAME )

C.............  Store sorted arrays
            DO I = 1, NINVORIS
                J = FOIDXA( I )
                INVORIS( I ) = INVORISA( J )
                INVORFP( I ) = INVORFPA( J )
                INVODSC( I ) = INVODSCA( J )
            END DO
            
            DO I = 1, NORISBLR
                J = OBIDXA( I )
                ORISBLR( I ) = ORISBLRA( J )
                OBSRCNT( I ) = OBSRCNTA( J )
            END DO

C.............  Build ORIS/boiler source start array
            OBSRCBG( 1 ) = 1
            DO I = 2, NORISBLR
                OBSRCBG( I ) = OBSRCBG( I-1 ) + OBSRCNT( I-1 )
            END DO

C.............  Create list of sources for each ORIS/boiler combo
            ALLOCATE( OBSRCNM( NOBLRSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OBSRCNM', PROGNAME )
            OBSRCNM = 0

            DO S = 1, NSRC

                CORS = CORIS ( S )
                BLID = CBLRID( S )

C.................  Skip missing ORIS or boiler IDs
                IF( CORS == ' ' .OR. BLID == ' ' ) CYCLE

C.................  Find combination in sorted list               
                J = FINDC( CORS // BLID, NORISBLR, ORISBLR )
                
                N = OBSRCBG( J )
                NS = N + OBSRCNT( J )
                DO
                    IF( OBSRCNM( N ) == 0 ) THEN
                        OBSRCNM( N ) = S
                        EXIT
                    ELSE
                        N = N + 1
                    END IF
                END DO
            
            END DO

C.............  Deallocate unneeded unsorted arrays
            DEALLOCATE( FOIDXA, OBIDXA, INVORFPA, OBSRCNTA, INVORISA,
     &                  ORISBLRA, INVODSCA )
                
            FIRSTORS = .FALSE.
        
        END IF   ! End boiler processing

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE GENUSLST
