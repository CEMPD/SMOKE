
        SUBROUTINE SETNONHAP( NRAWBP )

C**************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine integrates criteria and toxics pollutant 
C      emissions by creating NONHAP values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 11/02 by C. Seppanen
C      modified 1/2/06 by B. Baek
C
C**************************************************************************
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: POLVAL, IPOSCOD, NPCNT, CSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: INVSTAT, MXIDAT, INVDNAM, INVDVTS,
     &                      ITMSPC, ITEXPL, ITNAMA, NINVTBL

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: LNONHAP
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NPPOL, NCHARS, NEM, NDY

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  physical and mathematical constants
 
C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER(2)    CRLF
        INTEGER         INDEX1
        INTEGER         ENVINT
        LOGICAL         ENVYN
 
        EXTERNAL        CRLF, INDEX1, ENVINT, ENVYN
        
C.........  Pollutant names
        CHARACTER(11),      PARAMETER :: NONHAPCOMP  = 'NONHAP_TYPE'
        CHARACTER(6),       PARAMETER :: NONHAPREFIX = 'NONHAP'
        CHARACTER(4),       PARAMETER :: NOIEND = '_NOI'

C.........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (INOUT)      :: NRAWBP  ! no. raw records by pollutant

C.........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: VOCPOS( : )      ! location of VOC or TOG entry in srcs array
        INTEGER, ALLOCATABLE :: NVOCPOS( : )     ! location of VOC or TOG entry in srcs array for VOC_INV
        INTEGER, ALLOCATABLE :: VNMPOS( : )      ! position of VOC or TOG in poll name array

        INTEGER, ALLOCATABLE :: TMPIDX( : )      ! sorting index
        INTEGER, ALLOCATABLE :: TMPNPCNT( : )    ! number of pol per source
        INTEGER, ALLOCATABLE :: TMPPOSCOD( : )   ! tmp pollutant code array

        INTEGER, ALLOCATABLE :: NCRNOTOX( : )    ! number of sources with criteria but no toxics
        INTEGER, ALLOCATABLE :: NTOXNOCR( : )    ! number of sources with toxics but no criteria

        REAL,    ALLOCATABLE :: HAPEANN( : )     ! summed annual VOC/TOG HAPs emissions
        REAL,    ALLOCATABLE :: HAPEDAY( : )     ! summed average day VOC/TOG HAPs emissions

        REAL,    ALLOCATABLE :: TMPVAL( : )      ! tmp emissions array
        REAL,    ALLOCATABLE :: TMPPOLVAL( :,: ) ! tmp emissions array

        CHARACTER(IOVLEN3),ALLOCATABLE :: NONHAPNAM( : ) ! NONHAPS pol name array
        CHARACTER(3)      ,ALLOCATABLE :: NONHAPMOD( : ) ! NONHAPS modes array

C.........   Other local variables
        INTEGER  H,I,II,J,K,L,LL,L2,M,N,S   ! counters and indices
        
        INTEGER  IDXSIZE      ! size of sorting index
        INTEGER  IOS          ! I/O status
        INTEGER  CPOL         ! current pollutant number
        INTEGER  PPOL         ! previous pollutant number
        INTEGER  CPOLRAW      ! pollutant number in raw arrays
        INTEGER  CURRPOS      ! current position in POLVAL and IPOSCOD arrays
        INTEGER  NEWSRCPOS    ! current position in POLVAL and IPOSCOD arrays with VOC_INV poll
        INTEGER  NHAP         ! number of pollutants that needs NONHAP calculation in src array
        INTEGER  NNONHAP      ! number of pollutants that needs NONHAP calculation
        INTEGER  NONVPOS      ! position of NONHAP[VOC|TOG] in poll name array
        INTEGER  NVPOS        ! tmp location of VOC or TOG entry in srcs array for VOC_INV
        INTEGER  VPOS         ! tmp location of VOC or TOG entry in srcs array
        INTEGER  TMPPOS       ! tmp location of VOC or TOG entry in srcs array
        INTEGER  TOXPOS       ! position of NOI toxic in poll name array
        INTEGER  MXWARN       ! maximum number of warnings
        INTEGER  STIDX        ! starting index into pollutant array for current source
        INTEGER  ENDIDX       ! ending index into pollutant array
        INTEGER::NWARN = 0    ! current number of warnings
        
        LOGICAL  FNDPOL       ! true: found toxic pollutant to be processed
        LOGICAL  FNDVOC       ! true: found VOC entry in inventory
        LOGICAL  NFLAG        ! true: start non-HAP calculation
        LOGICAL  MFLAG        ! true: treat all sources integrated
        LOGICAL  EFLAG        ! true: error occured
        LOGICAL  LASTFLAG     ! true: entry is last for current source
        LOGICAL  NEEDSORT     ! true: need to resort pollutants for current source
        LOGICAL::PROCVOC=.FALSE. ! true: process VOC pollutants
        LOGICAL::PROCTOG=.FALSE. ! true: process TOG pollutants
        
        CHARACTER(3)       VOC_TOG  ! VOC or TOG
        CHARACTER(14)      NONHAP   ! name of NONHAP[VOC|TOG]
        CHARACTER(IOVLEN3) POLNAM   ! temporary pollutant name
        CHARACTER(IOVLEN3) CURPOL   ! temporary current pollutant name
        CHARACTER(IOVLEN3) TMPPOL   ! temporary current pollutant name
        CHARACTER(IOVLEN3) VOCPOL   ! temporary current pollutant name
        CHARACTER(300)     BUFFER   ! message buffer
        CHARACTER(256)     MESG     ! message buffer 
        
        CHARACTER(10) :: PROGNAME = 'SETNONHAP' ! program name

C***********************************************************************
C   begin body of subroutine SETNONHAP

C.........  Check setting of NONHAP[VOC|TOG] computation
        NFLAG = ENVYN ( 'SMK_NHAPEXCLUDE_YN', ' ', .FALSE., IOS )

C.........  Override NFLAG to treat all the sources as integrated sources.
        MFLAG = ENVYN ( 'SMK_INTEGRATE_ALL_YN', ' ', .FALSE., IOS )
        IF( MFLAG ) THEN
            NFLAG = .TRUE.
C.............  Allocate memory for source-based array from MODLISTS
            ALLOCATE( LNONHAP( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LNONHAP', PROGNAME )
            LNONHAP = .TRUE.   ! array
        END IF
        
        IF( .NOT. NFLAG ) RETURN

        CALL ENVSTR( NONHAPCOMP, MESG, ' ', VOC_TOG, IOS )
        MESG ='Processing NONHAP'// VOC_TOG // ' calculation.' 
        CALL M3MSG2( MESG )

        IF( IOS .NE. 0 ) THEN
            MESG = 'ERROR: NONHAP_TYPE environment variable ' //
     &             'is not defined for NONHAPVOC calculation'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( VOC_TOG /= 'VOC' .AND. VOC_TOG /= 'TOG' ) THEN
            MESG = 'ERROR: NONHAP_TYPE environment variable ' //
     &             'should be set to either VOC or TOG only.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
            
C.........  Set target processing NONHAP pollutant name (NONHAPVOC or NONHAPTOG)
        NONHAP = NONHAPREFIX // VOC_TOG

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS ) 
        NWARN = 0

C.........  Define size of arrays for NONHAP calculation
        NNONHAP = 0
        DO I = 1, MXIDAT

C.............  Counts total numbers of names for mobile emission processing modes
            IF( INVSTAT( I ) < 0 ) CYCLE

            L  = INDEX( INVDNAM( I ) , ETJOIN ) 

            IF( L > 0 ) THEN
                LL = LEN( INVDNAM( I ) )

                IF( INVDNAM( I )( L+2:LL ) == VOC_TOG ) THEN
                    NNONHAP = NNONHAP + 1
                END IF

C.............  Counts total numbers of pollutatns that need NONHAP calculation for VOC and TOG
            ELSE IF( INVDNAM( I ) == VOC_TOG ) THEN
                NNONHAP = NNONHAP + 1
 
            END IF

        END DO

C.........  Check that either VOC or TOG is in the inventory
        IF( NNONHAP == 0 ) THEN
            MESG = 'WARNING: No ' // NONHAP // ' calculattion was ' //
     &             ' processed because there are no ' // VOC_TOG //
     &             ' & ' // NONHAP // ' found in the INVTABLE file.'//
     &             CRLF() // BLANK5 // ':: Skipping ' // NONHAP // 
     &            ' calculation and NHAPEXCLUDE processing steps.'
            CALL M3MSG2( MESG )
            RETURN
        END IF

C.........  Determine maximum size for allocating memory
        ALLOCATE( VNMPOS( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNMPOS', PROGNAME )
        ALLOCATE( NONHAPNAM( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NONHAPNAM', PROGNAME )
        ALLOCATE( NONHAPMOD( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NONHAPMOD', PROGNAME )
        ALLOCATE( VOCPOS( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VOCPOS', PROGNAME )
        ALLOCATE( NVOCPOS( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NVOCPOS', PROGNAME )
        ALLOCATE( HAPEANN( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HAPEANN', PROGNAME )
        ALLOCATE( HAPEDAY( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HAPEDAY', PROGNAME )
        ALLOCATE( NCRNOTOX( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCRNOTOX', PROGNAME )
        ALLOCATE( NTOXNOCR( NNONHAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NTOXNOCR', PROGNAME )

        VNMPOS = 0
        NONHAPNAM = ' '
        NONHAPMOD = ' '
        NCRNOTOX = 0
        NTOXNOCR = 0 

C.........  Max number of souces with pol/act including a new species (VOC_INV)
        NRAWBP = NRAWBP + NSRC

C.........  Allocate memory for sorted inventory data
        ALLOCATE( TMPNPCNT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPNPCNT', PROGNAME )
        ALLOCATE( TMPPOSCOD( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOSCOD', PROGNAME )
        ALLOCATE( TMPPOLVAL( NRAWBP,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOLVAL', PROGNAME )
        ALLOCATE( TMPVAL( NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPVAL', PROGNAME )
        
        TMPNPCNT  = 0
        TMPPOSCOD = 0
        TMPPOLVAL = BADVAL3
        TMPVAL    = BADVAL3

C.........  Store the poisiont of pollutants for HAPs and NONHAPs
        NHAP = 0
        DO I = 1, MXIDAT

C.............  Skip if pol is activity pol as of VMT
            IF( INVSTAT( I ) < 0 ) CYCLE

            L  = INDEX( INVDNAM( I ) , ETJOIN ) 

C.............  Store positions of names for mobile emission processing modes
            IF( L > 0 ) THEN
                LL = LEN( INVDNAM( I ) )

                IF( INVDNAM( I )( L+2:LL ) == VOC_TOG ) THEN
                    NHAP = NHAP + 1
                    VNMPOS( NHAP ) = I
                    NONHAPMOD( NHAP ) = INVDNAM( I )( 1:L-1 )
                    NONHAPNAM( NHAP ) = INVDNAM( I )( 1:L+1 ) //
     &                             NONHAPREFIX // INVDNAM( I )( L+2:LL )

C.....................  look for NONHAP pol names in a list of pols
                    NONVPOS=INDEX1( NONHAPNAM( NHAP ), MXIDAT, INVDNAM )
                    IF( NONVPOS < 1 ) THEN
                        MESG = 'ERROR: Pollutant ' // TRIM
     &                         ( NONHAPNAM( NHAP ))//' is not found '
     &                         // 'in the INVTABLE file.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

C.............  Store position of pollutatns that need NONHAP calculation for VOC and TOG
            ELSE IF( INVDNAM( I ) == VOC_TOG ) THEN
                NHAP = NHAP + 1
                VNMPOS( NHAP ) = I     ! position of VOC or TOG pol in a list of pols
                NONHAPNAM( NHAP ) = NONHAPREFIX // INVDNAM( I )

C.................  Look for NONHAP pol names in a list of pols
                NONVPOS = INDEX1( NONHAPNAM( NHAP ), MXIDAT, INVDNAM )
                IF( NONVPOS < 1 ) THEN
                    MESG = 'ERROR: Pollutant ' // 
     &                     TRIM( NONHAPNAM( NHAP  ))//' is not found '
     &                     // 'in the INVTABLE file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

       END DO     ! eloop of a list of inv pols

C.........  Set up logical flags based on found pollutants (VOC or TOG)
        PROCVOC = .TRUE.
        IF( VOC_TOG == 'TOG' ) THEN
            PROCTOG = .TRUE.
            PROCVOC = .FALSE.
        END IF

C.........  Check if any toxics pollutants are to be subtracted
        FNDPOL = .FALSE.
        DO I = 1, MXIDAT
            IF( INVDVTS( I ) /= 'N' ) THEN
                IF( .NOT. PROCTOG .AND. INVDVTS( I ) == 'T' ) CYCLE
                IF( INVSTAT( I ) /= 0 ) THEN
                    FNDPOL = .TRUE.
                    EXIT
                END IF
            END IF
        END DO

        IF( .NOT. FNDPOL ) THEN
            MESG = 'WARNING: No pollutants were selected from ' // 
     &             'the inventory as being subratactable for ' //
     &             'the NHAPEXCLUDE processing.'
            CALL M3MSG2( MESG )
            RETURN
        END IF

        CURRPOS = 0
        NEWSRCPOS = 0

C.........  Loop through sources
        DO I = 1, NSRC

C.............  Initialize values for this source
            TMPPOS = 0
            TMPPOL = ' '
            TMPVAL = 0.0
            VOCPOS = 0
            NVOCPOS = 0
            FNDVOC = .FALSE.            
            HAPEANN = 0.
            HAPEDAY = 0.
            
            EFLAG = .FALSE.
            LASTFLAG = .FALSE.

            TMPNPCNT( I ) = NPCNT( I )

C.............  Process source if it is not integrated
            IF( ALLOCATED( LNONHAP ) ) THEN
            IF( .NOT. LNONHAP( I ) ) THEN

C.................  Loop through pollutants for this source
                DO J = 1, NPCNT( I )

C.....................  Increment current position in arrays            
                    CURRPOS = CURRPOS + 1
                    NEWSRCPOS = NEWSRCPOS + 1
                    
C.....................  Store pollutant for this position
                    CPOL = IPOSCOD( CURRPOS )

C.....................  Store to tmp arrays of IPOSCOD and POLVAL
                    TMPPOSCOD( NEWSRCPOS )   = IPOSCOD( CURRPOS )
                    TMPPOLVAL( NEWSRCPOS,: ) = POLVAL ( CURRPOS,: )

C.....................  If current pollutant is VOC(_INV) entry, save position and emissions
                    DO H = 1, NNONHAP
                        IF( CPOL == VNMPOS( H ) ) THEN
                            TMPPOS = CURRPOS
                            TMPPOL = INVDNAM( CPOL )
                            TMPVAL( : ) = POLVAL( CURRPOS,: )
                        END IF
                    END DO

C.....................  If pollutant is not part of VOC or TOG, cycle
                    IF( INVDVTS( CPOL ) == 'N' ) CYCLE

C.....................  Find pollutant position in raw list
                    CPOLRAW = INDEX1( INVDNAM( CPOL ), NINVTBL, ITNAMA )

C.....................  If pollutant is not a model species, set it to zero
                    IF( .NOT. ITMSPC( CPOLRAW ) ) THEN
                        POLVAL( CURRPOS,NEM ) = 0.0
                        POLVAL( CURRPOS,NDY ) = 0.0

C..................... Otherwise, if pollutant is not an explicit species, rename to NOI
                    ELSE IF( .NOT. ITEXPL( CPOLRAW ) ) THEN
                
C.........................  Create NOI name
                        POLNAM = INVDNAM( CPOL )
                        IF( LEN_TRIM( POLNAM ) > 11 ) THEN
                            POLNAM = POLNAM(1:12)
                        END IF
                        POLNAM = TRIM( POLNAM ) // NOIEND

C.........................  Find NOI name in pollutant names array                        
                        TOXPOS = INDEX1( POLNAM, MXIDAT, INVDNAM )

C.........................  If found, set pollutant for current source to NOI pollutant
                        IF( TOXPOS > 0 ) THEN
                            IPOSCOD( CURRPOS ) = TOXPOS
                            INVSTAT( TOXPOS ) = 2
                        ELSE
                            MESG = 'ERROR: Non-integrated toxic ' //
     &                             'pollutant ' // TRIM( POLNAM ) //
     &                             ' was not found in the ' //
     &                             'INVTABLE file.'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF
                     
                    END IF  ! check if pollutant is model or explicit species

                END DO  ! loop through pollutants
                
C.................  Skip rest of loop and add new poll VOC_INV since we're done with this source
                GOTO 999    ! Add new poll VOC_INV
                
            END IF
            END IF

C.............  Loop through pollutants for this source            
            DO J = 1, NPCNT( I )

C.................  Increment current position in arrays            
                CURRPOS   = CURRPOS + 1
                NEWSRCPOS = NEWSRCPOS + 1

C.................  Store pollutant for this position
                CPOL   = IPOSCOD( CURRPOS )
                POLNAM = INVDNAM( CPOL )

C.................  Store to tmp arrays of IPOSCOD and POLVAL
                TMPPOSCOD( NEWSRCPOS )   = IPOSCOD( CURRPOS )
                TMPPOLVAL( NEWSRCPOS,: ) = POLVAL ( CURRPOS,: )

C.................  If current pollutant is VOC(_INV) entry, save position and emissions
                DO H = 1, NNONHAP
                    IF( CPOL == VNMPOS( H ) ) THEN
                        TMPPOS = CURRPOS
                        TMPPOL = POLNAM
                        TMPVAL( : ) = POLVAL( CURRPOS,: )
                        VOCPOS( H )  = CURRPOS
                        NVOCPOS( H ) = NEWSRCPOS
                    END IF
                END DO

C.................  Check if this is last pollutant for this source
                IF( J == NPCNT( I ) ) THEN
                    LASTFLAG = .TRUE.

C.................  Otherwise, if not part of VOC or TOG, skip rest of loop
                ELSE IF( INVDVTS( CPOL ) == 'N' ) THEN
                    CYCLE
                    
                END IF

C.................  Check emission processing modes
                L = INDEX( POLNAM , ETJOIN )
                IF( L > 0 ) THEN
                    M = INDEX1( POLNAM( 1:L-1 ), NNONHAP, NONHAPMOD )
                ELSE
                    M = INDEX1( ' ', NNONHAP, NONHAPMOD )
C                    M = 1
                END IF

C.................  Sum toxic emissions for this source
C                   INVDVTS = 'V' => part of VOC and TOG
C                   INVDVTS = 'T' => part of TOG only
                IF( PROCVOC .AND. INVDVTS( CPOL ) == 'V' ) THEN
                    HAPEANN( M ) = HAPEANN( M ) + POLVAL( CURRPOS,NEM )
                    HAPEDAY( M ) = HAPEDAY( M ) + POLVAL( CURRPOS,NDY )
                    FNDVOC = .TRUE.
                END IF

                IF( PROCTOG ) THEN
                    HAPEANN( M ) = HAPEANN( M ) + POLVAL( CURRPOS,NEM )
                    HAPEDAY( M ) = HAPEDAY( M ) + POLVAL( CURRPOS,NDY )
                    FNDVOC = .TRUE.

                END IF

C.................  If this is not the last entry for source, cycle
C                   Otherwise, check conditions and subtract toxic
C                   emissions from criteria values
                IF( .NOT. LASTFLAG ) CYCLE

                LASTFLAG = .FALSE.

C.................  Loop through pollutants that need NONHAP calculation
                DO K = 1, NNONHAP

C.....................  Skip NONHAP calculation when summed HAPs is zero
                    IF( HAPEANN( K ) == 0.0 .AND. FNDVOC ) CYCLE
                    IF( HAPEDAY( K ) == 0.0 .AND. FNDVOC ) CYCLE

C.....................  Restore position of VOC or TOG pollutants
                    VPOS   = VOCPOS( K )
                    NVPOS  = NVOCPOS( K )
                    CURPOL = INVDNAM( VNMPOS( K ) )   ! VOC, EXH__VOC,,,

C....................  Check if this source had no criteria and no toxics
C                      This happens for activity only sources
                    IF( VPOS == 0 .AND. .NOT. FNDVOC ) CYCLE

C.....................  Format information for this source
                    CALL FMTCSRC( CSOURC( I ), NCHARS, BUFFER, L2 )

C.....................  Give warning if source has toxics but no criteria
                    IF( VPOS == 0 .AND. FNDVOC ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Found toxic emissions ' //
     &                             'but no ' // TRIM( CURPOL ) //
     &                             ' emissions for source:' // 
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        NTOXNOCR( K ) = NTOXNOCR( K ) + 1
                        EFLAG = .TRUE.
                    END IF

C.....................  Give warning if source has criteria but no toxics
                    IF( VPOS /= 0 .AND. .NOT. FNDVOC ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Found ' // TRIM( CURPOL )// 
     &                             ' emissions but no toxic '//
     &                             'emissions for source:' // 
     &                            CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        NCRNOTOX( K ) = NCRNOTOX( K ) + 1
                        EFLAG = .TRUE.
                    END IF

C.....................  Skip rest of loop if an error has occured
                    IF( EFLAG ) CYCLE

C.....................  Subtract toxic emissions from criteria emissions  
                    POLVAL( VPOS,NEM ) = 
     &                      POLVAL( VPOS,NEM ) - HAPEANN( K )
                    POLVAL( VPOS,NDY ) =
     &                      POLVAL( VPOS,NDY ) - HAPEDAY( K )

C.....................  Check that annual NONHAP value is not negative
                    IF( POLVAL( VPOS,NEM ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total annual toxic ' //
     &                         'emissions greater than annual ' //
     &                         TRIM( CURPOL )// ' emissions for source:'
     &                         // CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        POLVAL( VPOS,NEM ) = .0
                    END IF

C....................  Check that average day NONHAP value is not negative
                    IF( POLVAL( VPOS,NDY ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total average day ' //
     &                             'toxic emissions greater than ' //
     &                             'average day ' // TRIM( CURPOL ) // 
     &                             ' emissions for source:' // CRLF() //
     &                             BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF

                        POLVAL( VPOS,NDY ) = .0
                    END IF

C.....................  Rename VOC to NONHAPVOC
C.....................  increment a number of poll for NONHAP[VOC|TOG]
                    NONVPOS = INDEX1( NONHAPNAM( K ), MXIDAT, INVDNAM )
                    IPOSCOD  ( VPOS )   = NONVPOS
                    TMPPOSCOD( NVPOS )  = NONVPOS
                    TMPPOLVAL( NVPOS,: )= POLVAL( VPOS,: )
                    INVSTAT  ( NONVPOS )= 2

                END DO  !loop through pollutants that need NONHAP calculation

            END DO  ! loop through pollutants

C.............  Add a new poll VOC_INV
999         IF( TMPPOS > 0 ) THEN
                TMPNPCNT( I )  = NPCNT( I ) + 1
                VOCPOL = TRIM( TMPPOL ) // '_INV'
                NONVPOS = INDEX1( VOCPOL,MXIDAT,INVDNAM )

C.................  Check to see a new pol VOC_INV listed in INVTABLE file
                IF( NONVPOS < 1 ) THEN
                     MESG = 'ERROR: Pollutant '//TRIM( VOCPOL )
     &                       //' is not found in the INVTABLE file.'
                     CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NEWSRCPOS = NEWSRCPOS + 1
                TMPPOSCOD( NEWSRCPOS )   = NONVPOS
                TMPPOLVAL( NEWSRCPOS,: ) = TMPVAL ( : )
                INVSTAT  ( NONVPOS )     = 2

            END IF

        END DO  ! loop through sources

        NRAWBP = NEWSRCPOS

C.........  Update NPCNT, IPOSCOD and POLVAL including VOC_INV to retain org VOC
        DEALLOCATE( IPOSCOD, POLVAL )
        
        ALLOCATE( POLVAL( NRAWBP,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
        ALLOCATE( IPOSCOD( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )
        
        POLVAL = 0.0
        IPOSCOD = 0

C.........  Reinitialize arrays with updated arrays including VOC_INV
        CURRPOS = 0
        DO S = 1, NSRC

            NPCNT( S )   = TMPNPCNT( S )
 
            DO J = 1, NPCNT( S )
C................  Increment current position in arrays            
               CURRPOS   = CURRPOS + 1

C................  Store pollutant for this position
               IPOSCOD( CURRPOS ) = TMPPOSCOD( CURRPOS )
               POLVAL( CURRPOS,:) = TMPPOLVAL( CURRPOS,: )
               POLNAM = INVDNAM( IPOSCOD( CURRPOS) )
            ENDDO

        ENDDO
        
        DEALLOCATE( TMPNPCNT, TMPPOLVAL, TMPPOSCOD )

C.........  Loop through pollutants that need NONHAP calculation
        N = 0
        DO K = 1, NNONHAP
            IF( NCRNOTOX( K ) > 0 .OR. NTOXNOCR( K ) > 0 ) THEN
                N = N + 1
                MESG = 'During processing, the following number ' //
     &                 'of sources were encountered: '
                IF( N == 1 ) CALL M3MESG( MESG )

                CURPOL = INVDNAM( VNMPOS( K ) )   ! VOC, EXH__VOC,,,

                IF( NCRNOTOX( K ) > 0 ) THEN
                    WRITE( MESG,94010 ) BLANK5 // 'Sources with ' // 
     &                  TRIM( CURPOL ) // ' emissions but no toxic ' //
     &                  'emissions: ', NCRNOTOX( K )
                    CALL M3MESG( MESG )
                END IF

                IF( NTOXNOCR( K ) > 0 ) THEN
                    WRITE( MESG,94010 ) BLANK5 // 'Sources with toxic'//
     &                  ' emissions but no ' // TRIM( CURPOL ) //
     &                  ' emissions: ', NTOXNOCR( K )
                    CALL M3MESG( MESG )
                END IF

            END IF

        END DO

C.........  Sort POLVAL and IPOSCOD to put new pollutants (NONHAP and NOI) in correct order

C.........  Determine maximum size for sorting index, then allocate memory
        IDXSIZE = MAXVAL( NPCNT )
        
        ALLOCATE( TMPIDX( IDXSIZE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPIDX', PROGNAME )
        ALLOCATE( TMPPOLVAL( IDXSIZE,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOLVAL', PROGNAME )
        ALLOCATE( TMPPOSCOD( IDXSIZE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOSCOD', PROGNAME )

        STIDX = 1
        
C.........  Loop through sources        
        DO I = 1, NSRC

C.............  Set ending index for this source
            ENDIDX = STIDX + NPCNT( I ) - 1

C.............  Initialize sorting index and check if any sorting needs to be done
            PPOL = 0
            NEEDSORT = .FALSE.

            DO J = 1, NPCNT( I )
                CPOL = IPOSCOD( STIDX + J - 1 )

C.................  If current pollutant is lower than previous, need to resort                
                IF( CPOL < PPOL ) THEN
                    NEEDSORT = .TRUE.
                END IF
                
                PPOL = CPOL
                TMPIDX( J ) = J
            END DO

C.............  Make sure this source needs to be resorted
            IF( NEEDSORT ) THEN

C.................  Store current values in temporary arrays
                TMPPOLVAL( 1:NPCNT( I ),: ) = POLVAL ( STIDX:ENDIDX,: )
                TMPPOSCOD( 1:NPCNT( I ) )   = IPOSCOD( STIDX:ENDIDX )
                
C.................  Sort section of IPOSCOD corresponding to this source
                CALL SORTI1( NPCNT( I ), TMPIDX( 1:NPCNT( I ) ), 
     &                       IPOSCOD( STIDX:ENDIDX ) ) 

C.................  Loop through pollutants for current source        
                DO J = 1, NPCNT( I )
            
                    K = TMPIDX( J )
                    
                    POLVAL ( STIDX + J - 1,: ) = TMPPOLVAL( K,: )
                    IPOSCOD( STIDX + J - 1 )   = TMPPOSCOD( K )
            
                END DO

            END IF  ! check if pollutants need sorting

C.............  Increment starting index        
            STIDX = ENDIDX + 1
        
        END DO

C.........  Deallocate local memory
        DEALLOCATE( TMPIDX, TMPPOLVAL, TMPPOSCOD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SETNONHAP
       
