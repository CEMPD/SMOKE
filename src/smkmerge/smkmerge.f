
        PROGRAM SMKMERGE

C***********************************************************************
C  program SMKMERGE body starts at line 148
C
C  DESCRIPTION:
C      The purpose of this program is to merge the inventory or hourly
C      emissions files from the Temporal program with gridding matrices and
C      with optionally any combination of speciation matrices and 3 control
C      matrices (different types).  The program can operate on from 1 to 4
C      source categories (area, biogenic, mobile, or point sources), or any
C      combination of these.  If a layer fractions file is input, then the
C      output file is 3-d.  This program is not used for the MPS/MEPSE files
C      for CMAQ.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY:
     &          AFLAG, BFLAG, MFLAG, PFLAG,                     ! source flags
     &          AUFLAG, MUFLAG, PUFLAG,                         ! mult control flags
     &          ARFLAG, MRFLAG, PRFLAG,                         ! reac control flags
     &          APRJFLAG, MPRJFLAG, PPRJFLAG, PFACFLAG,         ! growth flags (pfac variable)
     &          AFLAG_BD, MFLAG_BD, PFLAG_BD,                   ! by-day hourly emis flags
     &          TFLAG, SFLAG, LFLAG,                            ! use temporal, spec, layers
     &          PINGFLAG, ELEVFLAG, EXPLFLAG,                   ! ping, elevated, expl. plume
     &          INLINEFLAG, SRCGRPFLAG, SGDEV,                  ! inline, source groups
     &          SUBSECFLAG, NGRPS, SUBOUTNAME,                  ! sub-sector source groups
     &          LMKTPON, LREPANY,                               ! mkt penetration, any reports
     &          CDEV, EDEV, GDEV,                               ! costcy, elev/ping, grid surg
     &          AENAME, ATNAME, AGNAME, ASNAME, ARNAME, AUNAME, ! area files
     &          BTNAME,                                         ! biogenic files
     &          MENAME, MTNAME, MGNAME, MSNAME, MRNAME, MUNAME, ! mobile files
     &          PENAME, PTNAME, PGNAME, PSNAME, PRNAME, PUNAME, ! point files
     &          PLNAME, PVNAME, PHNAME,
     &          NASRC, NMSRC, NPSRC, EMLAYS,                    ! no. of srcs, no. emis layers
     &          ANMSPC, BNMSPC, MNMSPC, PNMSPC, NMSPC,          ! no. species
     &          ANGMAT, MNGMAT,                                 ! no. gridding matrix entries
     &          ANSREAC, MNSREAC, PNSREAC,                      ! no. src w/ reac controls
     &          ARNMSPC, MRNMSPC, PRNMSPC,                      ! no. reac species
     &          AEMNAM, BEMNAM, MEMNAM, PEMNAM, EMNAM, NSMATV,  ! species names and length of EMNAM
     &          ANMAP, AMAPNAM, AMAPFIL,                        ! area map file
     &          MNMAP, MMAPNAM, MMAPFIL,                        ! mobile map file
     &          PNMAP, PMAPNAM, PMAPFIL,                        ! point map file
     &          VGRPCNT, IDVGP, GVNAMES,                        ! group count, ids, var names
     &          SIINDEX, SPINDEX, GVLOUT,                       ! EANAM & EMNAM idx, output pts
     &          A_EXIST, M_EXIST, P_EXIST,                      ! grp indices for inv emis
     &          AU_EXIST, MU_EXIST, PU_EXIST,                   ! grp indices for mult controls
     &          AR_EXIST, MR_EXIST, PR_EXIST,                   ! grp indices for reac controls
     &          AS_EXIST, BS_EXIST, MS_EXIST, PS_EXIST,         ! grp indices for spec matrices
     &          SDATE, STIME, NSTEPS, TSTEP, PVSDATE, PVSTIME,  ! episode information
     &          ASDATE, MSDATE, PSDATE,                         ! dates for by-day hrly emis
     &          BIOGFAC, BIOTFAC, GRDFAC, TOTFAC,               ! conversion factors
     &          AEMSRC, MEMSRC, PEMSRC,                         ! inv or hrly emissions
     &          AEISRC, MEISRC, PEISRC,                         ! inv only emissions
     &          AGMATX, MGMATX, PGMATX,                         ! gridding matrices
     &          ASMATX, MSMATX, PSMATX,                         ! speciation matrices
     &          ARINFO, MRINFO, PRINFO,                         ! reactivity matrices
     &          AEMGRD, BEMGRD, MEMGRD, PEMGRD, TEMGRD,         ! gridded emissions
     &          AEBCNY, BEBCNY, MEBCNY, PEBCNY,                 ! cnty total spec emissions
     &          AEUCNY, MEUCNY, PEUCNY,                         ! cnty total mult control emis
     &          AERCNY, MERCNY, PERCNY,                         ! cnty total reac control emis
     &          AECCNY, MECCNY, PECCNY,                         ! cnty total all-control emis
     &          LFRAC, EANAM, TONAMES,                          ! layer frac, pol/act names
     &          ISRCGRP, EMGGRD                        ! emis by grid cell and src group

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: ACRIDX, ACRREPEM, ACRPRJFC, ACRMKTPN,
     &                      MCRIDX, MCRREPEM, MCRPRJFC, MCRMKTPN,
     &                      PCRIDX, PCRREPEM, PCRPRJFC, PCRMKTPN,
     &                      ACRFAC, MCRFAC, PCRFAC,
     &                      ACUMATX, MCUMATX, PCUMATX

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: INDXH, NHRSRC, GRPGID, ELEVFLTR, ELEVSRC,
     &                     GROUPID

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTY, AICNY, MICNY, PICNY

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGFIPS, SRGFIPS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, OFFLAG

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF
        CHARACTER(10)   HHMMSS
        INTEGER         INDEX1
        INTEGER         WKDAY
        LOGICAL         USEEXPGEO

        EXTERNAL    CRLF, HHMMSS, INDEX1, WKDAY, USEEXPGEO

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(16), PARAMETER :: PROGNAME = 'SMKMERGE' ! program name

        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local temporary array for input and output variable names
        CHARACTER(IOVLEN3), ALLOCATABLE :: VARNAMES( : )
        CHARACTER(IOVLEN3), ALLOCATABLE :: INNAMES ( : )

C...........   Local allocatable arrays for creating list of all explicit srcs
        INTEGER, ALLOCATABLE :: TMPSRC( : )
        INTEGER, ALLOCATABLE :: TMPIDX( : )
        LOGICAL, ALLOCATABLE :: SRCFLG( : )
        REAL,    ALLOCATABLE :: TCUMATX( : )    ! tmp array to store control/proj matrix

C...........    Local allocatable array for tracking whether species have
C               already been processed for elevated sources.
        LOGICAL, ALLOCATABLE :: ELEV_SPCSET( : )

C...........    Local variables for array sizes
        INTEGER         APOLSIZ ! work area inventory emissions array size
        INTEGER         MPOLSIZ ! work mobile inventory emissions array size
        INTEGER         PPOLSIZ ! work point inventory emissions array size
        INTEGER         ASPCSIZ ! work area speciation matrix array size
        INTEGER         MSPCSIZ ! work mobile speciation matrix array size
        INTEGER         PSPCSIZ ! work point speciation matrix array size
        INTEGER         AMULSIZ ! work area multipl control matrix array size
        INTEGER         MMULSIZ ! work mobile multipl control matrix array size
        INTEGER         PMULSIZ ! work point multipl control matrix array size

C...........   Logical names and unit numbers (not in MODMERGE)
        INTEGER         LDEV

C...........   Other local variables

        INTEGER          I, J, K, L1, L2, M, N, NG, V, S, T ! counters and indices

        INTEGER          AJDATE        ! area-source Julian date for by-day
        INTEGER          DAY           ! day-of-week index (monday=1)
        INTEGER          IOS           ! tmp I/O status
        INTEGER          JDATE         ! Julian date (YYYYDDD)
        INTEGER          JTIME         ! time (HHMMSS)
        INTEGER          K1            ! tmp index for valid ar spc matrix
        INTEGER          K2            ! tmp index for valid mb spc matrix
        INTEGER          K3            ! tmp index for valid pt spc matrix
        INTEGER          K4            ! tmp index for valid ar reactvty matrix
        INTEGER          K5            ! tmp index for valid mb reactvty matrix
        INTEGER          KA, KB, KM, KP! tmp index to src-category species
        INTEGER          LDATE         ! Julian date from previous iteration
        INTEGER          MJDATE        ! mobile-source Julian date for by-day
        INTEGER          MXGRP         ! max no. of variable groups
        INTEGER          MXVARPGP      ! max no. of variables per group
        INTEGER          NGRP          ! actual no. of pollutant groups
        INTEGER          NMAJOR        ! no. elevated sources
        INTEGER          NPING         ! no. plum-in-grid sources
        INTEGER          NVPGP         ! tmp actual no. variables per group
        INTEGER          NTSRC         ! tmp actual no. of sources
        INTEGER          OCNT          ! tmp count output variable names
        INTEGER       :: PDAY = 0      ! previous iteration day no.
        INTEGER          PGID          ! previous iteration group ID no.
        INTEGER          PJDATE        ! point-source Julian date for by-day
        INTEGER      :: SRGNROWS = 0   ! no. rows in surrogates file
        INTEGER      :: SRGNCOLS = 0   ! no. cols in surrogates file

        REAL             RDUM
        REAL             F1, F2, FB    ! tmp conversion factors

        LOGICAL      :: INITELEV = .TRUE.   ! true: reintialize ELEVEMIS array

        CHARACTER(16)      SRGFMT           ! gridding surrogates format
        CHARACTER(16)   :: SRGGRDNM  = ' '  !  surrogates file grid name
        CHARACTER(300)     MESG    ! message buffer
        CHARACTER(IOVLEN3) LBUF    ! previous species or pollutant name
        CHARACTER(IOVLEN3) PBUF    ! tmp pollutant or emission type name
        CHARACTER(IOVLEN3) SBUF    ! tmp species or pollutant name
        CHARACTER(PLSLEN3) VBUF    ! pol to species or pol description buffer

C***********************************************************************
C   begin body of program SMKMERGE

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Retrieve control environment variables and set logical control
C           flags. Use a local module to pass the control flags.
        CALL GETMRGEV

C.........  Open input files and retrieve episode information
        CALL OPENMRGIN( SRGNROWS, SRGNCOLS, SRGGRDNM, SRGFMT )

C.........  Do setup for biogenic state and county reporting or source
C           apportionment
        IF( BFLAG .AND. ( LREPANY .OR. SRCGRPFLAG .OR. SUBSECFLAG ) ) THEN

C.............  Read gridding surrogates
            CALL RDSRG( .FALSE., GDEV, SRGFMT, SRGNROWS, SRGNCOLS )

C.........  If output grid is different from surrogates, write message
            IF ( OFFLAG ) THEN
                L1 = LEN_TRIM( SRGGRDNM )
                MESG = 'NOTE: gridding surrogates (for biogenic '//
     &                 'totals) extracted for output'// CRLF()//
     &                 BLANK10 //'grid from grid "' //
     &                 SRGGRDNM( 1:L1 ) // '"'
                CALL M3MSG2( MESG )
            END IF

        END IF

C.........  Create arrays of sorted unique pol-to-species
C.........  Create arrays of sorted unique pollutants
C.........  Create arrays of sorted unique species
        CALL MRGVNAMS

C.........  Determine units conversion factors
        CALL MRGUNITS

C.........  Read in any needed source characteristics
        CALL RDMRGINV

C.........  Do setup for state and county reporting
C.........  Do this even if there LREPANY is false, in order to allocate
C           memory for the state and county total arrays to ensure
C           MRGMULT will work.

C.........  Read the state and county names file and store for the
C           states and counties in the grid
C.........  For anthropogenic source categories, use FIPS list
C           from the inventory for limiting state/county list

        IF( AFLAG .OR. MFLAG .OR. PFLAG ) THEN
            IF( USEEXPGEO() ) THEN
                CALL RDGEOCODES( NINVIFIP, INVCFIP )
            ELSE
                CALL RDSTCY( CDEV, NINVIFIP, INVCFIP )
            END IF

C.........  Otherwise, for biogenic merge only, use list of codes from the
C           surrogates file needed for state and county totals
        ELSE
            IF( USEEXPGEO() ) THEN
                CALL RDGEOCODES( NSRGFIPS, SRGFIPS )
            ELSE
                CALL RDSTCY( CDEV, NSRGFIPS, SRGFIPS )
            END IF

        END IF

C.........  Allocate memory for fixed-size arrays by source category...

        CALL ALLOCMRG( MXGRP, MXVARPGP, AMULSIZ, MMULSIZ, PMULSIZ,
     &                 ASPCSIZ, MSPCSIZ, PSPCSIZ, APOLSIZ, MPOLSIZ,
     &                 PPOLSIZ )

C.........  Read in elevated sources and plume-in-grid information, if needed
C.........  Reset flag for PinG if none in the input file
        IF( PFLAG .AND. (ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG) ) THEN

            CALL RDPELV( EDEV, NPSRC, ELEVFLAG, NMAJOR, NPING )

            IF( ELEVFLAG .AND. NMAJOR .EQ. 0 ) THEN
                MESG = 'WARNING: No sources are major elevated ' //
     &                 'sources in input file, ' // CRLF() //
     &                 BLANK10 // 'so elevated source emissions ' //
     &                 'file will not be written.'
                CALL M3MSG2( MESG )
                ELEVFLAG = .FALSE.
            ELSE IF ( NMAJOR .EQ. 0 ) THEN
                NMAJOR = NPSRC
            END IF

            IF( PINGFLAG .AND. NPING .EQ. 0 ) THEN
                MESG = 'WARNING: No sources are PinG sources in ' //
     &                 'input file, so PinG ' // CRLF() // BLANK10 //
     &                 'emissions file will not be written.'
                CALL M3MSG2( MESG )
                PINGFLAG = .FALSE.
            END IF

C.............  Read stack group IDs
            IF ( .NOT. READ3( PVNAME, 'ISTACK', 1,
     &                        PVSDATE, PVSTIME, GRPGID ) ) THEN

                L2 = LEN_TRIM( PVNAME )
                MESG = 'Could not read "ISTACK" from file "' //
     &                 PVNAME( 1:L2 ) // '"'
                CALL M3EXIT( PROGNAME, SDATE, 0, MESG, 2 )

            END IF

C.............  Update elevated sources filter for elevated sources
            DO S = 1, NPSRC
                IF( GROUPID( S ) .GT. 0 ) ELEVFLTR( S ) = 1.
            END DO

        END IF

C.........  Create complete source list for explicit elevated sources
        IF( EXPLFLAG ) THEN

C.............  Allocate memory for temporary unsorted list and index
            ALLOCATE( SRCFLG( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCFLG', PROGNAME )
            ALLOCATE( TMPSRC( NHRSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPSRC', PROGNAME )
            ALLOCATE( TMPIDX( NHRSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPIDX', PROGNAME )
            SRCFLG = .FALSE.

C.............  Loop through hours in PLAY_EX file and determine all sources
C               that are listed for all hours
            JDATE  = SDATE
            JTIME  = STIME
            N      = 0
            DO T = 1, NSTEPS
                IF( .NOT. READ3( PHNAME, 'INDXH', ALLAYS3,
     &                           JDATE, JTIME, INDXH      ) ) THEN

                    MESG = 'Could not read INDXH from ' // PHNAME
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF   ! if read3() failed

C.................  Loop through this hour's indices, and store any new ones
                DO I = 1, NHRSRC

                    S = INDXH( I )

C.....................  Exit loop if done sources for this hour
                    IF ( S .EQ. 0 ) EXIT

C.....................  Store if not already
                    IF ( .NOT. SRCFLG( S ) ) THEN
                        SRCFLG( S ) = .TRUE.
                        N = N + 1
                        TMPSRC( N ) = S
                        TMPIDX( N ) = N
                    END IF

                END DO

                CALL NEXTIME( JDATE, JTIME, TSTEP )

            END DO
            NHRSRC = N   ! Reset to permit FINDs in case whole PLAY_EX not used

C.............  Sort list index
            CALL SORTI1( NHRSRC, TMPIDX, TMPSRC )

C.............  Store final sorted list
            DO I = 1, NHRSRC
                J = TMPIDX( I )
                ELEVSRC( I ) = TMPSRC( J )
            END DO

C.............  Deallocate temporary memory
            DEALLOCATE( SRCFLG, TMPSRC, TMPIDX )

        END IF

C.........  Read reactivity matrices
        IF( ARFLAG ) CALL RDRMAT( ARNAME, ANSREAC, ARNMSPC, ACRIDX,
     &                            ACRREPEM, ACRPRJFC, ACRMKTPN, ACRFAC )

        IF( MRFLAG ) CALL RDRMAT( MRNAME, MNSREAC, MRNMSPC, MCRIDX,
     &                            MCRREPEM, MCRPRJFC, MCRMKTPN, MCRFAC )

        IF( PRFLAG ) CALL RDRMAT( PRNAME, PNSREAC, PRNMSPC, PCRIDX,
     &                            PCRREPEM, PCRPRJFC, PCRMKTPN, PCRFAC )

C.........  Read gridding matrices (note, must do through subroutine because of
C           needing contiguous allocation for integer and reals)
        IF( AFLAG ) CALL RDGMAT( AGNAME, NGRID, ANGMAT, ANGMAT,
     &                           AGMATX(1), AGMATX( NGRID+1 ),
     &                           AGMATX( NGRID+ANGMAT+1 ) )

        IF( MFLAG ) CALL RDGMAT( MGNAME, NGRID, MNGMAT, MNGMAT,
     &                           MGMATX(1), MGMATX( NGRID+1 ),
     &                           MGMATX( NGRID+MNGMAT+1 ) )

        IF( PFLAG ) THEN

            PGMATX = 1.  ! initialize array b/c latter part not in file
            CALL RDGMAT( PGNAME, NGRID, NPSRC, 1,
     &                   PGMATX(1), PGMATX( NGRID + 1 ), RDUM )
        END IF

C.........  Build indicies for pollutant/species groups

        CALL BLDMRGIDX( MXGRP, MXVARPGP, NGRP )

C.........  Intialize state/county summed emissions to zero

        CALL INITSTCY

C.........  Read source group cross-reference file and assign sources to groups

        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
            CALL RDSRCGRPS( SGDEV, .FALSE., .FALSE. )
        END IF

C.........  Open NetCDF output files, open ASCII report files, and write headers

        CALL OPENMRGOUT( NGRP )

C.........  In case reactivity does not exist, initialize temporary arrays
C           for reactivity information anyway.  These are used even without
C           reactivity matrix inputs so that the code does not need even
C           more conditionals in the matrix multiplication step.
        IF( AFLAG ) ARINFO = 0.  ! array
        IF( MFLAG ) MRINFO = 0.  ! array
        IF( PFLAG ) PRINFO = 0.  ! array

C.........  Allocate memory for temporary list of species and pollutant names
        ALLOCATE( VARNAMES( MXVARPGP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VARNAMES', PROGNAME )
        ALLOCATE( INNAMES( MXVARPGP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INNAMES', PROGNAME )

C.........   Allocate memory for logical array for whether species set for elevated output
        ALLOCATE( ELEV_SPCSET( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELEV_SPCSET', PROGNAME )
        ELEV_SPCSET = .FALSE.  ! array

C..........  Allocate memory for tmp control/projection array
        IF( PFACFLAG ) THEN
            IF( AUFLAG ) NTSRC = NASRC
            IF( MUFLAG ) NTSRC = NMSRC
            IF( PUFLAG ) NTSRC = NPSRC
            ALLOCATE( TCUMATX( NTSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TCUMATX', PROGNAME )
            TCUMATX = 0.0
        END IF

C.........  Loop through processing groups (if speciation, this will be specia-
C           tion groups, but if no speciation, this will be pollutant groups,
C           for purposes of memory usage if many pollutants and/or species)
        PGID = IMISS3
        DO N = 1, NGRP

C.............  Set the number of variables per group
            NVPGP = VGRPCNT( N )

C.............  If pollutants in current group are different from those
C               in the previous group, read pollutant-specific control matrices
C.............  For reactivity matrices, read inventory emissions that will
C               be needed for getting ratios of inventory to hourly for applying
C               reactivity-based projection to hourly emissions
C.............  Note that only the pollutants in this group that are actually
C               in the control matrices are stored, and the index that says
C               which are valid is *U_EXIST and *A_EXIST
            IF( IDVGP( N ) .NE. PGID ) THEN

                IF( AUFLAG ) THEN
                    IF( PFACFLAG ) THEN
                        CALL RD3MASK( AUNAME,0,0, NASRC, 1, 1, 'pfac',
     &                                1, TCUMATX )
                        DO NG = 1, NVPGP
                            J = AU_EXIST( NG,1 )
                            ACUMATX( :,J ) = TCUMATX( : )
                        END DO
                    ELSE
                        CALL RD3MASK( AUNAME,0,0, NASRC, AMULSIZ, NVPGP,
     &                      GVNAMES( 1,N ), AU_EXIST( 1,N ), ACUMATX )
                    ENDIF
                ENDIF

                IF( MUFLAG ) THEN
                    IF( PFACFLAG ) THEN
                        CALL RD3MASK( MUNAME,0,0, NMSRC, 1, 1,'pfac',
     &                                1, TCUMATX )
                        DO NG = 1, NVPGP
                            J = MU_EXIST( NG,1 )
                            MCUMATX( :,J ) = TCUMATX( : )
                        END DO
                    ELSE
                        CALL RD3MASK( MUNAME,0,0, NMSRC, MMULSIZ, NVPGP,
     &                      GVNAMES( 1,N ), MU_EXIST( 1,N ), MCUMATX )
                    ENDIF
                ENDIF

                IF( PUFLAG )  THEN
                    IF( PFACFLAG ) THEN
                        CALL RD3MASK( PUNAME,0,0, NPSRC, 1, 1, 'pfac',
     &                                1 , TCUMATX )
                        DO NG = 1, NVPGP
                            J = PU_EXIST( NG,1 )
                            PCUMATX( :,J ) = TCUMATX( : )
                        END DO
                    ELSE
                        CALL RD3MASK( PUNAME,0,0, NPSRC, PMULSIZ, NVPGP,
     &                      GVNAMES( 1,N ), PU_EXIST( 1,N ), PCUMATX )
                    ENDIF
                ENDIF

                IF( ARFLAG )
     &              CALL RD3MASK( AENAME, 0, 0, NASRC, APOLSIZ, NVPGP,
     &                      GVNAMES( 1,N ), A_EXIST( 1,N ), AEISRC   )

                IF( MRFLAG )
     &              CALL RD3MASK( MENAME, 0, 0, NMSRC, MPOLSIZ, NVPGP,
     &                      GVNAMES( 1,N ), M_EXIST( 1,N ), MEISRC   )

                IF( PRFLAG )
     &              CALL RD3MASK( PENAME, 0, 0, NPSRC, PPOLSIZ, NVPGP,
     &                      GVNAMES( 1,N ), P_EXIST( 1,N ), PEISRC   )

            END IF

C.............  Loop through variables in current group...
            OCNT = 0
            LBUF = ' '
            INNAMES  = ' '  ! array
            VARNAMES = ' '  ! array
            DO V = 1, NVPGP  ! No. variables per group

                K1 = 0
                K2 = 0
                K3 = 0

C.................  Extract name of variable in group
                VBUF = GVNAMES( V,N )

C.................  For speciation...
                IF( SFLAG ) THEN

C.....................  Update list of output species names for message
                    SBUF = EMNAM( SPINDEX( V,N ) )
                    PBUF = EANAM( SIINDEX( V,N ) )
                    M = INDEX1( SBUF, OCNT, VARNAMES )

                    IF( M .LE. 0 .AND. SBUF .NE. LBUF ) THEN
                        OCNT = OCNT + 1
                        VARNAMES( OCNT ) = SBUF
                        LBUF = SBUF
                    END IF

C.....................  Set position for input of speciation matrix
                    IF( AFLAG ) K1 = AS_EXIST( V,N )
                    IF( MFLAG ) K2 = MS_EXIST( V,N )
                    IF( PFLAG ) K3 = PS_EXIST( V,N )

C.....................  Read speciation matrix for current variable and
C                       position
                    IF ( K1 .GT. 0 )
     &                    CALL RDSMAT( ASNAME, VBUF, ASMATX( 1,K1 ) )
                    IF ( K2 .GT. 0 )
     &                    CALL RDSMAT( MSNAME, VBUF, MSMATX( 1,K2 ) )
                    IF ( K3 .GT. 0 )
     &                    CALL RDSMAT( PSNAME, VBUF, PSMATX( 1,K3 ) )

C.................  For no speciation, prepare list of variables for output mesg
                ELSE

C.....................  Update list of pollutants names for message
                    PBUF = EANAM( SIINDEX( V,N ) )
                    M = INDEX1( PBUF, OCNT, VARNAMES )

                    IF( M .LE. 0 ) THEN
                        OCNT = OCNT + 1
                        VARNAMES( OCNT ) = PBUF
                    END IF

                END IF  ! end speciation or not

C.................  Set input variable names
                INNAMES ( V ) = TONAMES( SIINDEX( V,N ) )

            END DO      ! End variables in group loop

C.............  Write out message about data currently being processed
            CALL POLMESG( OCNT, VARNAMES )

C.............  Initializations before main time loop
            JDATE  = SDATE
            JTIME  = STIME
            LDATE  = 0
            DAY    = 1

C.............  Loop through output time steps
            DO T = 1, NSTEPS   ! at least once for time-independent

C.................  Reinitialize array for flagging elevated species as "set"
                ELEV_SPCSET = .FALSE.  ! array

C................. For time-dependent processing, write out a few messages...
                IF( TFLAG ) THEN

C.....................  Determine weekday index (Monday is 1)
                    DAY = WKDAY( JDATE )

C.....................  Write out message for new day.  Note, For time-
C                       independent, LDATE and JDATE will both be zero.
                    IF( JDATE .NE. LDATE ) THEN

                        CALL WRDAYMSG( JDATE, MESG )

                    END IF

C.....................  Write out files that are being used for by-day treatment
                    IF( DAY .NE. PDAY ) THEN

                        IF( AFLAG_BD ) THEN
                            MESG = '   with ATMP file ' // ATNAME( DAY )
                            CALL M3MSG2( MESG )
                        END IF

                        IF( MFLAG_BD ) THEN
                            MESG = '   with MTMP file ' // MTNAME( DAY )
                            CALL M3MSG2( MESG )
                        END IF

                        IF( PFLAG_BD ) THEN
                            MESG = '   with PTMP file ' // PTNAME( DAY )
                            CALL M3MSG2( MESG )
                        END IF

                        PDAY = DAY

                    END IF

C.....................  For new hour...
C.....................  Write to screen because WRITE3 only writes to LDEV
                    WRITE( *, 93020 ) HHMMSS( JTIME )

                END IF

C.................  Initialize source-category current dates
                AJDATE = JDATE
                MJDATE = JDATE
                PJDATE = JDATE

C.................  Reset the date for each source category when by-day
C                   processing is being done for that category
                IF( AFLAG_BD ) AJDATE = ASDATE( DAY )
                IF( MFLAG_BD ) MJDATE = MSDATE( DAY )
                IF( PFLAG_BD ) PJDATE = PSDATE( DAY )

C.................  If area sources, read inventory emissions for this time
C                   step for all area-source pollutants in current pol group
C.................  The *_EXIST are counters that point to the position in
C                   the source category emissions of the variables names
C                   in INNAMES. Data are stored in *EMSRC in the global order.
                IF( AFLAG ) THEN

C.................  If using map-formatted inventory for time-independent
                    IF( ANMAP .NE. 0 .AND. .NOT. TFLAG ) THEN

                        CALL RDMAPMASK( AENAME, ANMAP, AMAPNAM, AMAPFIL,
     &                               NASRC, APOLSIZ, NVPGP, VARNAMES(1),
     &                               INNAMES(1), A_EXIST(1,N), AEMSRC  )

C.................  If using hourly data
                    ELSE
                        CALL RDSETMASK( ATNAME( DAY ), AJDATE, JTIME,
     &                                NASRC, APOLSIZ, NVPGP, INNAMES(1),
     &                                A_EXIST( 1,N ), AEMSRC )
                    END IF
                END IF

C.................  If mobile sources, read inventory emissions or activities
C                   for this time step for all mobile-source pollutants in
C                   current pol group
                IF( MFLAG ) THEN

C.................  If using map-formatted inventory for time-independent
                    IF( MNMAP .NE. 0 .AND. .NOT. TFLAG ) THEN
                        CALL RDMAPMASK( MENAME, MNMAP, MMAPNAM, MMAPFIL,
     &                               NMSRC, MPOLSIZ, NVPGP, VARNAMES(1),
     &                                INNAMES(1), M_EXIST(1,N), MEMSRC )

C.................  If using hourly data
                    ELSE
                        CALL RDSETMASK( MTNAME( DAY ), MJDATE, JTIME,
     &                                NMSRC, MPOLSIZ, NVPGP, INNAMES(1),
     &                                M_EXIST( 1,N ), MEMSRC )
                    END IF
                END IF

C.................  If point sources, read inventory emissions for this time
C                   step for all point-source pollutants in current pol group
                IF( PFLAG ) THEN

C.................  If using map-formatted inventory for time-independent
                    IF( PNMAP .NE. 0 .AND. .NOT. TFLAG ) THEN
                        CALL RDMAPMASK( PENAME, PNMAP, PMAPNAM, PMAPFIL,
     &                               NPSRC, PPOLSIZ, NVPGP, VARNAMES(1),
     &                                INNAMES(1), P_EXIST(1,N), PEMSRC )

C.................  If using hourly data
                    ELSE
                        CALL RDSETMASK( PTNAME( DAY ), PJDATE, JTIME,
     &                                NPSRC, PPOLSIZ, NVPGP, INNAMES(1),
     &                                P_EXIST( 1,N ), PEMSRC )
                    END IF
                END IF

C.................  If layer fractions, read them for this time step
                IF( LFLAG ) THEN

                    IF( .NOT. READ3( PLNAME, 'LFRAC', ALLAYS3,
     &                               JDATE, JTIME, LFRAC      ) ) THEN

                        MESG = 'Could not read LFRAC from ' // PLNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF   ! if read3() failed

C.................  Otherwise, if explicit plume rise, read fractions and
C                   indices from the file
                ELSE IF ( EXPLFLAG ) THEN

                    IF( .NOT. READ3( PHNAME, 'INDXH', ALLAYS3,
     &                               JDATE, JTIME, INDXH      ) ) THEN

                        MESG = 'Could not read INDXH from ' // PHNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF   ! if read3() failed

                    IF( .NOT. READ3( PHNAME, 'LFRAC', ALLAYS3,
     &                               JDATE, JTIME, LFRAC      ) ) THEN

                        MESG = 'Could not read LFRAC from ' // PHNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF   ! if read3() failed

                END IF

C.................  Loop through variables in the current group
                LBUF = ' '
                DO V = 1, NVPGP

C.....................  Set species or pollutant/activity name for this
C                       iteration
                    IF( SFLAG ) THEN
                        SBUF = EMNAM( SPINDEX( V,N ) )
                        IF( AFLAG ) KA  = INDEX1( SBUF, ANMSPC, AEMNAM )
                        IF( BFLAG ) KB  = INDEX1( SBUF, BNMSPC, BEMNAM )
                        IF( MFLAG ) KM  = INDEX1( SBUF, MNMSPC, MEMNAM )
                        IF( PFLAG ) KP  = INDEX1( SBUF, PNMSPC, PEMNAM )
                    ELSE
                        SBUF = EANAM( SIINDEX( V,N ) )
                    END IF

C.....................  Set conversion factors
                    IF( SFLAG ) THEN
                        F1 = GRDFAC( SPINDEX( V,N ) )
                        F2 = TOTFAC( SPINDEX( V,N ) )
                    ELSE
                        F1 = GRDFAC( SIINDEX( V,N ) )
                        F2 = TOTFAC( SIINDEX( V,N ) )
                    END IF

C.....................  If area reactivity matrix applies, pre-compute
C                       source array of reactivity emissions & mkt pentrtn
                    IF( ARFLAG ) THEN
                        K1 = A_EXIST ( V,N )
                        K2 = AR_EXIST( V,N )
                        IF( K2 .GT. 0 ) THEN
                            CALL APPLREAC( NASRC, ANSREAC, K1, K2,
     &                             APRJFLAG, LMKTPON, AEISRC,AEMSRC,
     &                             ACRIDX, ACRREPEM, ACRPRJFC,
     &                             ACRMKTPN, ACRFAC, ARINFO )

                        ELSE
                            ARINFO = 0.  ! array
                        END IF
                    END IF

C.....................  Process for area sources...
                    IF( AFLAG ) THEN

                        K1 = A_EXIST ( V,N )
                        K2 = AU_EXIST( V,N )
                        K4 = AS_EXIST( V,N )
                        K5 = NGRID + ANGMAT + 1

C.............................  Apply valid matrices & store
                        CALL MRGMULT( NASRC, NGRID, 1, ANGMAT,
     &                         ANGMAT, K1, K2, K4, KA, F1, F2,
     &                         AEMSRC, ARINFO, ACUMATX, ASMATX,
     &                         AGMATX(1), AGMATX(NGRID+1),
     &                         AGMATX(K5), AICNY, AEMGRD, TEMGRD,
     &                         AEBCNY, AEUCNY, AERCNY, AECCNY )
                    END IF

C.....................  For biogenic sources, read gridded emissions,
C                       add to totals and store
                    IF( BFLAG ) THEN

                        K4 = BS_EXIST( V,N )

                        IF( K4 .GT. 0 ) THEN
                            CALL MRGBIO( SBUF, BTNAME, JDATE, JTIME,
     &                                   NGRID, BIOGFAC, BEMGRD,
     &                                   TEMGRD( 1,1 ) )


C.............................  Update country, state, & county totals
C.............................  Also convert the units from the gridded output
C                               units to the totals output units
                            IF( LREPANY .OR. SRCGRPFLAG ) THEN
                                FB = BIOTFAC / BIOGFAC
                                CALL GRD2CNTY( 0, KB, NCOUNTY,
     &                                         FB, BEMGRD, BEBCNY,
     &                                         SRCGRPFLAG, BIOGFAC,
     &                                         ISRCGRP, EMGGRD )

                            END IF
                        END IF

                    END IF

C.....................  If mobile reactivity matrix applies, pre-compute
C                       source array of reacvty emissions and mkt pntrtn
                    IF( MRFLAG ) THEN
                        K1 = M_EXIST ( V,N )
                        K2 = MR_EXIST( V,N )
                        IF( K2 .GT. 0 ) THEN
                            CALL APPLREAC( NMSRC, MNSREAC, K1, K2,
     &                             MPRJFLAG, LMKTPON, MEISRC,MEMSRC,
     &                             MCRIDX, MCRREPEM, MCRPRJFC,
     &                             MCRMKTPN, MCRFAC, MRINFO )

                        ELSE
                            MRINFO = 0.  ! array
                        END IF

                    END IF

C.....................  Process for mobile sources...
                    IF( MFLAG ) THEN

                        K1 = M_EXIST ( V,N )
                        K2 = MU_EXIST( V,N )
                        K4 = MS_EXIST( V,N )
                        K5 = NGRID + MNGMAT + 1

C.........................  Apply valid matrices & store

                        CALL MRGMULT( NMSRC, NGRID, 1, MNGMAT,
     &                         MNGMAT, K1, K2, K4, KM, F1, F2,
     &                         MEMSRC, MRINFO, MCUMATX, MSMATX,
     &                         MGMATX(1), MGMATX(NGRID+1),
     &                         MGMATX(K5), MICNY, MEMGRD, TEMGRD,
     &                         MEBCNY, MEUCNY, MERCNY, MECCNY )

                    END IF

C.....................  If reactivity matrix applies, pre-compute source
C                       array of reactivity emissions and market penetration
                    IF( PRFLAG ) THEN
                        K1 = P_EXIST ( V,N )
                        K2 = PR_EXIST( V,N )
                        IF( K2 .GT. 0 ) THEN
                            CALL APPLREAC( NPSRC, PNSREAC, K1, K2,
     &                             PPRJFLAG, LMKTPON, PEISRC,PEMSRC,
     &                             PCRIDX, PCRREPEM, PCRPRJFC,
     &                             PCRMKTPN, PCRFAC, PRINFO )
                        ELSE
                            PRINFO = 0.  ! array
                        END IF
                    END IF

C.....................  Process for point sources...
                    IF( PFLAG ) THEN

                        K1 = P_EXIST ( V,N )
                        K2 = PU_EXIST( V,N )
                        K4 = PS_EXIST( V,N )
                        K5 = NGRID + NPSRC + 1

C.........................  Apply valid matrices & store
                        CALL MRGMULT( NPSRC, NGRID, EMLAYS, NPSRC,
     &                         NPSRC, K1, K2, K4, KP, F1, F2,
     &                         PEMSRC, PRINFO, PCUMATX, PSMATX,
     &                         PGMATX(1), PGMATX(NGRID+1),
     &                         PGMATX(K5), PICNY, PEMGRD, TEMGRD,
     &                         PEBCNY, PEUCNY, PERCNY, PECCNY )

C.........................  Apply matrices for elevated and plume-in-grid
C                           outputs, if this pollutant is used for point srcs.
                        IF( K1. GT. 0 .AND.
     &                   (ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG) ) THEN

C.............................  Determine whether or not this species has been
C                               merged before to set initialization flag. This
C                               additional step is to handle the case where multiple
C                               pollutants feed the same species.
                            INITELEV = .TRUE.
                            IF ( SFLAG ) THEN
                                IF( ELEV_SPCSET( SPINDEX( V,N ) ) )
     &                              INITELEV = .FALSE.
                                ELEV_SPCSET( SPINDEX( V,N ) ) = .TRUE.
                            END IF

                            CALL MRGELEV( NPSRC, NMAJOR, NPING,
     &                                    K1, K2, K4, F1, INITELEV )

                        END IF

                    END IF

C.....................  Check the flag that indicates the entries for which
C                       we need to output the gridded data
                    IF( GVLOUT( V,N ) ) THEN

C.........................  Write out gridded data and Models-3 PinG file
                        CALL WMRGEMIS( SBUF, JDATE, JTIME )

C.........................  Write out ASCII elevated sources file
                        IF( ELEVFLAG ) THEN
                            CALL WMRGELEV( SBUF, NPSRC, NMAJOR,
     &                                     JDATE, JTIME        )
                        END IF

                        IF( SRCGRPFLAG ) THEN
                            CALL WRSRCGRPS( SBUF, JDATE,JTIME,.FALSE.,0)

                        ELSE IF( SUBSECFLAG ) THEN
                            DO K = 1, NGRPS
                                IF( .NOT. WRITESET( SUBOUTNAME( K ), SBUF, ALLFILES,
     &                                               JDATE, JTIME, EMGGRD(1,K) ) ) THEN
                                    MESG = 'Could not write "' // SBUF // '" ' //
     &                                      'to file "' // SUBOUTNAME(K) // '"'
                                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                                END IF
                            END DO

                        END IF

C.........................  Initialize gridded arrays
                        IF( AFLAG ) THEN
                            AEMGRD = 0.  ! array
                        ENDIF

                        IF( MFLAG ) THEN
                            MEMGRD = 0.  ! array
                        ENDIF

                        IF( PFLAG ) THEN
                            PEMGRD = 0.  ! array
                        ENDIF

                        TEMGRD = 0.      ! array

                        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                            EMGGRD = 0.  ! array
                        END IF
                    END IF

                END DO      ! End loop on variables in group

C.................  Write country, state, and county emissions (all that apply)
C.................  The subroutine will only write for certain hours and
C                   will reinitialize the totals after output
                IF( LREPANY ) THEN
                    CALL WRMRGREP( JDATE, JTIME, N )
                END IF

                LDATE = JDATE

                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

            END DO          ! End loop on time steps

        END DO   ! End of loop on pollutant/pol-to-spcs groups

C.........  Successful completion of program

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A )

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

94020   FORMAT( A, I4, 2X, 10 ( A, :, 1PG14.6, :, 2X ) )

94030   FORMAT( 8X, 'at time ', A8 )


        END PROGRAM SMKMERGE
