
        SUBROUTINE GENREACT( NSRC, NIPOL, BYEAR, PYEAR, ENAME,
     &                       RPOL, EINAM )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine processes the reactivity data ...
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C************************************************************************
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the speciation profiles
        USE MODSPRO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDC
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE

        EXTERNAL   CRLF, FINDC, INDEX1, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC   ! no. sources
        INTEGER     , INTENT (IN) :: NIPOL  ! no. inventory pollutants
        INTEGER     , INTENT (IN) :: BYEAR  ! base year of react. projections
        INTEGER     , INTENT (IN) :: PYEAR  ! projection year for reactivity
        CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name
        CHARACTER(LEN=IOVLEN3), INTENT (IN) :: RPOL! pol for procssing reactvity
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL )  ! all pollutant names

C...........   Local arrays allocated by subroutine arguments
        INTEGER          ISREA( NSRC )   ! reactivity control data table index
        REAL             EMIS ( NSRC )   ! base inventory emissions

C.........  Names for output variables in mass-based and mole-based files
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: MASSONAM( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: MOLEONAM( : )

C...........   Logical names and unit numbers
        INTEGER, SAVE :: RDEV         !  speciation profiles unit no.
        INTEGER       :: SDEV         !  supplement file unit no.

        CHARACTER*16     SNAME   ! logical name for mass-based react. cntrl mat
        CHARACTER*16     LNAME   ! logical name for mole-based react. cntrl mat

C...........   Other local variables

        INTEGER          I, J, K, N, P, S, V     ! counters and indices

        INTEGER          IDUM        ! dummy integer
        INTEGER          IOS         ! i/o error status
        INTEGER          SIC         ! tmp SIC code for report
        INTEGER          ITBL        ! position in full table of current profile
        INTEGER          LP          ! length of RPOL string
        INTEGER          NMSPC       ! number of model species
        INTEGER          NSREAC      ! number of srcs w/ reactivity controls
        INTEGER          NTBL        ! number of species in current profile

        REAL             EMREP       ! tmp replacement emissions

        CHARACTER(LEN=SPNLEN3) SPCODE ! tmp speciation profile code

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
        LOGICAL, SAVE :: MASSOUT  = .FALSE.  ! true: output mass-based spc facs
        LOGICAL, SAVE :: MOLEOUT  = .FALSE.  ! true: output mole-based spc facs
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: point source processing

        CHARACTER*4      OUTTYPE             ! speciation output type
        CHARACTER*300    MESG                ! message buffer

        CHARACTER*16  :: PROGNAME = 'GENREACT' ! program name

C***********************************************************************
C   begin body of subroutine GENREACT

        IF( FIRSTIME ) THEN 

C.............  Get environment variables that control subroutine behavior
C.............  Retrieve the type of speciation outputs (mass,mole,or all)
            MESG = 'Type of speciation outputs to create'  
            CALL ENVSTR( 'SPEC_OUTPUT', MESG, 'ALL', OUTTYPE, IOS )

C.............  Set flags that depend on the value of OUTTYPE
            CALL UPCASE( OUTTYPE )
            MASSOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
            MOLEOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
            MASSOUT = ( INDEX( OUTTYPE, 'MASS' ) .GT. 0 .OR. MASSOUT )
            MOLEOUT = ( INDEX( OUTTYPE, 'MOLE' ) .GT. 0 .OR. MOLEOUT )

            MESG = 'Speciation profiles needed to process ' //
     &             'reactivity packet...'
            CALL M3MSG2( MESG )

            RDEV = PROMPTFFILE( 
     &             'Enter logical name for SPECIATION PROFILES file',
     &             .TRUE., .TRUE., 'GSPRO', PROGNAME )

            CALL DSCSPROF( RDEV, NIPOL, EINAM )

C.............  Allocate memory for arrays of speciation tables and unique lists
C               using the maximum number of profile table entires per pollutant, 
C               MXSPFUL, which is from module MODSPRO
            ALLOCATE( INPRF( MXSPFUL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INPRF', PROGNAME )
            ALLOCATE( SPECID( MXSPFUL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPECID', PROGNAME )
            ALLOCATE( MOLEFACT( MXSPFUL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MOLEFACT', PROGNAME )
            ALLOCATE( MASSFACT( MXSPFUL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MASSFACT', PROGNAME )

C.............  Determine if source category is point sources, or other
            PFLAG = ALLOCATED( CSOURC )

            FIRSTIME = .FALSE.

        END IF

C.........  For pollutant in subroutine argument...

C.........  Determine which reactivity packet goes to each source
        CALL ASGNCNTL( NSRC, 1, 'REACTIVITY', RPOL, IDUM, ISREA )

C.........  Count up the number of sources that have reactivity data
        NSREAC = 0
        DO S = 1, NSRC

            IF( ISREA( S ) .GT. 0 ) NSREAC = NSREAC + 1

        END DO

C.........  Read speciation profiles file

        LP = LEN_TRIM( RPOL )
        MESG = 'Reading SPECIATION PROFILES file for ' // RPOL( 1:LP )
        CALL M3MSG2( MESG )

        CALL RDSPROF( RDEV, RPOL, MXSPFUL, NSPFUL, NMSPC,
     &                INPRF, SPECID, MOLEFACT, MASSFACT )

C.........  Ensure that profile(s) exist for this pollutant
        IF( NSPFUL .LE. 0 ) THEN
            MESG = 'No speciation profiles found for pollutant "' // 
     &             RPOL( 1:LP ) // '"! '
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Abridge profiles so that there is an array of unique profiles
        V = INDEX1( RPOL, NIPOL, EINAM )
        CALL PROCSPRO( NMSPC, SPCNAMES( 1,V ) )

C.........  Allocate memory for the compressed reactivity matrix.  Use data
C           structures for point sources, but this routine can be used for area
C           sources or mobile sources as well. 

        ALLOCATE( PCRIDX( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRIDX', PROGNAME )
        ALLOCATE( PCRREPEM( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRREPEM', PROGNAME )
        ALLOCATE( PCRPRJFC( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRPRJFC', PROGNAME )
        ALLOCATE( PCRMKTPN( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRMKTPN', PROGNAME )
        ALLOCATE( PCRCSCC( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRCSCC', PROGNAME )
        ALLOCATE( PCRSPROF( NSREAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PCRSPROF', PROGNAME )

        IF( MASSOUT ) THEN
            ALLOCATE( RMTXMASS( NSREAC, NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'RMTXMASS', PROGNAME )
        END IF

        IF( MOLEOUT ) THEN
            ALLOCATE( RMTXMOLE( NSREAC, NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'RMTXMOLE', PROGNAME )
        END IF

C.........  Allocate memory for names of output variables
        ALLOCATE( MASSONAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MASSONAM', PROGNAME )
        ALLOCATE( MOLEONAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLEONAM', PROGNAME )

C.........  Read emissions for current pollutant from the inventory file
        IF( .NOT. READ3( ENAME, RPOL, ALLAYS3, 0, 0, EMIS ) ) THEN
            MESG = 'Could not read "' // RPOL( 1:LEN_TRIM( RPOL ) ) //
     &             '" from inventory file ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Loop through all sources and store reactivity information for
C           those that have it
        N = 0
        DO S = 1, NSRC

            K = ISREA( S )       ! index to reactivity data tables

            IF( K .GT. 0 ) THEN 

C.................  For storing inventory emissions, use the value from the
C                   reactivity table only if it is greater than zero. Otherwise,
C                   store the original base-year emissions from the inventory. 
                IF( EMREPREA( K ) .GT. 0 ) THEN
                    EMREP = EMREPREA( K )
                ELSE
                    EMREP = EMIS( S )
                END IF

C.................  Make sure speciation profile is valid for the current source
C                   and pollutant
                SPCODE = CSPFREA( K )
                V = MAX( FINDC( SPCODE, NSPROF, SPROFN ), 0 )

                IF( V .EQ. 0 ) THEN

                    CALL WRITE_WARNING
                    CYCLE  ! To next iteration

C.................  Get indices to full speciation table
                ELSE

                    ITBL = IDXSPRO ( V )
                    NTBL = NSPECIES( V )

                END IF

C.................  Increment counter of reactivity sources
                N = N + 1

C.................  Store various parameters in reactivity matrix output
C                   structures, but double-check for array overflow
                IF( N .LE. NSREAC ) THEN

                    PCRIDX  ( N ) = S
                    PCRREPEM( N ) = EMREP
                    PCRPRJFC( N ) = PRJFCREA( K )
                    PCRMKTPN( N ) = MKTPNREA( K )
                    PCRCSCC ( N ) = CSCCREA ( K )
                    PCRSPROF( N ) = SPCODE

C.....................  Store speciation factors based on speciation profiles
C                       and indices retrieved for SPCODE
                    I = ITBL - 1
                    DO P = 1, NTBL   ! Loop over species for this profile

                        I = I + 1
                        J = IDXSSPEC( V,P )

                        IF( MASSOUT ) THEN
                            RMTXMASS( N,J )= MASSFACT( I )
                        END IF

                        IF( MOLEOUT ) THEN
                            RMTXMOLE( N,J )= MOLEFACT( I )
                        END IF

                    END DO

                END IF  ! End check for overflow

C NOTE: Add output of the information to temporary file on which to base control
C       reports

                SIC = IREASIC( K ) ! (use for report)

            END IF      ! End check if reactivity applies to this source

        END DO

C.........  Reset number of reactivity sources in case any were dropped
        NSREAC = N

C.........  Find position of pollutant in list
        V = INDEX1( RPOL, NIPOL, EINAM )

C.........  Set up for and open output reactivity matrices for current pollutant

        CALL OPENRMAT( 'POINT', ENAME, RPOL, MASSOUT, MOLEOUT, 
     &                 BYEAR, PYEAR, NSREAC, NMSPC, SPCNAMES( 1,V ), 
     &                 SDEV, SNAME, LNAME, MASSONAM, MOLEONAM )

C.........  Write reactivity matrices for current pollutant
        IF( MASSOUT ) THEN
            CALL WRRMAT( NSREAC, NMSPC, SDEV, SNAME, PCRIDX, PCRREPEM, 
     &                   PCRPRJFC, PCRMKTPN, RMTXMASS, 
     &                   PCRCSCC, PCRSPROF, MASSONAM )

        END IF

        IF( MOLEOUT ) THEN
            CALL WRRMAT( NSREAC, NMSPC, SDEV, LNAME, PCRIDX, PCRREPEM, 
     &                   PCRPRJFC, PCRMKTPN, RMTXMOLE, 
     &                   PCRCSCC, PCRSPROF, MOLEONAM )

        END IF

C.........  Deallocate memory used to generate reactivity matrices
        DEALLOCATE( INPRF, SPECID, MOLEFACT, MASSFACT )

        DEALLOCATE( PCRIDX, PCRREPEM, PCRPRJFC, PCRMKTPN, PCRCSCC, 
     &              PCRSPROF )

        IF( ALLOCATED( RMTXMASS ) ) DEALLOCATE( RMTXMASS )
        IF( ALLOCATED( RMTXMOLE ) ) DEALLOCATE( RMTXMOLE )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram writes a warning message for 
C               a bad speciation profile for either point sources or
C               other sources.
            SUBROUTINE WRITE_WARNING

C.............  Local variables
            INTEGER                L2

            CHARACTER*300          BUFFER   ! source information buffer
            CHARACTER*300          MESG     ! message buffer
            CHARACTER(LEN=FIPLEN3) CFIP     ! tmp (character) FIPS code
            CHARACTER(LEN=SRCLEN3) CSRC     ! tmp source chars string
            CHARACTER(LEN=SCCLEN3) TSCC     ! tmp 10-digit SCC

C----------------------------------------------------------------------

C.............  Retrieve SCC for this source
            TSCC = CSCC( S )   ! SCC

C.............  Generate point source part of message
            IF( PFLAG ) THEN

                CSRC = CSOURC( S )   ! source characteristics
                CALL FMTCSRC( CSRC, 7, BUFFER, L2 )

                BUFFER = BUFFER( 1:L2 ) // CRLF() // BLANK10 // 
     &                   'SCC: ' // TSCC // ' POL: ' // RPOL

                L2 = LEN_TRIM( BUFFER )

C.............  Generate area or mobile source part of message
            ELSE

                WRITE( CFIP, '(I6.6)' ) IFIP( S )

                BUFFER = 'FIP: ' // CFIP // ' SCC: ' // TSCC // 
     &                   ' POL: ' // RPOL

                L2 = LEN_TRIM( BUFFER )

            END IF

            MESG = 'WARNING: Speciation profile "' // SPCODE // 
     &             '" is not in profiles, but it was assigned'//
     &             CRLF() // BLANK10 // 'to source:' //
     &             CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &             CRLF() // BLANK5 // 
     &             'Source will be excluded from reactivity matrix.'
            CALL M3MESG( MESG )

            RETURN

            END SUBROUTINE WRITE_WARNING

        END SUBROUTINE GENREACT
