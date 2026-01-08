
        SUBROUTINE GENREACT( PYEAR, ENAME, RPOL, USEPOL )

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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC, CSCC, CIFIP, CISIC, CLINK

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: RMTXMASS, RMTXMOLE, PCRIDX, PCRREPEM,
     &                      PCRPRJFC, PCRMKTPN, PCRCSCC, PCRSPROF,
     &                      RPTDEV, EMREPREA, CSPFREA, PRJFCREA,
     &                      MKTPNREA, CSCCREA

C.........  This module contains the speciation profiles
        USE MODSPRO, ONLY: NSPFUL, NSPROF, SPROFN, SPECID, MOLEFACT,
     &                     MASSFACT, INPRF, MXSPFUL, SPCNAMES, 
     &                     IDXSPRO, IDXSSPEC, NSPECIES

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPOL, MXCHRS, EINAM, NSRC, CATDESC, BYEAR,
     &                     NCHARS, SC_BEGP, SC_ENDP, JSCC, CATEGORY, CRL

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       LOGICAL         ENVYN
C       INTEGER         FINDC
C       INTEGER         INDEX1
C       INTEGER         PROMPTFFILE
C       REAL            YR2DAY

C        EXTERNAL   CRLF, ENVYN, FINDC, INDEX1, PROMPTFFILE, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: PYEAR  ! projection year for reactivity
        CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name
        CHARACTER(IOVLEN3), INTENT (IN) :: RPOL! pol for procssing reactvity
        LOGICAL     , INTENT (IN) :: USEPOL( NIPOL ) ! true: pol valid for pkt

C...........   Local parameters
        INTEGER, PARAMETER :: NHEADER  = 15
        CHARACTER(15), PARAMETER :: HEADERS( NHEADER ) = 
     &                          ( / 'Source         ',
     &                              'Region         ',
     &                              'Plant          ',
     &                              'Char1          ',
     &                              'Char2          ',
     &                              'Char3          ',
     &                              'Char4          ',
     &                              'Link           ',
     &                              'Base SCC       ',
     &                              'Base Emis      ',
     &                              'New Base Emis  ',
     &                              'Proj Factor    ',
     &                              'New SCC        ',
     &                              'New Spc Profile',
     &                              'Mkt Pen Rate   '  / )

C...........   Local arrays allocated by subroutine arguments
        INTEGER          MCWID( MXCHRS )      !  max characteristic width
        INTEGER          HINDX( NHEADER )     !  header label index
        INTEGER          HCWID( NHEADER )     !  header label widths
        LOGICAL          LF   ( MXCHRS )      !  true: column should be output
        CHARACTER(20)    CHARS( MXCHRS )      !  source fields for output

C.........  Allocatable arrays
        INTEGER, ALLOCATABLE :: ISREA( : )   ! reactivity control data table index
        REAL   , ALLOCATABLE :: EMIS ( : )   ! base inventory emissions

        CHARACTER(IOVLEN3), ALLOCATABLE :: MASSONAM( : ) ! mass output species names
        CHARACTER(IOVLEN3), ALLOCATABLE :: MOLEONAM( : ) ! mole output species names

C...........   Logical names and unit numbers
        INTEGER, SAVE :: PDEV         !  speciation profiles unit no.
        INTEGER       :: SDEV         !  supplement file unit no.

        CHARACTER(16) :: RNAME = 'IOAPI_DAT' ! logical name for reading pols
        CHARACTER(16)    SNAME   ! logical name for mass-based react. cntrl mat
        CHARACTER(16)    LNAME   ! logical name for mole-based react. cntrl mat

C...........   Other local variables

        INTEGER          I, J, K, L, N, P, S, V     ! counters and indices

        INTEGER          IDUM        ! dummy integer
        INTEGER          IOS         ! i/o error status
        INTEGER          ITBL        ! position in full table of current profile
        INTEGER          NC          ! tmp number of src chars
        INTEGER          NCOLS       ! no. of output columns in report
        INTEGER          NCOUT       ! no. of src chars to put in report
        INTEGER          NMSPC       ! number of model species
        INTEGER          NSREAC      ! number of srcs w/ reactivity controls
        INTEGER          NTBL        ! number of species in current profile
        INTEGER          RDEV          ! Report unit number

        REAL             EMREP       ! tmp replacement emissions
        REAL             YFAC        ! tmp yr conversion factor

        CHARACTER(SPNLEN3) SPCODE ! tmp speciation profile code

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
        LOGICAL       :: LFLAG    = .FALSE.  ! true: link will be included in report
        LOGICAL, SAVE :: LAVEDAY  = .FALSE.  ! true: use average day emissions
        LOGICAL, SAVE :: MASSOUT  = .FALSE.  ! true: output mass-based spc facs
        LOGICAL, SAVE :: MOLEOUT  = .FALSE.  ! true: output mole-based spc facs
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: point source processing

        CHARACTER(4)     OUTTYPE             ! speciation output type
        CHARACTER(IOVLEN3) VNAM              ! tmp average day var name
        CHARACTER(256)   BUFFER              ! string buffer for building output fmt
        CHARACTER(256)   HDRSTR              ! string for part of header line
        CHARACTER(256)   MESG                ! message buffer
        CHARACTER(256)   OUTFMT              ! output format for report

        CHARACTER(16) :: PROGNAME = 'GENREACT' ! program name

C***********************************************************************
C   begin body of subroutine GENREACT

        IF( FIRSTIME ) THEN 

C.............  Get environment variables that control subroutine behavior
C.............  Retrieve the type of speciation outputs (mass,mole,or all)
            MESG = 'Type of speciation outputs to create'  
            CALL ENVSTR( 'SPEC_OUTPUT', MESG, 'ALL', OUTTYPE, IOS )

C.............  Get environment variables that control program behavior
            MESG = 'Use annual or average day emissions'
            LAVEDAY = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )

C.........  Open reports file
            RPTDEV( 3 ) = PROMPTFFILE( 
     &                     'Enter logical name for REACTIVITY REPORT', 
     &                    .FALSE., .TRUE., CRL // 'REACREP', PROGNAME )
            RDEV = RPTDEV( 3 )

C.............  Set flags that depend on the value of OUTTYPE
            CALL UPCASE( OUTTYPE )
            MASSOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
            MOLEOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
            MASSOUT = ( INDEX( OUTTYPE, 'MASS' ) .GT. 0 .OR. MASSOUT )
            MOLEOUT = ( INDEX( OUTTYPE, 'MOLE' ) .GT. 0 .OR. MOLEOUT )

            MESG = 'Speciation profiles needed to process ' //
     &             'reactivity packet...'
            CALL M3MSG2( MESG )

            PDEV = PROMPTFFILE( 
     &             'Enter logical name for SPECIATION PROFILES file',
     &             .TRUE., .TRUE., 'GSPRO', PROGNAME )

            CALL DSCSPROF( PDEV, NIPOL, EINAM )

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

C...........  Determine if source category is point sources, or other
            PFLAG = ASSOCIATED( CISIC )

C...........  Allocate memory for source-based arrays and initialize
            ALLOCATE( ISREA( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ISREA', PROGNAME )
            ALLOCATE( EMIS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
            ISREA = 0 
            EMIS  = 0.

C...........  Write output report header
            IF( RDEV .GT. 0 ) THEN
                WRITE(RDEV,93000) 'Processed as '// CATDESC// ' sources'
                WRITE(RDEV,93390) 'Base inventory year ', BYEAR
                IF ( PYEAR .GT. 0 ) THEN
                    WRITE(RDEV,93390) 'Reactivity packet projected ' //
     &                                'year ', PYEAR
                END IF

                WRITE( RDEV, 93000 ) 
     &                  'Controls applied with /REACTIVITY/ packet '//
     &                  'for pollutant "' // TRIM( RPOL ) // '".'

                IF( LAVEDAY ) THEN
                    WRITE(RDEV,93000)'Average day data basis in report'
                ELSE
                    WRITE(RDEV,93000)'Annual total data basis in report'
                END IF

            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  For pollutant in subroutine argument...

C.........  Determine which reactivity packet goes to each source
        CALL ASGNCNTL( NSRC, 1, 'REACTIVITY', USEPOL, RPOL, 
     &                 IDUM, ISREA )

C NOTE:
C.........  Please note that the indexing for MCWID, HCWID, and HINDX
C           is confusing. The main confusing aspect is that the
C           +1's are added to index past the "Source" column, which
C           is not a part of the MCWID (which is just for the
C           source characteristics).  The +2's are added to index
C           past the "Source" column and to the column after the
C           last NCOUT column (last column of source chars).

C.........  Initialize valid columns
        LF = .FALSE.  ! array
        DO I = 1, NCHARS
            LF( I ) = .TRUE.
        END DO

C.........  Count up the number of sources that have reactivity data
        NSREAC = 0
        N      = 0
        MCWID   = 0  ! array        
        DO S = 1, NSRC

            IF( ISREA( S ) .GT. 0 ) THEN
                NSREAC = NSREAC + 1

C............... Determine maximum width of report columns for source info
                CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP, 
     &                         LF, NC, CHARS )

                DO J = 1, NC
                    MCWID( J ) = MAX( LEN_TRIM( CHARS(J) ), MCWID(J) )
                    IF( MCWID( J ) .GT. 0 ) N = MAX( N, J )
                END DO

            END IF

        END DO

C........  If SCC is not included in source chars, then add it to report
        IF( JSCC .EQ. 0 ) THEN
            N = N + 1
            MCWID( N ) = SCCLEN3
        END IF

        NCOUT = N

C.........  Determine source characteristic header indices
        HINDX = 0    ! array
        SELECT CASE( CATEGORY )

C........  Area source
        CASE( 'AREA' )
            NCOLS = 3
            HINDX( 1 ) = 1
            HINDX( 2 ) = 2
            HINDX( 3 ) = 9

C........  Mobile source
       CASE( 'MOBILE' )
            NCOLS = 3
            HINDX( 1 ) = 1
            HINDX( 2 ) = 2
            HINDX( 3 ) = 9

C...........  Determine if link will be included in reports
            IF( ALLOCATED( CLINK ) ) THEN
                DO S = 1, NSRC
                    IF( CLINK( S ) .NE. ' ' ) THEN
                        LFLAG = .TRUE.
                        EXIT
                    END IF
                END DO
            END IF
            IF( LFLAG ) THEN
                HINDX( 4 ) = 8
                NCOLS = NCOLS + 1
            END IF

C........  Point source
        CASE( 'POINT' )
            NCOLS = 1
            HINDX( 1 ) = 1
            DO J = 1, NCOUT
                HINDX( J+1 ) = J+1
                IF( JSCC .EQ. J ) HINDX( J+1 ) = 9  ! Ensure SCC is labeled as such
                NCOLS = NCOLS + 1
            END DO

C............  If SCC is not included in source chars, then add it to report
            IF( JSCC .EQ. 0 ) HINDX( NCOUT ) = 9

        END SELECT
    
C.........  Reset max characteristics widths based on headers
        HCWID( 1 ) = 6
        DO J = 1, NCOUT
            HCWID( 1+J ) = LEN_TRIM( HEADERS( HINDX(1+J) ))
            MCWID( J ) = MAX( MCWID(J), HCWID( 1+J ) )
            HCWID( 1+J ) = MCWID( J )
        END DO

C.........  Create remaining header widths - must be consistent with format
C           statments. Maximum is 15 because of HEADERS string size.
        NCOLS = NCOLS + 6
        HCWID( NCOUT + 2 ) = MAX( LEN_TRIM( HEADERS( 10 ) ), 10 )
        HCWID( NCOUT + 3 ) = MAX( LEN_TRIM( HEADERS( 11 ) ), 10 )
        HCWID( NCOUT + 4 ) = MAX( LEN_TRIM( HEADERS( 12 ) ), 7 )
        HCWID( NCOUT + 5 ) = MAX( LEN_TRIM( HEADERS( 13 ) ), SCCLEN3 )
        HCWID( NCOUT + 6 ) = MAX( LEN_TRIM( HEADERS( 14 ) ), SPNLEN3 )
        HCWID( NCOUT + 7 ) = MAX( LEN_TRIM( HEADERS( 15 ) ), 10 )

C.........  Set remaining header indices
        I = 9
        DO J = NCOUT+2, NCOLS
            I = I + 1
            HINDX( J ) = I
        END DO

C.........  Write column headers
        HDRSTR = ' ' // HEADERS( HINDX( 1 ) )( 1:HCWID( 1 ) ) // ';'
        DO J = 2, NCOLS-1
            BUFFER = HDRSTR
            HDRSTR = TRIM( BUFFER ) // ' ' // 
     &               HEADERS( HINDX( J ) )( 1:HCWID( J ) ) // ';'
        END DO
        BUFFER = HDRSTR
        HDRSTR = TRIM( BUFFER ) // ' ' // 
     &           HEADERS( HINDX( NCOLS ) )( 1:HCWID( NCOLS ) )

        WRITE( RDEV, 93000 ) TRIM( HDRSTR )

C.........  Write units headers
        HDRSTR = '       ;'
        DO I = 1, NCOUT
            BUFFER = HDRSTR
            HDRSTR = TRIM( BUFFER ) // REPEAT( ' ',MCWID( I )+1 ) // ';'
        END DO 
        BUFFER = HDRSTR
        HDRSTR = TRIM( BUFFER ) // ' [tons/day];'//
     &           '    [tons/day];            ;           ;' //
     &           '                ; [frac/year]'

        WRITE( RDEV, 93000 ) TRIM( HDRSTR )
        WRITE( RDEV, 93000 ) REPEAT( '-', LEN_TRIM( HDRSTR ) )

C.........  Develop format for output report. Use buffer to prevent run-time
C           error with Portland Group compilers.
        OUTFMT = '(I7,"; "'
        DO J = 1, NCOUT
            BUFFER = OUTFMT
            WRITE( OUTFMT, '(A,I2.2,A)' ) TRIM( BUFFER ) // ',A', 
     &             MCWID( J ), ',"; "'
        END DO
        BUFFER = OUTFMT
        WRITE( OUTFMT,'(7(A,:,I2.2))' ) TRIM( BUFFER ) // ' E',
     &         HCWID( NCOUT+2 ), '.4,"; ",E', 
     &         HCWID( NCOUT+3 ), '.4,"; ",F',
     &         HCWID( NCOUT+4 ), '.4,"; ",A',
     &         HCWID( NCOUT+5 ), '  ,"; ",A',
     &         HCWID( NCOUT+6 ), '  ,"; ",F',
     &         HCWID( NCOUT+7 ), '.4)'

C.........  Read speciation profiles file

        MESG = 'Reading SPECIATION PROFILES file for ' // TRIM( RPOL )
        CALL M3MSG2( MESG )

        CALL RDSPROF( PDEV, RPOL, NMSPC )

C.........  Ensure that profile(s) exist for this pollutant
        IF( NSPFUL .LE. 0 ) THEN
            MESG = 'No speciation profiles found for pollutant "' // 
     &             TRIM( RPOL ) // '"! '
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
            RMTXMASS = 0.
        END IF

        IF( MOLEOUT ) THEN
            ALLOCATE( RMTXMOLE( NSREAC, NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'RMTXMOLE', PROGNAME )
            RMTXMOLE = 0.
        END IF

C.........  Allocate memory for names of output variables
        ALLOCATE( MASSONAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MASSONAM', PROGNAME )
        ALLOCATE( MOLEONAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLEONAM', PROGNAME )

C.........  Read emissions for current pollutant from the inventory file
C.........  Note that average day read will not work if RPOL is
C           more than 13 characters - must use BLDENAMS routine to
C           do this correctly.
        IF( LAVEDAY ) THEN
            VNAM = AVEDAYRT // RPOL( 1:MIN( LEN_TRIM( RPOL ), 13 ) )
        ELSE
            VNAM = RPOL
        END IF

        CALL RDMAPPOL( NSRC, 1, 1, VNAM, EMIS )

C.........  Loop through all sources and store reactivity information for
C           those that have it
        YFAC = YR2DAY( BYEAR )  ! for loop below
        N = 0
        DO S = 1, NSRC

            K = ISREA( S )       ! index to reactivity data tables

            IF( K .GT. 0 ) THEN 

C................. Adjust annual emissions values
                IF( .NOT. LAVEDAY ) EMIS( S ) = EMIS( S ) * YFAC

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
                SPCODE = ADJUSTR( CSPFREA( K ) )
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

C.................  Format source characteristic information for report
                CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP, 
     &                         LF, NC, CHARS )
                
C..............  Write report entry: source number, source characteristics,
C                current base-year emissions, replacement base-year emissions,
C                projection factor, future-year SCC, future-year profile number,
C                market penetration of new SCC/speciation.
                WRITE( RDEV, OUTFMT ) S,
     &               ( CHARS( J )( 1:MCWID(J) ), J=1,NCOUT ), EMIS( S ), 
     &                 EMREP, PRJFCREA( K ), CSCCREA ( K ), SPCODE, 
     &                 MKTPNREA( K )
 
            END IF      ! End check if reactivity applies to this source

        END DO

C.........  Reset number of reactivity sources in case any were dropped
        NSREAC = N

        IF( NSREAC .LE. 0 ) THEN
            MESG = 'No sources assigned reactivity controls.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Find position of pollutant in list
        V = INDEX1( RPOL, NIPOL, EINAM )

C.........  Set up for and open output reactivity matrices for current pollutant

        CALL OPENRMAT( ENAME, RPOL, MASSOUT, MOLEOUT, 
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

c        DEALLOCATE( PCRIDX, PCRREPEM, PCRPRJFC, PCRMKTPN, PCRCSCC )
c     &              PCRSPROF )

c        IF( ALLOCATED( RMTXMASS ) ) DEALLOCATE( RMTXMASS )
c        IF( ALLOCATED( RMTXMOLE ) ) DEALLOCATE( RMTXMOLE )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93390   FORMAT( A, I4.4 )


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

            CHARACTER(300)     BUFFER   ! source information buffer
            CHARACTER(300)     MESG     ! message buffer
            CHARACTER(FIPLEN3) CFIP     ! tmp (character) FIPS code
            CHARACTER(SRCLEN3) CSRC     ! tmp source chars string
            CHARACTER(SCCLEN3) TSCC     ! tmp 10-digit SCC

C----------------------------------------------------------------------

C.............  Retrieve SCC for this source
            TSCC = CSCC( S )   ! SCC

C.............  Generate point source part of message
            IF( PFLAG ) THEN

                CSRC = CSOURC( S )   ! source characteristics
                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                BUFFER = BUFFER( 1:L2 ) // CRLF() // BLANK10 // 
     &                   'SCC: ' // TSCC // ' POL: ' // RPOL

                L2 = LEN_TRIM( BUFFER )

C.............  Generate area or mobile source part of message
            ELSE

                CFIP = CIFIP( S )

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
