
        PROGRAM SPCMAT

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from spcpmat.F 1/99 by M. Houyoux
C
C***********************************************************************
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
C************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the speciation profiles
        USE MODSPRO

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         ENVYN   
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE

        EXTERNAL        CRLF, ENVYN, GETFLINE, PROMPTFFILE, PROMPTMFILE

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C...........   LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C.........   Speciation matrices:

        REAL, ALLOCATABLE :: MASSMATX( :,: )    !  mass-speciation coefficients
        REAL, ALLOCATABLE :: MOLEMATX( :,: )    !  mole-speciation coefficients

C.........  Inventory pollutants actually in the inventory
        LOGICAL               , ALLOCATABLE :: SPCOUT( : ) ! true: output spcs
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: SINAM ( : ) ! spro search pols

C.........  Names for output variables in mass-based and mole-based files
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: MASSONAM( :,: )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: MOLEONAM( :,: )

C.........  Unit numbers and logical file names
        INTEGER         KDEV    !  unit no. for optional pol to pol conversion
        INTEGER         LDEV    !  unit number for log file
        INTEGER         RDEV    !  unit number for speciation profiles file
        INTEGER         SDEV    !  unit number  for ASCII inventory file
        INTEGER         XDEV    !  unit no. for cross-reference file

        CHARACTER*16    ANAME   !  logical name for additive control matrix
        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    SNAME   !  logical name for mass spec matrix output file
        CHARACTER*16    LNAME   !  logical name for mole spec matrix output file

C.........   Other local variables
        INTEGER          I, J, L1, L2, L3, L4, V !  counters and indices

        INTEGER          IOS               ! i/o status
        INTEGER          NINVARR           ! number inventory variables to input
        INTEGER          NMSPC             ! number of model species

        LOGICAL       :: EFLAG   = .FALSE. !  error flag
        LOGICAL       :: KFLAG   = .FALSE. !  if pol to pol convert file or not
        LOGICAL       :: MASSOUT = .TRUE.  !  true: output mass-based matrix
        LOGICAL       :: MOLEOUT = .TRUE.  !  true: output mole-based matrix
        LOGICAL       :: DEFREPRT= .TRUE.  !  true: report default spc profiles
        LOGICAL       :: MULTIPRO= .TRUE.  !  true: multiple profs for pollutant

        CHARACTER*4            OUTTYPE   !  output type from the environment
        CHARACTER*300          MESG      !  message buffer 
        CHARACTER(LEN=IOVLEN3) CBUF      !  smat output name temporary buffer 
        CHARACTER(LEN=IOVLEN3) IBUF      !  pollutant name temporary buffer 
        CHARACTER(LEN=IOVLEN3) SBUF      !  species name temporary buffer 

        CHARACTER*16  :: PROGNAME = 'SPCMAT' ! program name

C***********************************************************************
C   begin body of program SPCMAT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control program behavior
C.........  Retrieve the whether to prompt for and use pollutant conversion file
        KFLAG = ENVYN( 'POLLUTANT_CONVERSION', 
     &                 'Use pollutant-to-pollutant conversion file',
     &                 .FALSE., IOS )

        MESG = 'Type of speciation outputs to create'  
        CALL ENVSTR( 'SPEC_OUTPUT', MESG, 'ALL', OUTTYPE, IOS )

C.........  Set flags that depend on the value of OUTTYPE
        CALL UPCASE( OUTTYPE )
        MASSOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
        MOLEOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
        MASSOUT = ( INDEX( OUTTYPE, 'MASS' ) .GT. 0 .OR. MASSOUT )
        MOLEOUT = ( INDEX( OUTTYPE, 'MOLE' ) .GT. 0 .OR. MOLEOUT )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.......   Get file names and units; open input files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file',
     &           .TRUE., .TRUE., 'GSREF', PROGNAME )

        RDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'GSPRO', PROGNAME )

        IF( KFLAG ) 
     &  KDEV = PROMPTFFILE( 
     &           'Enter logical name for POLLUTANT CONVERSION file',
     &           .TRUE., .TRUE., 'GSCNV', PROGNAME )

C.........  Get header description of inventory file 
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:LEN_TRIM( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Store additional pollutant arrays needed for speciation
            ALLOCATE( SINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SINAM', PROGNAME )
            ALLOCATE( SPCOUT( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPCOUT', PROGNAME )

            DO I = 1, NIPOL
                SINAM ( I ) = EINAM( I )     ! Initialization
                SPCOUT( I ) = .TRUE.         ! Initialize assuming output
            END DO
 
        END IF

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'CSCC'
        IVARNAMS( 2 ) = 'CSOURC'

C.........  Set inventory variables to read for specific source categories
        IF( CATEGORY .EQ. 'AREA' ) THEN
            NINVARR = 2

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
            NINVARR = 4
            IVARNAMS( 3 ) = 'IRCLAS'  ! ??????
            IVARNAMS( 4 ) = 'CLINK'   ! ??????

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
            NINVARR = 2
        END IF

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........   Read the speciation cross-reference file

        CALL RDSREF( XDEV )

C.........  Read the pollutant to pollutant conversion file, if any
C.........  Resulting tables are passed via MODSPRO
        IF ( KFLAG ) THEN

            CALL RDSCONV( KDEV, NIPOL, EINAM, SINAM )

        END IF

C.........  Scan speciation profiles file to get all of the pollutant-species
C           combinations that are valid for the pollutants in the inventory.
C.........  The species names are sorted in ABC order for each pollutant, and
C           and the pollutants are in the same order as SINAM.
C.........  Also retrieve the maximum number of species per pollutant and 
C           maximum number of profile entries per pollutant.

        CALL DSCSPROF( RDEV, NIPOL, SINAM )

C.........  Give warning if some pollutants won't be speciated, and keep track
C           of which ones don't get species.
        J = 0
        DO I = 1, NIPOL

            IF( SPCNAMES( 1,I ) .EQ. ' ' ) THEN
                L1   = LEN_TRIM( EINAM( I ) )
                J    = J + 1
                MESG = 'WARNING: No speciation profiles found ' //
     &                 'for pollutant "' // EINAM( I )( 1:L1 ) // '"' //
     &                 CRLF( )// BLANK10// 'Pollutant ignored!'
                CALL M3MSG2( MESG )

                SPCOUT( I ) = .FALSE.
            END IF

        END DO

C.........  Make sure at least one pollutant will be speciated
        IF( J .EQ. NIPOL ) THEN
            MESG = 'No speciation profiles for pollutants in inventory!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for speciation factors by source using the maximum
C           number of species per pollutant, MXSPEC. Also, initialize.
        IF( MASSOUT ) THEN
            ALLOCATE( MASSMATX( NSRC, MXSPEC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MASSMATX', PROGNAME )
        END IF
        IF( MOLEOUT ) THEN
            ALLOCATE( MOLEMATX( NSRC, MXSPEC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MOLEMATX', PROGNAME )
        END IF

C.........  Allocate memory for arrays of speciation tables and unique lists
C           using the maximum number of profile table entires per pollutant, 
C           MXSPFUL, which is from module MODSPRO
        ALLOCATE( INPRF( MXSPFUL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INPRF', PROGNAME )
        ALLOCATE( SPECID( MXSPFUL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPECID', PROGNAME )
        ALLOCATE( MOLEFACT( MXSPFUL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLEFACT', PROGNAME )
        ALLOCATE( MASSFACT( MXSPFUL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MASSFACT', PROGNAME )

C.........  Allocate memory for names of output variables
        ALLOCATE( MASSONAM( MXSPEC, NIPOL  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MASSONAM', PROGNAME )
        ALLOCATE( MOLEONAM( MXSPEC, NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLEONAM', PROGNAME )

C.........  Open speciation matrix file(s).  Depending on MASSOUT and MOLEOUT,
C           the mass-based and/or mole-based files will be set up and opened.

        CALL OPENSMAT( ENAME, MASSOUT, MOLEOUT, MXSPEC, SPCNAMES, 
     &                 SNAME, LNAME, MASSONAM, MOLEONAM )

C.........  Loop through inventory pollutants to create speciation factors for
C           each species used for each pollutant. In some cases, a pollutant
C           will have no species (e.g., CO), so the factors will simply be a 
C           converstion factor from tons to grams or moles.

        DO V = 1, NIPOL

            IBUF = EINAM( V ) 

C.............  Only process for pollutants we know have at least one species,
C               otherwise, go to end of loop
            IF( SPCOUT( V ) ) THEN

C.................  Write message stating the pollutant being processed
                L1 = LEN_TRIM( IBUF )
                MESG = 'Processing pollutant "' // IBUF( 1:L1 ) // '"'

                IF( EINAM( V ) .NE. SINAM( V ) ) THEN
                    L1 = LEN_TRIM( MESG )
                    L2 = LEN_TRIM( SINAM( V ) )
                    MESG = MESG( 1:L1 ) // ' using pollutant "' //
     &                     SINAM( V )( 1:L2 ) // '" for profiles'
                ENDIF
                CALL M3MSG2( MESG )

            ELSE

C.................  Write message stating the pollutant being skipped
                L1 = LEN_TRIM( IBUF )
                MESG = 'Skipping pollutant "' // IBUF( 1:L1 ) // '"'
                CALL M3MSG2( MESG )
                CYCLE

            END IF

C.............  Read speciation profiles file

            CALL M3MSG2( '     Reading SPECIATION PROFILES file...' )

            CALL RDSPROF( RDEV, SINAM( V ), MXSPFUL, NSPFUL, NMSPC,
     &                    INPRF, SPECID, MOLEFACT, MASSFACT )

C.............  Initilialize multiple profiles and default reporting to true
            MULTIPRO = .TRUE.
            DEFREPRT = .TRUE.

C.............  When the number of profile table entries is the same as the
C               number of model species, then we know that there is only
C               one profile used for all sources, so do simple processing. The
C               one exception is when there is a pollutant-to-pollutant 
C               conversion, then we must still do the standard processing.
            IF( NSPFUL .EQ. NMSPC ) THEN

                L1 = LEN_TRIM( SINAM( V ) )
                L2 = LEN_TRIM( IBUF )
 
                IF( NMSPC .EQ. 1 ) THEN
                    MESG = '     NOTE: "' // IBUF( 1:L2 ) // 
     &                     '" only has a unit conversion ' //
     &                     'using profile "' // INPRF( 1 ) // '"'
                ELSE
                    MESG = '     NOTE: "' // IBUF( 1:L2 ) // 
     &                     '" is split for all sources ' //
     &                      'using profile "' // INPRF( 1 ) // '"'
                END IF

C.................  If there is no pollutant-to-pollutant conversion, then
C                   set speciation matrices using the one profile and continue
                IF( IBUF .EQ. SINAM( V ) ) THEN

                    IF( MASSOUT ) THEN
                        DO J = 1, NMSPC
                            MASSMATX( :,J ) = MASSFACT( J )
                        END DO
                    END IF

                    IF( MOLEOUT ) THEN
                        DO J = 1, NMSPC
                            MOLEMATX( :,J ) = MOLEFACT( J )
                        END DO
                    END IF

                    MULTIPRO = .FALSE.  ! no multiple profiles

C.................  Otherwise, need to continue so that pollutant-to-pollutant 
C                   conversion is done (so don't reset multipro)
                ELSE

                    L3 = LEN_TRIM( MESG )
                    MESG = MESG( 1:L3 ) // CRLF() // BLANK10 // 
     &                     'and a pollutant conversion to "' // 
     &                     SINAM( V )( 1:L1 ) // '"'

                    DEFREPRT = .FALSE.  ! no default reporting

                END IF  ! End EINAM .EQ. SINAM or not

                CALL M3MSG2( MESG )

            END IF      ! End single profile or nots

C.............  If this pollutant has multiple profiles...
            IF( MULTIPRO ) THEN

C.................  Abridge profiles so that there is an array of unique profiles
                CALL PROCSPRO( NMSPC, SPCNAMES( 1,V ) )

C.................  Assign speciation profile and populate speciation matrices
C                   for all sources for this pollutant.
                CALL ASGNSPRO( MASSOUT, MOLEOUT, DEFREPRT, NSRC, 
     &                         IBUF, MASSMATX, MOLEMATX )

C.................  Deallocate memory for unique profiles arrays
                DEALLOCATE( SPROFN, IDXSPRO, NSPECIES, IDXSSPEC )

            END IF      ! End multi-profile processing

C.............  Write out the speciation matrix for current pollutant

            IF( MASSOUT ) THEN

                MESG= '     Writing MASS-BASED SPECIATION MATRIX...'
                CALL M3MSG2( MESG )

                DO J = 1, NMSPC

                    CBUF = MASSONAM( J,V )
                    SBUF = SPCNAMES( J,V )
                    L1 = LEN_TRIM( CBUF )
                    L2 = LEN_TRIM( IBUF )
                    L3 = LEN_TRIM( SBUF )
                    L4 = LEN_TRIM( SNAME )
                    IF( .NOT. 
     &                  WRITE3( SNAME, CBUF, 0, 0, MASSMATX(1,J) )) THEN

                        EFLAG = .TRUE.

                        MESG = '     Could not write "' // 
     &                    IBUF( 1:L2 ) // '"-to-"' // SBUF( 1:L3 ) // 
     &                    '" speciation factor using name "' //
     &                    CBUF( 1:L1 ) // '" to file "' // 
     &                    SNAME( 1:L4 ) // '"'

                        CALL M3MSG2( MESG )
                        CYCLE

                    ELSE
                        MESG = BLANK10 // IBUF( 1:L2 ) // '-to-' // 
     &                         SBUF( 1:L3 ) // ' written to ' // 
     &                         SNAME( 1:L4 ) // ' as variable ' //
     &                         CBUF( 1:L1 ) 
                        CALL M3MSG2( MESG )

                    END IF

                END DO ! End write out of model species

            END IF    ! End mass-based output

            IF( MOLEOUT ) THEN

                MESG= '     Writing MOLE-BASED SPECIATION MATRIX...'
                CALL M3MSG2( MESG )

                 DO J = 1, NMSPC

                    CBUF = MOLEONAM( J,V )
                    IF( .NOT. 
     &                  WRITE3( LNAME, CBUF, 0, 0, MOLEMATX(1,J) )) THEN

                        EFLAG = .TRUE.

                        IBUF = EINAM( V )
                        SBUF = SPCNAMES( J,V )
                        L1 = LEN_TRIM( CBUF )
                        L2 = LEN_TRIM( IBUF )
                        L3 = LEN_TRIM( SBUF )
                        L4 = LEN_TRIM( LNAME )
                        MESG = 'Could not write "' // IBUF( 1:L2 ) // 
     &                         '"-to-"' // SBUF( 1:L3 ) // 
     &                         '" speciation factor using name "' //
     &                         CBUF( 1:L1 ) // '" to file "' // 
     &                         LNAME( 1:L4 ) // '"'

                        CALL M3MSG2( MESG )
                        CYCLE
 
                    END IF

                END DO ! End write out of model species

C.............  Write out file of speciation profiles used per source

! NOTE: future

C.............  Write out file of non-unity pollutant-to-pollutant conversion
C               factors for only those source with one such factor for a 
C               pollutant. NOTE: file contains source ID for easy matching to
C               inventory using QA programs.

! NOTE: future
C.............  Write speciation profile codes that were applied per source
C               to a file

            END IF    ! End mole-based output

        END DO     ! End loop through inventory pollutants

C.........  Check error flag for problems and end
        IF( EFLAG ) THEN

            MESG = 'Problem creating speciation matrices'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

c92000   FORMAT( 5X, A )

c92010   FORMAT ( 5X , A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

c93000   FORMAT( A )

c93010   FORMAT( A16 )

c93030   FORMAT( I5, I6, I4, 1X, A10, 1X, I5, F6.3, 15X, 3 F6.3 ) ! for ASREF


C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM SPCMAT

