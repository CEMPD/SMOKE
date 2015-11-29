
        SUBROUTINE RDTAG( FDEV, NOPOL )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the tagging file for any source category.  It
C     allocates memory (locally) for reading the unsorted x-refs. It sorts the
C     x-refs for processing. It allocates memory for the appropriate x-ref 
C     tables and populates the tables (passed via modules).
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened... MORE
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/2009 by M. Houyoux
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, CMACTA, ISPTA, 
     &                     CTAGNA, CSPCTAGNA, CISICA

C.........  This module contains the speciation profiles
        USE MODSPRO, ONLY: SPCLIST, NSPCALL, MXSPEC, SPCNAMES

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP

C.........  This module contains the tagging arrays 
        USE MODTAG, ONLY: MXTAG, NTAGSALL, TAGSPECIES, TAGNUM, TAGNAME

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         BLKORCMT
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        LOGICAL         USEEXPGEO

        EXTERNAL  BLKORCMT, CRLF, ENVINT, ENVYN, FINDC, GETFLINE,
     &            INDEX1, STR2INT, USEEXPGEO


C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! tagging file unit no.
        INTEGER, INTENT (IN) :: NOPOL  ! number of output pollutants

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 7
        INTEGER, PARAMETER :: NSPECIAL = 30
        CHARACTER(1), PARAMETER :: SPECIAL( NSPECIAL ) = 
     &                          ( / '~', '!', '@', '#', '$', '%', '^',
     &                              '&', '*', '(', ')', '-', '=', '+',
     &                              '{', '}', '[', ']', '|', '\', '<',
     &                              '>', ',', '.', '?', '/', ';', ':',
     &                              '"', '`' / )

        INTEGER, ALLOCATABLE :: SPCTAGCNT( : )  ! count of tags per SPCLIST species
        INTEGER, ALLOCATABLE :: POS1     ( : )  ! first position of species in sorted SPC-TAG list

C...........   Arrays for error and warning counters   
        INTEGER         IWRN(3)      ! warning counter
        INTEGER         IERR(5)      ! error counter

C...........   Array for parsing list-formatted inputs
        CHARACTER(50)          SEGMENT( MXCOL )

C...........   Other local variables
        INTEGER         C, I, J, K, L, N, P, V    !  counters and indices

        INTEGER         IDUM    !  dummy integer
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         JSPC    !  position of CSPC in SPCLIST
        INTEGER         MXERR   !  max no. errors of each type
        INTEGER         MXWARN  !  max no. warnings of each type
        INTEGER         NLINES  !  number of lines
        INTEGER         NXREF   !  number of valid x-ref entries
        INTEGER         POS2    !  tmp second position counter

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found
        LOGICAL      :: PFLAG = .FALSE.   !  true: species-specific record skipped
        LOGICAL      :: SKIPREC = .FALSE. !  true: record skipped in x-ref file
 
        CHARACTER(5)       SPOS     !  temporary pol code or position
        CHARACTER(IOVLEN3) CBUF     !  generic character string 
        CHARACTER(FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(MACLEN3) CMCT     !  temporory MACT code
        CHARACTER(SICLEN3) CSIC     !  temporary SIC
        CHARACTER(ALLLEN3) CSRCALL  !  buffer for source char, incl pol
        CHARACTER(IOVLEN3) CSPC     !  temporary species name
        CHARACTER(IOVLEN3+TAGLEN3+1 ) CSPCTAG  ! tmp species-tag
        CHARACTER(TAGLEN3) CTAG     !  temporary tag name
        CHARACTER(PLTLEN3) PLT      !  tmp plant ID
        CHARACTER(IOVLEN3) PSPC     !  previous iteration tmp species name
        CHARACTER(IOVLEN3) PTAG     !  previous iteration tmp tag name
        CHARACTER(SCCLEN3) SCCZERO  !  buffer for zero SCC
        CHARACTER(SICLEN3) SICZERO  !  buffer for zero SIC
        CHARACTER(SCCLEN3) TSCC     !  temporary SCC

        CHARACTER(512)     LINE     !  line buffer
        CHARACTER(512)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDTAG' ! program name

C***********************************************************************
C   begin body of subroutine RDTAG

C.........  Write status message
        MESG = 'Reading Tagging file...'
        CALL M3MSG2( MESG )

C.........  Get the count for max errors and warnings
        MXWARN = ENVINT( WARNSET , ' ', 100, I )
        MXERR  = ENVINT( ERRSET  , ' ', 100, I )
        IWRN = 0     ! array
        IERR = 0     ! array

C.........  Set up zero strings for SCC code and SIC code
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SICZERO = REPEAT( '0', SICLEN3 )

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Tagging file' )

C.........  Allocate memory for unsorted data used in all source categories 
        ALLOCATE( CTAGNA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTAGNA', PROGNAME )
        ALLOCATE( CSPCTAGNA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPCTAGNA', PROGNAME )
        ALLOCATE( ISPTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPTA', PROGNAME )
        ALLOCATE( CSCCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        ALLOCATE( CMACTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMACTA', PROGNAME )
        ALLOCATE( CISICA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CISICA', PROGNAME )
        ALLOCATE( INDXTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        CTAGNA    = ' ' ! array
        CSPCTAGNA = ' ' ! array
        CSCCTA    = ' ' ! array
        CSRCTA    = ' ' ! array
        CMACTA    = ' ' ! array
        CISICA    = ' ' ! array
        ISPTA     = 0   ! array
        INDXTA    = 0   ! array

C.........  Allocate memory for list of tagged species and initialize
        ALLOCATE( TAGSPECIES( NSPCALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TAGSPECIES', PROGNAME )
        TAGSPECIES = ' '   ! array
        NTAGSALL   = 0

C.........  Put file read pointer at top of file
        REWIND( FDEV )

C.........  Initialize character strings
        SEGMENT = ' '  ! array

C.........  Read lines and store unsorted data for the source category of 
C           interest
        IREC   = 0
        N      = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading tagging file at line', IREC
                CALL M3MSG2( MESG )
                CYCLE
            END IF

C.............  Skip blank lines or comments
            IF( BLKORCMT( LINE ) ) CYCLE

            CALL PARSLINE( LINE, MXCOL, SEGMENT )

            CFIP   = SEGMENT( 1 )
            TSCC   = SEGMENT( 2 )
            CSPC   = SEGMENT( 3 )
            CTAG   = SEGMENT( 4 )                    
            CMCT   = SEGMENT( 5 )
            CSIC   = SEGMENT( 6 )
            PLT    = SEGMENT( 7 )

C.............  Check that tag name is <= TAGLEN3 (8 character limit)
            L = LEN( TRIM( SEGMENT( 4 ) ) ) 
            IF ( L > TAGLEN3 ) THEN
                EFLAG = .TRUE.
                IF ( IERR(1) .LE. MXERR ) THEN
                    WRITE( MESG,94010 ) 'ERROR: The tag "'//
     &                  TRIM(SEGMENT(4))//'" provided at line', IREC,
     &                 'is longer than the ', TAGLEN3, 
     &                 '-character limit'
                    CALL M3MSG2( MESG )
                END IF
                IERR(1) = IERR(1) + 1
            END IF

C.............  Check that species name  + tag is <= IOVLEN3 (16 character limit)
            L = LEN( TRIM( CSPC ) // TRIM( CTAG ) )
            IF ( L > IOVLEN3 ) THEN
                EFLAG = .TRUE.
                IF ( IERR(2) .LE. MXERR ) THEN
                    WRITE( MESG,94010 ) 'ERROR: The tag "'//TRIM(CTAG)//
     &                 '" provided at line', IREC, 'is longer than '//
     &                 'the', IOVLEN3, '-character limit when '//
     &                 'combined with species "'// TRIM( CSPC )// '"'
                    CALL M3MSG2( MESG ) 
                END IF
                IERR(2) = IERR(2) + 1
            END IF 

C.............  Check no special characters in tags
C.............  Intentially do not use the error counter here so that users
C               will get all of the problems on the first run.
            DO K = 1, NSPECIAL
                L = INDEX( CTAG, SPECIAL( K ) )

                IF ( L > 0 ) THEN
                    EFLAG = .TRUE.
                    IF ( IERR(3) .LE. MXERR ) THEN
                      WRITE( MESG,94010 ) 'ERROR: The tag "'//TRIM(CTAG)
     &                 //'" contains the special character "' //
     &                 SPECIAL(N) // '" at line', IREC, ', which '//
     &                 'is not permitted in I/O API variable names.'
                      CALL M3MSG2( MESG ) 
                    END IF
                    IERR(3) = IERR(3) + 1
                END IF
            END DO

C.............  Skip all point entries for nonpoint sectors
            IF( PLT == '-9' ) PLT = ' '   ! reset -9 to blank
            IF( CATEGORY /= 'POINT' .AND. PLT /= ' ' ) CYCLE 

C.............  Post-process x-ref information to scan for '-9', pad
C               with zeros, compare SCC version master list, and compare
C               pollutant/emission type name with master list.
C.............  JSPC is getting set to 0 since CBUF = ' '
            CBUF = ' '
            CALL FLTRXREF( CFIP, CSIC, TSCC, CBUF, CMCT, 
     &                     IDUM, IDUM, JSPC, PFLAG, SKIPREC )
     
C.................  Filter the case where the species code is not present
            IF( CSPC .EQ. ' ' .OR. CSPC == '-9' ) THEN
                EFLAG = .TRUE.
                IF( IERR(4) .LE. MXERR ) THEN
                    WRITE( MESG, 94010 ) 
     &                'ERROR: Bad tagging entry at line', IREC,
     &                'because of missing species.'
                    CALL M3MSG2( MESG )
                END IF
                IERR(4) = IERR(4) + 1
            END IF

C.............  Check to see if FIPS code matches inventory.  If no match,
C               then give a warning and skip the record from processing.
            IF ( USEEXPGEO() .OR. 
     &           CFIP( FIPEXPLEN3+4:FIPEXPLEN3+6) /= '000' ) THEN
                J = FINDC( CFIP, NINVIFIP, INVCFIP )
                IF ( J .LE. 0 ) THEN
                    IF( IWRN(1) .LE. MXWARN ) THEN 
                        WRITE( MESG, 94010 ) 
     &                   'WARNING: State/county FIPS code "'//CFIP//
     &                   '" at line', IREC, 'does not match inventory.'
                        CALL M3MSG2( MESG )
                        IWRN(1) = IWRN(1) + 1
                    END IF
                    SKIPREC = .TRUE.
                    CALL M3MSG2( MESG ) 
                END IF
            END IF

            IF( SKIPREC ) CYCLE  ! Skip this record

C.................  Set position of species in species list
            JSPC = INDEX1( CSPC, NSPCALL, SPCLIST )
            IF ( JSPC .LE. 0 ) THEN
                IF( IWRN(2) .LE. MXWARN ) THEN
                    WRITE( MESG, 94010 )'WARNING: Species "'//TRIM(CSPC)
     &                 //'" at line',IREC, 'is not a valid species '//
     &                 'for the pollutants and species'// CRLF() //
     &                 BLANK10 // 'provided by the inventory and ' //
     &                 'speciation profiles'
                    CALL M3MSG2( MESG )
                END IF
                IWRN(2) = IWRN(2) + 1
                CYCLE  ! skip
            END IF

C.................  Write species position to character string
            WRITE( SPOS, '(I5.5)' ) JSPC  

C.................  Check for bad tag code
            IF( CTAG == ' ' .OR. CTAG == '-9' ) THEN
                EFLAG = .TRUE.
                IF( IERR(5) .LE. MXERR ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'ERROR: Missing tag in tagging file at line',
     &                   IREC
                    CALL M3MSG2( MESG )
                END IF
                IERR(5) = IERR(5) + 1
            END IF

C.............  If SIC is defined, make sure SCC is not and fill SCC
C               with SIC value and special identifier
            IF( CSIC /= SICZERO .AND. TSCC /= SCCZERO ) THEN

                IF( IWRN(3) .LE. MXWARN ) THEN 
                    WRITE( MESG,94010 ) 'WARNING: Both SCC and SIC ' //
     &                 'values are given at line', IREC, '.'//CRLF()//
     &                 BLANK10 // 'Only the SCC will be used ' //
     &                 'for this tagging entry.'
                    CALL M3MSG2( MESG )
                    IWRN(3) = IWRN(3) + 1
                END IF
                CSIC = SICZERO

            END IF

C.............  If error found, no point in storing anything, so skip
            IF ( EFLAG ) CYCLE

C.............  Increment count of valid x-ref entries and check it
            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

C.................  Store case-specific fields from tagging file
            CSRCALL = ' '
            SELECT CASE( CATEGORY )

            CASE( 'AREA' )

                CALL BLDCSRC( CFIP, TSCC, CHRBLNK3,
     &                        CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                        CHRBLNK3, POLBLNK3, CSRCALL   )

                CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // CMCT // CSIC // SPOS

            CASE( 'MOBILE' )

                CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,
     &                        CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                        POLBLNK3, CSRCALL )

                CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // CMCT // CSIC // SPOS

            CASE( 'POINT' )

C.................  Store sorting criteria as right-justified in fields
                CALL BLDCSRC( CFIP, PLT, CHRBLNK3,
     &                        CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                        CHRBLNK3, POLBLNK3, CSRCALL   )

                CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC // 
     &                        CMCT // CSIC // SPOS

            END SELECT

C.................  Store case-indpendent fields from tagging file
            INDXTA( N )   = N
            ISPTA ( N )   = JSPC    ! Save index to master species list
            CSCCTA( N )   = TSCC
            CMACTA( N )   = CMCT
            CISICA( N )   = CSIC
            CTAGNA( N )   = CTAG
            CSPCTAGNA( N )= TRIM(CSPC)//'-'//CTAG  ! join using dash because it has been filtered out of species and tag names since it's invalid for I/O API

C.............  Add species to TAGSPECIES list, if needed
            IF( NTAGSALL == 0 ) THEN
                NTAGSALL = 1
                TAGSPECIES( 1 ) = CSPC
            ELSE
                K = INDEX1( CSPC, NTAGSALL, TAGSPECIES )
                IF ( K <= 0 ) THEN
                    NTAGSALL = NTAGSALL + 1
                    TAGSPECIES( NTAGSALL ) = CSPC
                END IF
            END IF

        END DO      ! End of loop on I for reading in speciation x-ref file


C.........  Reset number of tagging file entries in case some were dropped
        NXREF = N

C.........  Write error if no tagging entries match inventory
        IF( .NOT. EFLAG .AND. NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No tagging file entries match a source!'
            CALL M3MSG2( MESG )

        ELSE IF( NXREF .GT. NLINES ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing tagging data was', NLINES,
     &             CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )

        END IF

C.......  Check for errors reading XREF file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading tagging file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Processing tagging file...'
        CALL M3MSG2( MESG )

C.........  Sort speciation cross-reference entries.  Unlike other cross-
C           reference processing, the position of the species is included
C           using "SPOS" rather than the position of the pollutant.  Since
C           all entries must be species-specific, there will be no entries
C           that are SPOS = 0.
        CALL SORTIC( NXREF, INDXTA, CSRCTA )

        CALL XREFTBL( 'TAGGING', NXREF )

C See Madeleine's document to make sure I've put in all of the checks

C.........  Steps to get the maximum number of tags for any species (MXTAG)
C.........  Reset the sorting index
        DO I = 1, NXREF
            INDXTA( I ) = I
        END DO

C.........  Sort the data based on species then tag
        CALL SORTIC( NXREF, INDXTA, CSPCTAGNA )

        ALLOCATE( SPCTAGCNT( NSPCALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCTAGCNT', PROGNAME )
        ALLOCATE( POS1( NSPCALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POS1', PROGNAME )
        SPCTAGCNT = 0  ! array
        POS1      = 0  ! array

C.........  Count and store the number of tags for each species, SPCTAGCNT(:)
        PSPC = ' '
        PTAG = ' '
        DO I = 1, NXREF
            J = INDXTA( I )

            CSPCTAG = CSPCTAGNA( J ) 
            K = INDEX( CSPCTAG, '-' )
            L = LEN_TRIM( CSPCTAG )
            CSPC = CSPCTAG( 1:K-1 )
            CTAG = CSPCTAG( K+1:L )

            IF( CSPC /= PSPC ) THEN
                 N = INDEX1( CSPC, NSPCALL, SPCLIST )
                 PTAG = ' '
            END IF
            IF( CTAG /= PTAG ) SPCTAGCNT( N ) = SPCTAGCNT( N ) + 1
            IF( POS1( N ) == 0 ) POS1( N ) = I

	    PSPC = CSPC
	    PTAG = CTAG
        END DO

C.........  Get the maximum number of tags for any species, MXTAG
        MXTAG = MAXVAL( SPCTAGCNT )  ! array

C.........  Deallocate TAGNAME and re-allocate based on MXTAG
        DEALLOCATE( TAGNUM, TAGNAME )
        ALLOCATE( TAGNUM( MXSPEC, NOPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TAGNUM', PROGNAME )
        ALLOCATE( TAGNAME( 0:MXTAG, MXSPEC, NOPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TAGNAME', PROGNAME )
        TAGNUM  = 0   ! array
        TAGNAME = ' '   ! array

C.........  Look up species names in SPCLIST and use count of tags 
C           per species to assign TAGNUM and populate TAGNAME
        DO V = 1, NOPOL
            DO P = 1, MXSPEC
                IF ( SPCNAMES( P,V ) == ' ' ) CYCLE
                N = INDEX1( SPCNAMES( P,V ), NSPCALL, SPCLIST )
                TAGNUM( P,V ) = SPCTAGCNT( N )

C.................  If skip to next species if no tags
                IF( TAGNUM( P,V ) == 0 ) CYCLE

C.................  Loop through list of spc/tags and populate TAGNAME
                C = 0
                I = POS1( N )
                PTAG = ' '
                DO
                    J = INDXTA( I )
                    CSPCTAG = CSPCTAGNA( J )
                    K = INDEX( CSPCTAG, '-' )
                    L = LEN_TRIM( CSPCTAG )
                    CSPC = CSPCTAG( 1:K-1 )
                    CTAG = CSPCTAG( K+1:L )
                    IF ( CTAG /= PTAG .OR. CSPC /= PSPC ) THEN
                        C = C + 1
                        TAGNAME( C,P,V ) = CTAG

                        PTAG = CTAG
                	PSPC = CSPC
                    END IF 

                    I = I + 1

C.....................  Stop trying to collect the tag names if got them all
                    IF ( C == SPCTAGCNT( N ) ) EXIT

                END DO

            END DO
        END DO

C.........  Deallocate other temporary unsorted arrays
        DEALLOCATE( CSCCTA, ISPTA, CMACTA, CSRCTA, CTAGNA, INDXTA )
        DEALLOCATE( SPCTAGCNT, POS1 )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of tagging file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDTAG
