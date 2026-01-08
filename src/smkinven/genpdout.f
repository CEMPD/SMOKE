
        SUBROUTINE GENPDOUT( FDEV, CDEV, ODEV, RDEV,TZONE, SDATE, STIME, 
     &                       NSTEPS, INSTEP, OUTSTEP, NVAR, NVSP, 
     &                       MXPDSRC, TYPNAM, FNAME, EAIDX, SPIDX, CFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads and writes the day-specific or hour-specific
C      emissions.  It also write a report file for CEM formatted hour-specific
C      data in which the sources have been matched by ORIS ID.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C      09/2025 by HT UNC-IE:  Use M3UTILIO; 
C                             Removed MESG FMT 93000 and 9400 which were not used;
C                             replace NAMLEN3 with IOULEN3 for ONAME
C
C*************************************************************************
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
        USE M3UTILIO

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, CSCC, CPDESC, CORIS

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVSCC, NINVORIS, INVORIS, INVSCC,
     &                      IORSMTCH, SCCDESC, INVODSC, INVORFP,FIREFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NCHARS, SC_BEGP, SC_ENDP, NSRC,
     &                     NCOMP, VAR_FORMULA

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT, NPDPT, NPDPTP, CODEA, CIDXA, IDXSRC,
     &                      SPDIDA, EMISVA, DYTOTA, LPDSRC,
     &                      PDEMOUT, PDTOTL, NUNFDORS, UNFDORS

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NORIS, ORISLST, ORISFIP, ORISDSC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
c       INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
c       CHARACTER(2) CRLF
c       LOGICAL      ENVYN
c       INTEGER      FINDC
c       INTEGER      INDEX1

c       EXTERNAL     CRLF, ENVYN, FINDC, INDEX1

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV      ! hour-specific file unit no.
        INTEGER     , INTENT (IN) :: CDEV      ! SCC desc file unit no.
        INTEGER     , INTENT (IN) :: ODEV      ! ORIS desc file unit no.
        INTEGER     , INTENT (IN OUT) :: RDEV  ! Report file unit no.
        INTEGER     , INTENT (IN) :: TZONE     ! output time zone
        INTEGER     , INTENT (IN) :: SDATE     ! Julian starting date in TZONE
        INTEGER     , INTENT (IN) :: STIME     ! start time of data in TZONE
        INTEGER     , INTENT (IN) :: NSTEPS    ! no. time steps
        INTEGER     , INTENT (IN) :: INSTEP    ! expected data time step HHMMSS
        INTEGER     , INTENT (IN) :: OUTSTEP   ! output time step HHMMSS
        INTEGER     , INTENT (IN) :: NVAR      ! no. period-specific variables
        INTEGER     , INTENT (IN OUT) :: NVSP  ! no. period-spec special vars
        INTEGER     , INTENT (IN) :: MXPDSRC   ! maximum period-specific sources
        CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'
        CHARACTER(*), INTENT (IN) :: FNAME     ! logical file name
        INTEGER     , INTENT (IN) :: EAIDX( NIPPA ) ! index to EANAM
        INTEGER     , INTENT (IN) :: SPIDX( MXSPDAT ) ! index to SPDATNAM
        LOGICAL     , INTENT (IN) :: CFLAG     ! CEM processing

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: EASTAT( : )    ! true: act/pol present in data
        CHARACTER(SCCLEN3), ALLOCATABLE :: ELECSCC( : )

C.........  Local arrays
        INTEGER         SPSTAT( MXSPDAT )     ! true: special data variable used
        LOGICAL         LFG( 9 )          ! true: source characteristic is valid

        CHARACTER(15)   CHRHDRS( NCHARS )     ! Source characteristics headers
        CHARACTER(50)   CHARS( 9 )
        CHARACTER(40)   LABEL( 2 )

C...........   Other local variables
        INTEGER          I, J, K, L, N, S, T

        INTEGER          INVFMT               ! format code of files in list
        INTEGER          IOS                  ! i/o status
        INTEGER          JDATE                ! tmp Julian date
        INTEGER          JTIME                ! tmp time HHMMSS
        INTEGER          NELECSCC             ! number electric generating SCCs
        INTEGER          NPDSRC               ! number of day/hour-spec sources
        INTEGER          OWID                 ! width of ORIS desc field
        INTEGER          PDEMDIM              ! dim for PDEMOUT
        INTEGER          PWID                 ! width of plant desc field
        INTEGER          WID                  ! width of field

        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific
        LOGICAL       :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL       :: LASTSTEP = .FALSE.  ! true: last time step when calling wrpdemis
        LOGICAL, SAVE :: OFLAG    = .FALSE.  ! true: PFLAG & hourly
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: create hourly profiles
        LOGICAL, SAVE :: SFLAG    = .FALSE.  ! true: create daily totals

        CHARACTER(300) :: MESG = ' '         ! message buffer
        CHARACTER(256) :: FMTBUF = ' '       ! format buffer
        CHARACTER(256) :: FMTBUFB = ' '      ! format buffer B

        CHARACTER(FIPLEN3) CFIP          ! tmp co/st/cy code
        CHARACTER(SDSLEN3) BUFFER        ! tmp SCC description
        CHARACTER(IOULEN3) ONAME         ! output file name
        CHARACTER(SCCLEN3) TSCC          ! tmp SCC value
        CHARACTER(ORSLEN3) CORS          ! tmp ORIS ID
        CHARACTER(DSCLEN3) PDSC          ! tmp plant DSC
        CHARACTER(DSCLEN3) ODSC          ! tmp ORIS plant DSC
         
        CHARACTER(16) :: PROGNAME = 'GENPDOUT' !  program name

C***********************************************************************
C   begin body of program GENPDOUT

C.........  For the first time the routine is called...
        IF( FIRSTIME ) THEN

C.............  Get environment variable for creating daily data from the hourly
            MESG = 'Use daily totals only from hourly data file'
            SFLAG = ENVYN( 'HOURLY_TO_DAILY', MESG, .FALSE., IOS )

C.............  Get environment variable for creating hourly profiles from the
C               hourly data
            MESG = 'Create hourly profiles from hourly data'
            PFLAG = ENVYN( 'HOURLY_TO_PROFILE', MESG, .FALSE., IOS )

            IF( SFLAG .AND. PFLAG ) THEN
                MESG = 'WARNING: Ignoring HOURLY_TO_PROFILE "Y" ' //
     &                 'value because HOURLY_TO_DAILY is set to "Y"'
                CALL M3MSG2( MESG )
                PFLAG = .FALSE.

            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  Perform case-specific settings
        OFLAG = .FALSE.
        SELECT CASE( TYPNAM )
        CASE( 'day' ) 
            DFLAG = .TRUE.

        CASE( 'hour' )
            IF( PFLAG ) OFLAG = .TRUE.
            DFLAG = .FALSE.

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Allocate memory for logical status array for pol/act even though
C           it does not need to be set because EAIDX has already been 
C           determined.
        ALLOCATE( EASTAT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EASTAT', PROGNAME )
        EASTAT = 0  ! array
        SPSTAT = 0  ! array

C.........  Allocate memory for reading data
        ALLOCATE( MXPDPT( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MXPDPT', PROGNAME )
        ALLOCATE( NPDPT ( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPDPT', PROGNAME )
        ALLOCATE( CODEA ( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CODEA', PROGNAME )
        ALLOCATE( CIDXA ( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CIDXA', PROGNAME )
        ALLOCATE( IDXSRC( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRC', PROGNAME )
        ALLOCATE( SPDIDA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDIDA', PROGNAME )
        ALLOCATE( EMISVA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISVA', PROGNAME )
        ALLOCATE( DYTOTA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DYTOTA', PROGNAME )
        ALLOCATE( LPDSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LPDSRC', PROGNAME )

C.........  Initialize arrays
        MXPDPT = 0        ! array
        NPDPT  = 0        ! array
        CODEA  = 0        ! array
        CIDXA  = 0        ! array
        IDXSRC = 0        ! array
        SPDIDA = 0        ! array
        EMISVA = BADVAL3  ! array
        DYTOTA = BADVAL3  ! array
        LPDSRC = .FALSE.  ! array

C.........  Initialize special arrays for fires processing
        IF ( FIREFLAG ) THEN
            ALLOCATE( NPDPTP( NSTEPS,NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NPDPTP', PROGNAME )
            NPDPTP = 0        ! array
        END IF

C.........  Message before reading the input file (list of files)
        MESG = 'Reading ' // TYPNAM // '-specific data...'
        CALL M3MSG2( MESG )

C.........  Loop through input files and actually read the data
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, INVFMT, 
     &                 EASTAT, SPSTAT )

C.........  Determine the actual number of day-specific or hour-specific sources
C.........  If fire-data processing, use a different approach to determine the
C           sources, since the fire data are so sparse.  This prevents the PDAY
C           files from containing a ton of empty records at the end of
C           each time step (since the PDAY files are sparsely stored with a
C           source index that changes as needed for each time step)
        IF ( FIREFLAG ) THEN
            NPDSRC = MAXVAL( NPDPTP ) ! array
        ELSE
            NPDSRC = 0
            DO S = 1, NSRC
                IF( LPDSRC( S ) ) NPDSRC = NPDSRC + 1
            END DO
        END IF

C.........  Make sure that that actual number of sources over all sources does
C           not exceed the maximum number of sources over all hours
C        IF( NPDSRC .GT. MXPDSRC .AND. .NOT. FIREFLAG ) THEN
C            WRITE( MESG,94010 ) 'INTERNAL ERROR: Actual number of ' //
C     &             TYPNAM // 'sources, NPDSRC=', NPDSRC, CRLF() // 
C     &             BLANK10 // 'dimensioned number, MXPDSRC =', MXPDSRC,
C     &             '. Fix by ensuring ALL period-specific' // CRLF() //
C     &             BLANK10 // 'sources in file for ANY day or hour ' //
C     &             'have at least one entry ' // CRLF() // BLANK10 //
C     &             'for the SAME day or hour.'
C            CALL M3MSG2( MESG )
C            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        IF( NPDSRC .EQ. 0 ) THEN

            MESG = 'No period-specific sources found in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Open day-specific or hour-specific output file, depending on value
C           of TYPNAM
        CALL OPENPDOUT( NPDSRC, NVAR, TZONE, SDATE, STIME, OUTSTEP, 
     &                  INVFMT, TYPNAM, OFLAG, EAIDX, SPSTAT, 
     &                  ONAME, RDEV )

C.........  Allocate memory for daily or hourly output arrays.  Allocate 
C           memory as one block which will be separated into an integer section
C           and a real section when WRPDEMIS is called.  This permits
C           writing with a single WRITE3 statement.
        ALLOCATE( PDEMOUT( NPDSRC,NVSP+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PDEMOUT', PROGNAME )
        ALLOCATE( PDTOTL( NPDSRC,NVSP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PDTOTL', PROGNAME )

C.........  Loop through time steps and output emissions and other data

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS
        
            LASTSTEP = ( T .EQ. NSTEPS ) 
            CALL WRPDEMIS( DFLAG, JDATE, JTIME, T, NPDSRC, NVAR, NVSP, 
     &                     ONAME, OFLAG, CFLAG, EAIDX, SPIDX, LASTSTEP,
     &                     PDEMOUT( 1,1 ), PDEMOUT( 1,2 ), EFLAG, MXPDSRC )

            CALL NEXTIME( JDATE, JTIME, OUTSTEP )

        END DO     !  End of loop over time steps

C.............  Abort if error found
        IF ( EFLAG ) THEN
            MESG = 'Problem with input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate local memory
        DEALLOCATE( EASTAT )

C.........  Deallocate global memory
        DEALLOCATE( MXPDPT, NPDPT, CODEA, CIDXA, IDXSRC, SPDIDA, EMISVA, 
     &              DYTOTA, LPDSRC, PDEMOUT, PDTOTL )

C.........  Exit from subroutin if not writing CEM report...
        IF ( INVFMT .NE. CEMFMT ) RETURN

C.........  Read SCC descriptions
        CALL RDSCCDSC( CDEV )

C.........  Read ORIS descriptions
        CALL RDORSDSC( ODEV )

C.........  Get maximum width of description fields
        PWID = 29    ! Size of header
        DO I = 1, NINVORIS
            L = LEN_TRIM( INVODSC( I ) )
            PWID = MAX( PWID, L )
        END DO

        OWID = 24    ! Size of header
        DO I = 1, NORIS
            L = LEN_TRIM( ORISDSC( I ) )
            OWID = MAX( OWID, L )
        END DO

C.........  Create report of inventory ORIS IDs that matched the CEM data...

C.........  Create header:
        WRITE( FMTBUF, '(A,I3.3,A,I3.3,A)' ) 
     &      '("Inventory ORIS IDs that matched the CEM data",/,A,/,' //
     &      '"ORIS ID; Region; ",A', PWID, '"; ",A', OWID, ')'

        WID = 21 + PWID + OWID
        LABEL(1) = 'Inventory Plant Description'
        LABEL(2) = 'ORIS Plant Description'
        WRITE( RDEV, FMTBUF ) REPEAT( '-', WID ), LABEL(1), LABEL(2)

C.........  Create content
        WRITE( FMTBUF, '(A,I1,A,I2.2,A,I3.3,A,I3.3,A)' ) 
     &         '(A', MAX(7,ORSLEN3), ',"; ", A', FIPLEN3, 
     &         ', "; ", A', PWID, ', "; ", A', OWID, ')'

        DO I = 1, NINVORIS

            IF( IORSMTCH( I ) ) THEN
                CORS = INVORIS( I )
                CFIP = INVORFP( I )
                PDSC = INVODSC( I )
                J = INDEX1( CORS, NORIS, ORISLST )
                IF ( J .GT. 0 ) THEN
                    ODSC = ORISDSC( J )
                ELSE
                    ODSC = 'NOT AVAILABLE'
                END IF

                WRITE( RDEV,FMTBUF ) CORS, CFIP, PDSC, ODSC
            END IF

        END DO

C.........  Create report that lists CEM ORISs that were not in the 
C           inventory....

C.........  Create header:
        WRITE( FMTBUF, '(A,I3.3,A,I3.3,A)' ) 
     &      '(2/,"CEM ORIS IDs that did not match the inventory",/,A,/,'
     &      // '"ORIS ID; Region; ",A', OWID, ')'

        WID = 21 + OWID
        LABEL(1) = 'ORIS Plant Description'
        WRITE( RDEV, FMTBUF ) REPEAT( '-', WID ), LABEL(1)

C.........  Create content
        WRITE( FMTBUF, '(A,I1,A,I2.2,A,I3.3,A)' ) 
     &         '(A', MAX(7,ORSLEN3), ',"; ", A', FIPLEN3, 
     &         ', "; ", A', OWID, ')'

        WRITE( FMTBUFB, '(A,I1,A,I2.2,A,I3.3,A)' ) 
     &         '(A', MAX(7,ORSLEN3), ',"; ", A', FIPLEN3,  
     &         ', "; ", A', OWID, ')'

        DO I = 1, NUNFDORS 

            CORS = UNFDORS( I )

            J = INDEX1( CORS, NORIS, ORISLST )

            IF ( J .GT. 0 ) THEN
                CFIP = ORISFIP( J )
                ODSC = ORISDSC( J )
                WRITE( RDEV, FMTBUF ) CORS, CFIP, ODSC
            ELSE
                CFIP = ' '
                ODSC = 'NOT AVAILABLE'
                WRITE( RDEV, FMTBUFB ) CORS, CFIP, ODSC
            END IF

        END DO

C.........  Create list of powerplant SCCs...

C.........  Create header:
        CHRHDRS( 2 ) = 'Plant ID'
        IF ( NCHARS .GE. 3 ) CHRHDRS( 3 ) = 'Char 1'
        IF ( NCHARS .GE. 4 ) CHRHDRS( 4 ) = 'Char 2'
        IF ( NCHARS .GE. 5 ) CHRHDRS( 5 ) = 'Char 3'
        IF ( NCHARS .GE. 6 ) CHRHDRS( 6 ) = 'Char 4'
        IF ( NCHARS .GE. 7 ) CHRHDRS( 7 ) = 'Char 5'

        WRITE( FMTBUF, '(A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A)' ) 
     &      '(2/,"Power generating sources in inventory that ' //
     &      'did not match CEM",/,A,/,"Region; ",A', PLTLEN3,
     &      ',"; ",', NCHARS-2, '(A',CHRLEN3,',:,"; "),A',
     &       SCCLEN3, ',"; ",A', DSCLEN3, ')'

        WID = FIPLEN3 + PLTLEN3 + (NCHARS-2)*CHRLEN3 + SCCLEN3 + DSCLEN3
        LABEL(1) = 'SCC'
        LABEL(2) = 'Plt Name'
        WRITE( RDEV, FMTBUF ) REPEAT( '-', WID ), 
     &       ( CHRHDRS( I ), I=2, NCHARS ), LABEL(1), LABEL(2)

C.........  Create content
        IF( ALLOCATED( ELECSCC ) ) DEALLOCATE( ELECSCC )
        ALLOCATE( ELECSCC( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELECSCC', PROGNAME )
        ELECSCC = ' '  ! array

        K = 0
        DO I = 1, NINVSCC

            BUFFER = SCCDESC( I )
            CALL UPCASE( BUFFER )
            J = INDEX( BUFFER, 'ELECTRIC GENERATION' )
            IF ( J .GT. 0 ) THEN
                K = K + 1
                ELECSCC( K ) = INVSCC( I )
            END IF

        END DO
        NELECSCC = K

C.........  Set logical array for setting valid source characeristics columns
        LFG( 1:NCHARS ) = .TRUE.   ! array
        IF( NCHARS .LE. 8 ) LFG( NCHARS+1:9 ) = .FALSE.  ! array

C.........  Create output format
        WRITE( FMTBUF, 93042 ) FIPLEN3, PLTLEN3, NCHARS-2, 
     &                         CHRLEN3, SCCLEN3, DSCLEN3

C.........  Create report of inventory ORIS IDs with powerplant SCCs that were
C           not in the CEM database.
        DO S = 1, NSRC

            TSCC = CSCC( S )
            CORS = CORIS( S )
            PDSC = CPDESC( S )

            I = FINDC( TSCC, NELECSCC, ELECSCC ) 
            IF ( I .LE. 0 ) CYCLE

            I = FINDC( CORS, NINVORIS, INVORIS )

            CALL PARSCSRC( CSOURC( S ), NCHARS, SC_BEGP,
     &                     SC_ENDP, LFG, N, CHARS )

C.............  If source is in list of inventory ORIS IDs
            IF ( I .GT. 0 ) THEN

C.................  Check to see if a match was found in the CEM data
                IF( .NOT. IORSMTCH( I ) ) THEN
                    WRITE( RDEV, FMTBUF ) 
     &                   ( CHARS( I ), I = 1,NCHARS ), TSCC, PDSC 

                END IF

C.............  If not in list, the power generating SCC didn't have ORIS ID
C               in the inventory.
            ELSE
                WRITE( RDEV, FMTBUF ) 
     &                   ( CHARS( I ), I = 1,NCHARS ), TSCC, PDSC

            END IF

        END DO
 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93042   FORMAT( '( A', I2.2,', "; ", A', I2.2, ', "; ",', I2.2, '(A',
     &          I2.2, ',"; "),', 'A', I2.2, ', "; ", A', I2.2,')' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GENPDOUT
