
        SUBROUTINE RDVMIX( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the vehicle mix file, sorts it,  calls the appropriate routines
C     to group the sorted data
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened... MORE
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/2000 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: NVTYPE, CVTYPLST, VMTMIXA

C...........   This module is for cross reference tables
        USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         GETNLIST
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL  CHKINT, CRLF, ENVYN, FINDC, GETNLIST, GETFLINE, 
     &            INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV              ! VMT mix file unit no.

C...........   Local parameters
        INTEGER, PARAMETER :: NPRECOLL = 3  ! list-directed no. cols before mix
        INTEGER, PARAMETER :: NPRECOLF = 6  ! fixed format no. cols before mix
        INTEGER, PARAMETER :: FIXPREWD = 22 ! fixed format width before VMT mix
        INTEGER, PARAMETER :: FIXVMTWD = 5  ! width of VMT mix columns

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: VIDX ( : )    ! vmix index to master list

        REAL, ALLOCATABLE :: VMIX ( : )    ! temporary vehicle mix

        CHARACTER(20), ALLOCATABLE :: SEGMENT( : )  ! Segments of parsed lines

        CHARACTER(VTPLEN3), ALLOCATABLE :: VTNAMES( : ) ! veh type names

C...........   Other local variables
        INTEGER         I, J, K, L, L2, N, P1, P2, V    !  counters and indices

        INTEGER         CYID    !  tmp county ID
        INTEGER         FIP     !  temporary country/state/county code
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         LINTYPE !  temporary source category code
        INTEGER         MXTCOL  !  maximum number of table columns
        INTEGER         NLINES  !  number of lines
        INTEGER         NV      !  no. vehicle types in VMIX file
        INTEGER         NREF    !  number of x-ref entries before filtering
        INTEGER         NXREF   !  number of valid x-ref entries
        INTEGER         ROAD    !  temporary road class code
        INTEGER         STID    !  tmp state ID

        LOGICAL      :: EFLAG = .FALSE.   !  true: error occurred
        LOGICAL      :: FFLAG = .FALSE.   !  true: fixed format
        LOGICAL      :: VFLAG = .FALSE.   !  true: header found

        CHARACTER          ARTP     !  tmp area code
        CHARACTER(4)       FCTP     !  tmp facility code
        CHARACTER(10)      FIPFMT   !  format to write co/st/cy to string
        CHARACTER(10)      RWTFMT   !  format to write rdway type to string
        CHARACTER(20)      FFORMAT  !  format description
        CHARACTER(300)     LINE     !  line buffer
        CHARACTER(300)     MESG     !  message buffer

        CHARACTER(LNKLEN3) CLNK     !  temporary link code
        CHARACTER(ALLLEN3) CSRCALL  !  buffer for source char, incl pol/act
        CHARACTER(FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(SCCLEN3) TSCC     !  temporary SCC
        CHARACTER(RWTLEN3) CRWT     !  roadway type no.
        CHARACTER(VIDLEN3) VIDZERO ! zero vehicle type

        CHARACTER(16) :: PROGNAME = 'RDVMIX' ! program name

C***********************************************************************
C   begin body of subroutine RDVMIX

C.............  Get environment variable values for this routine

        MESG = 'Indicator for VMT mix with fixed-column EMS-95 format'
        FFLAG = ENVYN( 'SMK_VMTMIX_FIXFMT', MESG, FFLAG, IOS )

        FFORMAT = 'list-directed'
        IF( FFLAG ) FFORMAT = 'fixed'

C.........  Set up formats
        WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3
        WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3

C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        VIDZERO = REPEAT( '0', VIDLEN3 )

C.........  Write status message
        MESG = 'Reading VMT mix file...'
        CALL M3MSG2( MESG )

C.........  Set up constants for loop...

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'VMT mix file' )

C.........  Set the number of table columns for list-directed format
        MXTCOL = NPRECOLL + NVTYPE

C.........  Allocate local memory
        ALLOCATE( VMIX( NVTYPE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VMIX', PROGNAME )
        ALLOCATE( SEGMENT( MXTCOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
        SEGMENT = ' '  ! array

C.........   First pass through file.  Count the number of non-blank lines
C.........   Also, scan for the header of vehicle types
        IREC = 0
        NREF = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading VMT mix file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

C.............  Look for header lines
            IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
                IF( LINE( 2:7 ) .EQ. 'VTYPES' ) THEN

C.....................  Determine number of vehicle types in file
                    VFLAG = .TRUE.
                    L2 = LEN_TRIM( LINE )
                    NV = GETNLIST( L2, LINE )
                    NV = NV - 1                ! subtract field for #VTYPES

C.....................  Allocate memory for vehicle type names and index
                    ALLOCATE( VTNAMES( NV ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'VTNAMES', PROGNAME )
                    ALLOCATE( VIDX( NV ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'VIDX', PROGNAME )
                    VTNAMES = ' '  ! array
                    VIDX = 0       ! array

C.....................  Store vehicle type names
                    CALL PARSLINE( LINE( 8:L2 ), NV, VTNAMES )

C.....................  Get position numbers for list of valid vehicle types
                    DO V = 1, NV

                        K = INDEX1( VTNAMES( V ), NVTYPE, CVTYPLST )
                        IF( K .LE. 0 ) THEN
                            L = LEN_TRIM( VTNAMES( V ) )
                            MESG = 'WARNING: Vehicle type "' // 
     &                             VTNAMES( V )( 1:L ) // '" in ' //
     &                             'vehicle mix file is not found in '//
     &                             'mobile codes file.'
                            CALL M3MSG2( MESG )

                        ELSE
                            VIDX( V ) = K

                        END IF

                    END DO

                    CYCLE   ! to head of read loop

                END IF

            END IF

            NREF = NREF + 1
         
        END DO         ! End first pass through file

C.........  Write error if header was not found
        IF( .NOT. VFLAG ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: #VTYPES header entry not found in ' //
     &             'vehicle mix file'
            CALL M3MSG2( MESG )

        ELSE IF( NV .LT. NVTYPE ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Only', NV, 'vehicle mix ' //
     &             'fractions are available in VMT mix file, but ' //
     &             CRLF() // BLANK10 // 'the mobile codes file ' //
     &             'indicates', NVTYPE, 'vehicle types are to be ' //
     &             'modeled.'
            CALL M3MSG2( MESG )

        END IF

        REWIND( FDEV )

C.........  Check for errors from reading either format
        IF( EFLAG ) THEN
            MESG = 'Problem reading vehicle mix file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Report format for XREF file
        ELSE
            L = LEN_TRIM( FFORMAT )
            MESG = 'NOTE: File read in as ' // FFORMAT( 1:L ) // 
     &             ' format.'
            CALL M3MSG2( MESG )

        END IF

C.........  Allocate memory for unsorted data used in all source categories
        ALLOCATE( VMTMIXA( NV, NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VMTMIXA', PROGNAME )
        ALLOCATE( INDXTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( CSCCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

C.........  In initialize before loop
        IREC = 0
        N    = 0
C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading VMT mix file at line',IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines and header lines
            IF( LINE .EQ. ' ' ) CYCLE
            IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

C.............  Separate line of data or read as fixed format...
C.............  Fixed format. NOTE - the hour field is ignored in this format,
C               and only the first hour of data will be read
            IF( FFLAG ) THEN

                STID = STR2INT( LINE( 1:2 ) )
                CYID = STR2INT( LINE( 3:5 ) )
                ARTP = ADJUSTL( LINE(  6:6  ) )
                FCTP = ADJUSTL( LINE(  7:10 ) )
                ROAD = STR2INT( ARTP // FCTP ) 
                FIP = 1000*STID + CYID

                WRITE( CFIP,FIPFMT ) FIP
                WRITE( CRWT,RWTFMT ) ROAD
                CLNK = ADJUSTL( LINE( 11:20 ) )

                P1   = FIXPREWD + 1
                P2   = FIXPREWD + FIXVMTWD

                DO J = 1, NV
                    V = VIDX( J )
                    IF( V .LE. 0 ) CYCLE

                    VMIX( V ) = STR2REAL( LINE( P1:P2 ) )
                    P1 = P1 + FIXVMTWD
                    P2 = P2 + FIXVMTWD

                    IF( VMIX( V ) .LT. 0.0 ) THEN
                        EFLAG = .TRUE.
                        K = NPRECOLF + J
                        CALL BAD_VMIX( IREC, K )
                        CYCLE
                    END IF

                END DO

C.............  List-directed format
            ELSE

                CALL PARSLINE( LINE, MXTCOL, SEGMENT )

C.................  Set source characteristics
                CFIP = SEGMENT( 1 )   ! country/state/county code
                CRWT = SEGMENT( 2 )   ! roadway code
                CLNK = SEGMENT( 3 )   ! link ID

C.................  Store VMT mix value and check for bad value
                DO J = 1, NV
                    K = NPRECOLL + J
                    V = VIDX( J )
                    IF( V .LE. 0 ) CYCLE

                    VMIX( V ) = STR2REAL( SEGMENT( K ) )

                    IF( VMIX( V ) .LT. 0.0 ) THEN
                        EFLAG = .TRUE.
                        CALL BAD_VMIX( IREC, K )
                        CYCLE
                    END IF

                END DO

            END IF

C.............  Replace any -9 fields with zeros
            CALL FLTRNEG( CFIP )
            CALL PADZERO( CFIP )
            CALL FLTRNEG( CRWT )
            CALL FLTRNEG( CLNK )

C.............  Increment table counter 
            N = N + 1
            IF( N .GT. NREF ) CYCLE  ! Ensure no overflow

C.............  Convert roadway type to internal SCC
            TSCC = CRWT // VIDZERO
            CALL PADZERO( TSCC )
    
C.............  Build source description
            CSRCALL = ' '
            CALL BLDCSRC( CFIP, RWTBLNK3, CLNK, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    POLBLNK3, CSRCALL )

C.............  Store unsorted table information
            INDXTA( N ) = N
            CSCCTA( N ) = TSCC
            CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC 

C.............  Use NVTYPE here, because of NVTYPE < NV, the unused fractions
C               will not be in the VMIX array.
            VMTMIXA( 1:NVTYPE, N ) = VMIX( 1:NVTYPE )

        END DO      ! End of loop on I for reading in temporal x-ref file

C.........  Set actual number of cross-reference entries
        NXREF = N

C.........  Check for errors reading cross-reference file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading VMT mix file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL M3MSG2( 'Processing VMT mix file...' )

C.........  Sort temporal cross-reference entries. Since CPOS was used in 
C           building CSRCTA, and CPOS will equal "0" when the x-ref entry is
C           not pol/act-specific, the non-pol/act-specific entries will
C           always appear first.  This is necessary for the table-generating
C           subroutines.

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

        CALL XREFTBL( 'VMTMIX', NXREF )

C.........  Deallocate local and global temporal memory needs
        DEALLOCATE( SEGMENT, VTNAMES, VIDX, INDXTA, CSRCTA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of VMT mix file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram writes out the error message for
C               invalid vehicle mixes
            SUBROUTINE BAD_VMIX( IREC, ICOL )

C.............  Subroutine arguments 
            INTEGER, INTENT (IN) :: IREC  ! record number
            INTEGER, INTENT (IN) :: ICOL  ! column number

C----------------------------------------------------------------------

            WRITE( MESG,94010 ) 'ERROR: Invalid vehicle mix at line', 
     &             IREC, 'column', ICOL, 'of VMT mix file.'
            CALL M3MESG( MESG )

            RETURN

C************   SUBPROGRAM FORMAT  STATEMENTS   *************************

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE BAD_VMIX

        END SUBROUTINE RDVMIX
