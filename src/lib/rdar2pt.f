
        SUBROUTINE RDAR2PT( FDEV, CDEV, LDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the area-to-point data in ASCII format and populates 
C     the MODAR2PT arrays.  Sets up the cross-references and 
C     calls the cross-reference assignments routine, which in
C     turn populates the by-source arrays that contains the
C     information needed for assigning the point source locations
C     to area sources.
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
C
C****************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables

C.........  This module is for cross reference tables
        USE M3UTILIO

        USE MODXREF, ONLY: INDXTA, CFIPTA, CSCCTA, CSRCTA, IARPTA

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVSCC, INVSCC, SCCDESC  ! Note - needed for reporting only

C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT, ONLY: MXROWA2P, NTABLA2P, NAR2PT, AR2PTABL,
     &                      NA2PSCC, A2PSCC, AR2PT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)   CRLF
C       LOGICAL        ENVYN
        LOGICAL        USEEXPGEO
C        EXTERNAL       CRLF, ENVYN, USEEXPGEO
        EXTERNAL     USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: FDEV   ! area-to-point factors unit no.
        INTEGER , INTENT (IN) :: CDEV   ! SCC descriptions unit no.
        INTEGER , INTENT (IN) :: LDEV   ! log file unit no.

C.............  Parameters
        INTEGER, PARAMETER :: NFIELDS = 10  ! no. input fields
        INTEGER, PARAMETER :: FBEG( NFIELDS ) = 
     &                      ( / 1 , 8, 11, 34, 39, 63,
     &                          86, 97, 114, 126 / )
        INTEGER, PARAMETER :: FEND( NFIELDS ) = 
     &                      ( / 6 , 9, 32, 37, 61, 84,
     &                          93, 109, 124, 139 / )

C...........   Local allocatable arrays...
C...........   For processing:
        INTEGER      , ALLOCATABLE :: IDXA2P  ( :,: )   ! sorting index
        INTEGER      , ALLOCATABLE :: NUMSCC  ( : )     ! no. SCCs per table
        TYPE( AR2PT ), ALLOCATABLE :: UNSRTA2P( :,: )   ! unsorted tables
        CHARACTER(FIPLEN3), ALLOCATABLE :: LOCFIP  ( : )   ! tmp FIPS codes
        CHARACTER(SCCLEN3), ALLOCATABLE :: AR2PTSCC( :,: ) ! SCCs per table

C...........   For reporting:
        INTEGER      , ALLOCATABLE :: SCCIDX  ( : )       ! sorting index
        INTEGER      , ALLOCATABLE :: SCCSECTN( : )       ! section no. for SCC
        CHARACTER(SCCLEN3), ALLOCATABLE :: INVSCCA( : ) ! unsorted SCCs


C...........   Local arrays
        CHARACTER(32) SEGMENT( NFIELDS )

C...........   Other local variables
        INTEGER         I, J, K, L, N     ! counters and indices

        INTEGER      :: CNT = 0        ! tmp record counter 
        INTEGER      :: FIPCNT = 0     ! counter for no. entries per SCC/FIPs
        INTEGER         IOS            ! i/o status
        INTEGER      :: MXSCC = 0      ! max no. of SCCs per table 
        INTEGER      :: NXREF = 0      ! number of cross-ref entries 

        REAL         :: SUMTEST = 0.   ! value to check that factors sum to 1.

        LOGICAL      :: EFLAG  = .FALSE.  ! true: error found
        LOGICAL      :: WFLAG             ! true: convert lat-lons to Western hemisphere

        CHARACTER(256)        MESG    !  message buffer
        CHARACTER(512)     :: LINE    !  input line
        CHARACTER(FIPLEN3) :: CFIP    !  tmp char co/st/cy code
        CHARACTER(FIPLEN3) :: LFIP    !  previous co/st/cy FIPS code 
        CHARACTER(SCCLEN3) :: TSCC    !  tmp SCC code

        CHARACTER(16) :: PROGNAME = 'RDAR2PT' ! program name

C***********************************************************************
C   begin body of subroutine RDAR2PT

C.........  Rewind area-to-point file
        IF ( FDEV .NE. 0 ) THEN
            REWIND( FDEV )
        ELSE
            MESG = 'Area-to-point factors file is not opened!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Loop through the input file to determine how many header
C           lines there are and the maximum number of entries in 
C           a single table
        CALL READ_ARTOPT_FILE( 'COUNT', MXSCC )

C.........  Allocate memory for unsorted input tables
        ALLOCATE( IDXA2P( MXROWA2P, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXA2P', PROGNAME )
        ALLOCATE( UNSRTA2P( MXROWA2P, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNSRTA2P', PROGNAME )
        ALLOCATE( AR2PTSCC( MXSCC, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTSCC', PROGNAME )
        ALLOCATE( NUMSCC( NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NUMSCC', PROGNAME )
        IDXA2P = 0                ! array
        UNSRTA2P%FIP   = ' '      ! array
        UNSRTA2P%LAT   = BADVAL3  ! array
        UNSRTA2P%LON   = BADVAL3  ! array
        UNSRTA2P%ALLOC = 1.       ! array
        UNSRTA2P%NAME  = ' '      ! array
        AR2PTSCC = ' '            ! array
        NUMSCC = 0                ! array

C.........  Allocate memory for sorted input tables
        ALLOCATE( LOCFIP( MXROWA2P ), STAT=IOS )  ! local array
        CALL CHECKMEM( IOS, 'LOCFIP', PROGNAME )
        ALLOCATE( NAR2PT( NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NAR2PT', PROGNAME )
        ALLOCATE( AR2PTABL( MXROWA2P, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTABL', PROGNAME )
        NAR2PT         = 0        ! array
        AR2PTABL%FIP   = ' '      ! array
        AR2PTABL%LAT   = BADVAL3  ! array
        AR2PTABL%LON   = BADVAL3  ! array
        AR2PTABL%ALLOC = 1.       ! array
        AR2PTABL%NAME  = ' '      ! array

C.........  Check if lat-lons should be converted to western hemisphere
        MESG = 'Western hemisphere flag'
        WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .FALSE., IOS )

C.........  Read and store contents of the file (both the SCC 
C           entries as well as the tables).
C.........  In this call, populate INXA2P, UNSRTA2P, AR2PTSCC,
C           NUMSCC, and NAR2PT
        CALL READ_ARTOPT_FILE( 'STORE', MXSCC )

C.........  Sort index for finding sorted order for each table.
C.........  Compute expected number of cross-reference entries
        NXREF = 0
        DO N = 1, NTABLA2P

C.............  Transfer FIPS codes to local array
            LOCFIP = UNSRTA2P(:,N)%FIP    ! array

C.............  Sort for current table
            CALL SORTIC( NAR2PT(N), IDXA2P(1,N), LOCFIP )

C.............  Compute maximum expected x-ref entries
            LFIP = ' '
            DO I = 1, NAR2PT( N )
                J = IDXA2P( I,N )
                CFIP = UNSRTA2P( J,N )%FIP
                IF( LFIP .NE. CFIP ) THEN
                    NXREF = NXREF + NUMSCC( N )
                END IF
                LFIP = CFIP
            END DO          ! end loop through rows in table

        END DO              ! end loop through tables

C........  Allocate memory for cross-reference arrays
        ALLOCATE( INDXTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( CFIPTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFIPTA', PROGNAME )
        ALLOCATE( CSCCTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        ALLOCATE( IARPTA( NXREF,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IARPTA', PROGNAME )
        INDXTA = 0   ! array
        CFIPTA = ' ' ! array
        CSCCTA = ' ' ! array
        CSRCTA = ' ' ! array
        IARPTA = 0   ! array

C.........  Store sorted input tables.
C.........  Build arrays for giving to cross-reference routines for
C           grouping the cross-references.  Although, for now, there
C           will be only the FIPS//SCC group.
        CNT  = 0
        DO N = 1, NTABLA2P

C...........  Store sorted input tables
            DO I = 1, NAR2PT( N ) 
                J = IDXA2P( I,N )
                AR2PTABL( I,N ) = UNSRTA2P( J,N )
            END DO
            
C...........  Loop through SCCs for this section of the file
C             and create cross-referencing arrays
            DO K = 1, NUMSCC( N )

C...............  Reset number of entries per FIPS code
                LFIP = ' '

C...............  Loop through rows for this section of the file
                DO I = 1, NAR2PT( N )

                    CFIP = AR2PTABL( I,N )%FIP
                    TSCC = AR2PTSCC( K,N )

C..................  If this row is a new FIPS code
                    IF( CFIP .NE. LFIP ) THEN
                        CNT = CNT + 1  ! count of cross-reference entries

C.......................  If count of entries is w/i dimensioned array
                        IF( CNT .LE. NXREF ) THEN
                            INDXTA( CNT )   = CNT
                            CFIPTA( CNT )   = CFIP
                            CSCCTA( CNT )   = TSCC
                            IARPTA( CNT,1 ) = N
                            IARPTA( CNT,2 ) = I
                            CSRCTA( CNT )   = CFIP // TSCC
                        END IF

C.......................  Reset number of entries per FIPS code
                        FIPCNT = 0

                    END IF

C....................  Increment the number of entries per FIPS code
                    FIPCNT = FIPCNT + 1

C....................  Store the number of entries for this FIPS/SCC
                    IF( CNT .LE. NXREF ) IARPTA( CNT,3 ) = FIPCNT

C....................  Store previous FIPS code for next iteration
                    LFIP = CFIP

                END DO   ! end loop through rows of current table
            END DO       ! end loop through SCCs of current table

        END DO           ! end loop through tables

C.........  Check if count exceeded maximum expected
        IF( CNT .GT. NXREF ) THEN
            WRITE( MESG,94010 ) 
        ELSE
            NXREF = CNT
        END IF

C.........  Sort area-to-point factors cross-reference
        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Check that fractions for each FIPS/SCC combo sums to 1.
        DO I = 1, NXREF

            J    = INDXTA( I )
            CFIP = CFIPTA( J )
            TSCC = CSCCTA( J ) 

C............  Loop through entries for current FIPS/SCC and compute
C              sum of allocation factors
            N   = IARPTA( J,1 )
            CNT = IARPTA( J,2 ) - 1
            SUMTEST = 0.
            DO K = 1, IARPTA( J,3 )
                CNT = CNT + 1
                SUMTEST = SUMTEST + AR2PTABL( CNT,N )%ALLOC
            END DO

C............  If sum of allocation factors is outside allowable 
C              range, then write error
            IF( SUMTEST .GT. 1.001 .OR.
     &          SUMTEST .LT. 0.999      ) THEN
    
                EFLAG = .TRUE.
                WRITE( MESG, 94020 ) 'ERROR: Sum of factors for '//
     &                 'Co/St/Cy= '// CFIP//'and SCC= '// TSCC//
     &                 CRLF()// BLANK10 // 'is not equal to 1. '//
     &                 'Value is ', SUMTEST, '.'
                CALL M3MSG2( MESG )

            END IF

        END DO

C.........  Abort if errors found
        IF( EFLAG ) THEN
            MESG = 'Area-to-point factors file has bad data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Call cross-referencing routine, which will also populate
C           the ARPT09 array (full FIPs and full SCC matches only)
        CALL XREFTBL( 'AR2PT', NXREF )

C.........  Report to the log file the SCCs, SCC descriptions, and
C           the section of the input file set for each SCC.
C.........  To create this report, we will create a fake INVSCC array
C           for the MODLISTS module, which will allow us to read
C           and assign the SCC descriptions for this report...

C.........  Sum SCC count
        NINVSCC = 0
        DO N = 1, NTABLA2P
            NINVSCC = NINVSCC + NUMSCC( N )
        END DO

        NA2PSCC = NINVSCC
        
C.........  Allocate temporary "inventory" SCC list from MODLISTS
        ALLOCATE( INVSCC( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSCC', PROGNAME )
        ALLOCATE( A2PSCC( NA2PSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'A2PSCC', PROGNAME )

C.........  Allocate unsorted arrays
        ALLOCATE( SCCIDX( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCIDX', PROGNAME )
        ALLOCATE( INVSCCA( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSCCA', PROGNAME )
        ALLOCATE( SCCSECTN( NINVSCC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCCSECTN', PROGNAME )

C.........  Store unsorted SCCs
        J = 0
        DO N = 1, NTABLA2P
            DO K = 1, NUMSCC( N )
                J = J + 1
                SCCIDX  ( J ) = J
                INVSCCA ( J ) = AR2PTSCC( K,N )
                SCCSECTN( J ) = N
            END DO
        END DO

C.........  Sort SCC list, needed for reading SCC descriptions
        CALL SORTIC( NINVSCC, SCCIDX, INVSCCA )
        DO J = 1, NINVSCC
            INVSCC( J ) = INVSCCA( SCCIDX( J ) )
        END DO
        
        A2PSCC = INVSCC

C.........  Retrieve SCC descriptions
        CALL RDSCCDSC( CDEV )

C.........  Write header of report
        MESG = 'NOTE: The area-to-point factors file ' //
     &         'included the following SCCs and '// CRLF()//
     &         BLANK10 // 'section numbers'
        CALL M3MESG( MESG )
        WRITE( LDEV, 94010 ) ' '

        MESG = 'SCC                  Section   SCC Description'
        WRITE( LDEV, 94010 ) BLANK10 // TRIM( MESG )
        WRITE( LDEV, 94010 ) BLANK10 // REPEAT( '-', 70 )

C.........  Loop through SCCs and list them, section no.'s, and descriptions
        DO J = 1, NINVSCC
            K = SCCIDX( J )
            WRITE( LDEV,94675 ) INVSCC( J ), SCCSECTN( K ), 
     &                          TRIM( SCCDESC( J ) )
        END DO
        WRITE( LDEV, 94010 ) ' '

C.........  Deallocate borrowed MODLISTS arrays
        DEALLOCATE( INVSCC, SCCDESC )

C.........  Deallocate local memory
        DEALLOCATE( IDXA2P, NUMSCC, UNSRTA2P, AR2PTSCC )
        DEALLOCATE( SCCIDX, INVSCCA, SCCSECTN )

C.........  Deallocate cross-reference sorting arrays
        DEALLOCATE( INDXTA, CFIPTA, CSCCTA, CSRCTA, IARPTA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT(  A, 1X, F10.7, 1X, A )

94675   FORMAT( 10X, A20, 4X, I2.2, 5X, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram counts entries in or reads
C               the contents of the area-to-point factors file
            SUBROUTINE READ_ARTOPT_FILE( STATUS, MXSCC )

C...............   EXTERNAL FUNCTIONS and their descriptions:
            LOGICAL         CHKINT
            LOGICAL         CHKREAL
            LOGICAL         BLKORCMT
C           CHARACTER(2)    CRLF
            INTEGER         GETFLINE
            INTEGER         GETNLIST
C           INTEGER         STR2INT
C           REAL            STR2REAL

C            EXTERNAL CHKINT, CHKREAL, CRLF, GETFLINE, GETNLIST,
C     &               STR2INT, STR2REAL, BLKORCMT
        EXTERNAL     CHKINT, CHKREAL, GETFLINE, GETNLIST, BLKORCMT

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: STATUS  ! call status: COUNT|STORE
            INTEGER, INTENT (IN OUT)  :: MXSCC   ! count the max SCC per section

C.............  Local variables
            INTEGER          I, L, L1, L2, K, N, NS   ! counters and indices

            INTEGER          CNT        ! table record counter
            INTEGER          CNY        ! tmp county code
            INTEGER          COU        ! tmp country code
            INTEGER          IOS        ! i/o status
            INTEGER          IREC       ! record number
            INTEGER       :: LINSCC = 0 ! number of SCCs on current header line
            INTEGER, SAVE :: NLINES     ! no. lines in input file
            INTEGER       :: NSCC = 0   ! SCC counter
            INTEGER          NTBL       ! no. tables
            INTEGER          STA        ! tmp state code

            REAL             ALLOC      ! tmp allocation factor
            REAL             LAT        ! tmp latitude
            REAL             LON        ! tmp longitude
            
            CHARACTER(FIPLEN3) CFIP     ! tmp co/st/cy FIPS code
            CHARACTER(SCCLEN3) SEGSCC   ! segment SCC

            LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
            LOGICAL       :: PREVPKT  = .FALSE. ! true: previous line was a packet
            LOGICAL       :: THISPKT  = .FALSE. ! true: this line is a packet

C......................................................................

C.............  If the first time the internal subprogram is called
            IF( FIRSTIME ) THEN

C................  Determine the number of lines in the input file
                NLINES = GETFLINE( FDEV, "ARTOPNT file" )

                FIRSTIME = .FALSE.

            END IF

C............  Loop through lines of input file and
            IREC    = 0
            NTBL    = 0
            PREVPKT = .FALSE.
            DO I = 1, NLINES

                READ( FDEV, 93000, END=2002, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                 'ERROR: System I/O error', IOS, 'reading ' // 
     &                 'area-to-point factors file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Skip comment and blank lines
                IF( BLKORCMT ( LINE ) ) CYCLE

C................  If line is a header line, add number of SCCs for
C                  the current table
                THISPKT = .FALSE.
                L = INDEX( LINE, '/LOCATIONS/' )
                IF( L .GT. 0 ) THEN

C...................  If previous line was not another packet
                    IF( .NOT. PREVPKT ) THEN

C.......................  Update the number of tables
                        NTBL = NTBL + 1

C.......................  Initialize the SCC count
                        L1 = 12               ! based on LOCATIONS packet
                        L2 = LEN_TRIM( LINE )
                        LINSCC = GETNLIST( L2 - L1 + 1, LINE( L1:L2 ) )
                        NSCC = LINSCC
                        NS   = 0

C.......................  Initialise count of entries per table
                        CNT = 0
                        
C...................  Otherwise add to the SCC count
                    ELSE
                        L1 = 12               ! based on LOCATIONS packet
                        L2 = LEN_TRIM( LINE )
                        LINSCC = GETNLIST( L2 - L1 + 1, LINE( L1:L2 ) )
                        NSCC = NSCC + LINSCC

                    END IF

C....................  Store maximum SCC value
                    MXSCC = MAX( MXSCC, NSCC )

C....................  Give error if no SCCs are provided in the packet
                    IF( LINSCC .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: No SCCs provided '//
     &                         'in /LOCATIONS/ packet at line', IREC
                        CALL M3MSG2( MESG )
                    END IF

C....................  Set packet status 
                    THISPKT = .TRUE.

                END IF

C...............  Make sure that header has appeared at least once
                IF( NTBL .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: No header line found '//
     &                     'before line', IREC
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

C...............  Increment table counter and keep track of maximum
                IF( .NOT. THISPKT ) THEN
                    CNT = CNT + 1
                    MXROWA2P = MAX( MXROWA2P, CNT )
                END IF

C................  End read if just counting
                IF( STATUS .EQ. 'STORE' ) THEN

C....................  If header line, store SCCs and count
                    IF( THISPKT ) THEN
                        NUMSCC( NTBL ) = NSCC
                        
                        CALL PARSLINE( LINE( L1:L2 ), LINSCC, SEGMENT )

                        DO N = 1, LINSCC
                            NS = NS + 1

                            L = LEN_TRIM( SEGMENT( N ) )
                            IF( L .GT. SCCLEN3 ) THEN
                                EFLAG = .TRUE.
                                WRITE( MESG,94010 ) 'ERROR: SCC "'//
     &                            SEGMENT( 1:L ) // '" at line', IREC,
     &                            'has', L, 'characters, which '//
     &                            'exceeds the maximum of', SCCLEN3
                                CALL M3MSG2( MESG )
                                CYCLE
                            END IF

                            SEGSCC = TRIM( SEGMENT( N ) )
                            CALL PADZERO( SEGSCC )
                            AR2PTSCC( NS,NTBL ) = SEGSCC

                        END DO

C....................  Otherwise, parse and store table
                    ELSE

                        DO N = 1, NFIELDS
                            SEGMENT(N)= ADJUSTL( LINE(FBEG(N):FEND(N)) )
                        END DO

C........................  Check that FIPS code is an integer
                        IF( .NOT. USEEXPGEO() .AND. 
     &                      .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Country/state'//
     &                             '/county code is not an integer '//
     &                             'at line', IREC
                            CALL M3MSG2( MESG )
                        END IF

C........................  Check that reals are reals
                        IF( .NOT. CHKREAL( SEGMENT( 8 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Longitude is '//
     &                        'not a floating point value at line', IREC
                            CALL M3MSG2( MESG )
                        END IF
                        IF( .NOT. CHKREAL( SEGMENT( 9 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Latitude is '//
     &                        'not a floating point value at line', IREC
                            CALL M3MSG2( MESG )
                        END IF
                        IF( .NOT. CHKREAL( SEGMENT( 10 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Allocation '//
     &                        'factor is not a floating point value '//
     &                        'at line', IREC
                            CALL M3MSG2( MESG )
                        END IF

                        IF( EFLAG ) CYCLE

C........................  Store total count of packet for this table
                        NAR2PT( NTBL ) = CNT

C........................  Store contents of this row
                        CFIP  = SEGMENT( 1 )
                        CALL PADZERO( CFIP )
                        LON   = STR2REAL( SEGMENT( 8 ) )
                        IF( WFLAG .AND. LON > 0 ) LON = -LON  ! convert to western hemisphere
                        LAT   = STR2REAL( SEGMENT( 9 ) )
                        ALLOC = STR2REAL( SEGMENT( 10 ) )

                        IDXA2P        ( CNT,NTBL ) = CNT
                        UNSRTA2P( CNT,NTBL )%FIP   = CFIP
                        UNSRTA2P( CNT,NTBL )%LAT   = LAT
                        UNSRTA2P( CNT,NTBL )%LON   = LON
                        UNSRTA2P( CNT,NTBL )%ALLOC = ALLOC
                        UNSRTA2P( CNT,NTBL )%NAME  = TRIM( SEGMENT(5) )
                    END IF
                END IF

C...............  Store packet status for next iteration
                PREVPKT = THISPKT

            END DO   ! end of loop through lines

            NTABLA2P = NTBL

C.............  Abort if errors found
            IF( EFLAG ) THEN
                MESG = 'Problem reading area-to-point factors file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Rewind file for next call
            REWIND( FDEV )

            RETURN

C.............  Abort with errors...
2002        WRITE( MESG,94010 ) 'ERROR: Unexpected end of area-to-'//
     &             'point factors file at line', IREC
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C************   SUBPROGRAM FORMAT  STATEMENTS   *************************

C...........   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_ARTOPT_FILE

C----------------------------------------------------------------------

        END SUBROUTINE RDAR2PT
