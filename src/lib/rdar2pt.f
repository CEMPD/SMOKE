
        SUBROUTINE RDAR2PT( FDEV )

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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables

C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2    CRLF
        EXTERNAL       CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: FDEV   ! iventory table unit no.

C.............  Parameters
        INTEGER, PARAMETER :: NFIELDS = 11  ! no. input fields
        INTEGER, PARAMETER :: FBEG( NFIELDS ) = 
     &                      ( / 1 , 7, 10, 22, 32, 52,
     &                          66, 78, 92, 108, 122 / )
        INTEGER, PARAMETER :: FEND( NFIELDS ) = 
     &                      ( / 6 , 8, 20, 30, 50, 64,
     &                          76, 90, 105, 120, 127 / )

C...........   Local allocatable arrays
        INTEGER      , ALLOCATABLE :: IDXA2P  ( :,: )   ! sorting index
        INTEGER      , ALLOCATABLE :: LOCFIP  ( : )     ! tmp FIPS codes
        INTEGER      , ALLOCATABLE :: NUMSCC  ( : )     ! no. SCCs per table
        TYPE( AR2PT ), ALLOCATABLE :: UNSRTA2P( :,: )   ! unsorted tables
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE :: AR2PTSCC( :,: ) ! SCCs per table

C...........   Local arrays
        CHARACTER*32 SEGMENT( NFIELDS )

C...........   Other local variables
        INTEGER         I, J, K, L, N     ! counters and indices

        INTEGER      :: CNT            ! tmp record counter 
        INTEGER      :: FIP            ! tmp co/st/cy FIPS code 
        INTEGER      :: FIPCNT         ! counter for no. entries per SCC/FIPs
        INTEGER         IOS            ! i/o status
        INTEGER      :: LFIP           ! previous co/st/cy FIPS code 
        INTEGER      :: MXSCC          ! max no. of SCCs per table 
        INTEGER      :: NXREF          ! number of cross-ref entries 

        LOGICAL      :: EFLAG  = .FALSE.  ! true: error found

        CHARACTER*10              FIPFMT  !  format to write co/st/cy to string
        CHARACTER*256             MESG    !  message buffer
        CHARACTER*512          :: LINE    !  input line
        CHARACTER(LEN=FIPLEN3) :: CFIP    !  tmp char co/st/cy code
        CHARACTER(LEN=SCCLEN3) :: TSCC    !  tmp SCC code

        CHARACTER*16 :: PROGNAME = 'RDAR2PT' ! program name

C***********************************************************************
C   begin body of subroutine RDAR2PT

C.........  Rewind area-to-point file
        IF ( FDEV .NE. 0 ) THEN
            REWIND( FDEV )
        ELSE
            MESG = 'Area-to-point factors file is not opened!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set up formats
        WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

C.........  Loop through the input file to determine how many header
C           lines there are and the maximum number of entries in 
C           a single table
        CALL READ_ARTOPT_FILE( 'COUNT' )

C.........  Allocate memory for unsorted input tables
        ALLOCATE( IDXA2P( MXROWA2P, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXA2P', PROGNAME )
        ALLOCATE( UNSRTA2P( MXROWA2P, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNSRTA2P', PROGNAME )
        ALLOCATE( AR2PTSCC( MXSCC, NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTSCC', PROGNAME )
        ALLOCATE( NUMSCC( NTABLA2P ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTSCC', PROGNAME )
        IDXA2P = 0                ! array
        UNSRTA2P%FIP   = 0        ! array
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
        AR2PTABL%FIP   = 0        ! array
        AR2PTABL%LAT   = BADVAL3  ! array
        AR2PTABL%LON   = BADVAL3  ! array
        AR2PTABL%ALLOC = 1.       ! array
        AR2PTABL%NAME  = ' '      ! array

C.........  Read and store contents of the file (both the SCC 
C           entries as well as the tables).
C.........  In this call, populate INXA2P, UNSRTA2P, AR2PTSCC,
C           NUMSCC, and NAR2PT
        CALL READ_ARTOPT_FILE( 'STORE' )

C.........  Sort index for finding sorted order for each table.
C.........  Compute expected number of cross-reference entries
        NXREF = 0
        DO N = 1, NTABLA2P

C.............  Transfer FIPS codes to local array
            LOCFIP = UNSRTA2P(:,N)%FIP    ! array

C.............  Sort for current table
            CALL SORTI1( NAR2PT(N), IDXA2P(1,N), LOCFIP )

C.............  Compute maximum expected x-ref entries
            LFIP = -9
            DO I = 1, NAR2PT( N )
                J = IDXA2P( I,N )
                FIP = UNSRTA2P( J,N )%FIP
                IF( LFIP .NE. FIP ) THEN
                    NXREF = NXREF + NUMSCC( N )
                END IF
                LFIP = FIP
            END DO          ! end loop through rows in table

        END DO              ! end loop through tables

C........  Allocate memory for cross-reference arrays
        ALLOCATE( INDXTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( IFIPTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIPTA', PROGNAME )
        ALLOCATE( CSCCTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NXREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
        ALLOCATE( IARPTA( NXREF,3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IARPTA', PROGNAME )
        INDXTA = 0   ! array
        IFIPTA = 0   ! array
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
                LFIP = -9

C...............  Loop through rows for this section of the file
                DO I = 1, NAR2PT( N )

                    FIP  = AR2PTABL( I,N )%FIP
                    TSCC = AR2PTSCC( K,N )

C..................  If this row is a new FIPS code
                    IF( FIP .NE. LFIP ) THEN
                        CNT = CNT + 1  ! count of cross-reference entries

C.......................  If count of entries is w/i dimensioned array
                        IF( CNT .LE. NXREF ) THEN
                            INDXTA( CNT )   = CNT
                            IFIPTA( CNT )   = FIP
                            CSCCTA( CNT )   = TSCC
                            IARPTA( CNT,1 ) = N
                            IARPTA( CNT,2 ) = I

                            WRITE( CFIP,FIPFMT ) FIP
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
                    LFIP = FIP

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

C.........  Call cross-referencing routine, which will also populate
C           the ARPT09 array (full FIPs and full SCC matches only)
        CALL XREFTBL( 'AR2PT', NXREF )

C.........  Deallocate local memory
        DEALLOCATE( IDXA2P, NUMSCC, UNSRTA2P, AR2PTSCC )

C.........  Deallocate cross-reference sorting arrays
        DEALLOCATE( INDXTA, IFIPTA, CSCCTA, CSRCTA, IARPTA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram counts entries in or reads
C               the contents of the area-to-point factors file
            SUBROUTINE READ_ARTOPT_FILE( STATUS )

C...............   EXTERNAL FUNCTIONS and their descriptions:
            LOGICAL         CHKINT
            LOGICAL         CHKREAL
            CHARACTER*2     CRLF
            INTEGER         GETFLINE
            INTEGER         GETNLIST
            INTEGER         STR2INT
            REAL            STR2REAL

            EXTERNAL CHKINT, CHKREAL, CRLF, GETFLINE, GETNLIST,
     &               STR2INT, STR2REAL

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: STATUS  ! call status: COUNT|STORE

C.............  Local variables
            INTEGER          I, L, L1, L2, K, N, NS   ! counters and indices

            INTEGER          CNT        ! table record counter
            INTEGER          CNY        ! tmp county code
            INTEGER          COU        ! tmp country code
            INTEGER          FIP        ! tmp co/st/cy FIPS code
            INTEGER          IOS        ! i/o status
            INTEGER          IREC       ! record number
            INTEGER          LINSCC     ! number of SCCs on current header line
            INTEGER, SAVE :: NLINES     ! no. lines in input file
            INTEGER          NSCC       ! SCC counter
            INTEGER          NTBL       ! no. tables
            INTEGER          STA        ! tmp state code

            REAL             ALLOC      ! tmp allocation factor
            REAL             LAT        ! tmp latitude
            REAL             LON        ! tmp longitude

            LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
            LOGICAL       :: PREVPKT  = .FALSE. ! true: previous line was a packet
            LOGICAL       :: THISPKT  = .FALSE. ! true: this line is a packet

C......................................................................

C.............  If the first time the internal subprogram is called
            IF( FIRSTIME ) THEN

C................  Determine the number of lines in the input file
                NLINES = GETFLINE( FDEV )

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

C.................  Skip comment lines
                IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

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

                            AR2PTSCC( NS,NTBL ) = TRIM( SEGMENT( N ) )

                        END DO

C....................  Otherwise, parse and store table
                    ELSE

                        DO N = 1, NFIELDS
                            SEGMENT(N)= ADJUSTL( LINE(FBEG(N):FEND(N)) )
                        END DO

C........................  Check that FIPS code is an integer
                        IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
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
                        FIP   = STR2INT ( SEGMENT( 1 ) )
                        LON   = STR2REAL( SEGMENT( 8 ) )
                        LAT   = STR2REAL( SEGMENT( 9 ) )
                        ALLOC = STR2REAL( SEGMENT( 10 ) )

                        IDXA2P        ( CNT,NTBL ) = CNT
                        UNSRTA2P( CNT,NTBL )%FIP   = FIP
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
