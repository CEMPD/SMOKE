
       SUBROUTINE RDIDAMB ( FDEV, NRAWIN, NRAWBV, MXIDAT, WKSET, 
     &                      INVDNAM, NRAWOUT, EFLAG, NDROP, VDROP )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       This subroutine reads the IDA-formatted mobile files.  It can
C       be called multiple times for multiple files.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C...........   Modules for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module is for mobile-specific data
        USE MODMOBIL

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
	
C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        INTEGER         ENVINT
        INTEGER         INDEX1
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         STR2INT
        REAL            STR2REAL
        CHARACTER*2     CRLF

        EXTERNAL        CHKINT, CHKREAL, ENVINT, CRLF, INDEX1, FIND1, 
     &                  FINDC, PADZERO, STR2INT, STR2REAL

C...........   Subroutine arguments. Note that number and amount of dropped
C              VMT is initialied in calling program.

        INTEGER     , INTENT (IN) :: FDEV        ! file unit number
        INTEGER     , INTENT (IN) :: NRAWIN      ! total raw record count
        INTEGER     , INTENT (IN) :: NRAWBV      ! total rec x data var count
        INTEGER     , INTENT (IN) :: MXIDAT      ! max no of inventory pols/act
        INTEGER     , INTENT (IN) :: WKSET       ! weekly profile interpretation
        CHARACTER(*), INTENT (IN) :: INVDNAM( MXIDAT ) ! inv pol/actvty names
        INTEGER     , INTENT(OUT) :: NRAWOUT     ! valid raw record count
        INTEGER     , INTENT(OUT) :: NDROP       ! number of records dropped
        REAL        , INTENT(OUT) :: VDROP( MXIDAT ) ! sum of VMT dropped
        LOGICAL     , INTENT(OUT) :: EFLAG       ! error flag

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: MXVARFIL = 113 ! maximum data variables in file
        INTEGER, PARAMETER :: MBOTWIDE = 26  ! total width of all data fields
        INTEGER, PARAMETER :: MBNONDWD = 22  ! width of non-data fields

C...........   Local parameters, dependent
        INTEGER, PARAMETER :: LINSIZ  = MBNONDWD + MXVARFIL * MBOTWIDE

C...........   Local parameter arrays...
C...........   Start and end positions in the file format of the first set
C              of pollutant fields.
        INTEGER  :: ISINIT( NMBPPOL3 ) = ( / 23,36 / )

        INTEGER  :: IEINIT( NMBPPOL3 ) = ( / 35,48 / )

C...........   Local arrays
        INTEGER          IS( NMBPPOL3 )  ! start position for each pol char
        INTEGER          IE( NMBPPOL3 )  ! end position for each pol char

C...........   Counters of total number of input records
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER, SAVE :: NSRCVAR = 0 ! cumulative source x pollutants count

C...........   Local variables

        INTEGER         I, J, K, L, V     ! indices

        INTEGER         CYID           ! tmp county code
        INTEGER         ES             ! valid raw record counter
        INTEGER         FIP            ! tmp country/state/county
        INTEGER         ICC            ! tmp country code
        INTEGER         IDROP          ! local dropped raw record counter
        INTEGER         INY            ! inventory year
        INTEGER         IOS            ! I/O status
        INTEGER         IREC           ! record counter
        INTEGER         ISPEED         ! tmp speed, integer
        INTEGER         IVT            ! tmp vehicle type code
        INTEGER, SAVE:: MXWARN         ! maximum number of warnings
        INTEGER         NVAR           ! number of data variables
        INTEGER, SAVE:: NWARN =0       ! number of warnings in this routine
        INTEGER         SS             ! counter for sources
        INTEGER         RWT            ! roadway type
        INTEGER         SCCLEN         ! length of SCC string 
        INTEGER         STID           ! tmp state code
        INTEGER         TPF            ! tmp temporal adjustments setting

        REAL            DANN      ! tmp annual-ave data value
        REAL         :: DOZN = 0. ! tmp ozone-season-ave data value
        REAL            VMTR      ! tmp vehicle miles traveled (10^6 VMT/yr)

        LOGICAL    :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL       INVALID          ! tmp error flag for current record
        LOGICAL    :: SFLAG  = .FALSE. ! true: blank or zero speeds  found
        LOGICAL    :: OLDFMT = .FALSE. ! true: old mobile IDA format used

        CHARACTER*2   SPEEDCHR       ! tmp speed, character
        CHARACTER*20  VIDFMT         ! vehicle type ID format
        CHARACTER*20  RWTFMT         ! roadway type number format
        CHARACTER*300 MESG           ! message buffer

        CHARACTER(LEN=LINSIZ)  LINE  ! read buffer for a line
        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIPS code
        CHARACTER(LEN=VIDLEN3) CIVT  ! tmp vehicle type ID
        CHARACTER(LEN=LNKLEN3) CLNK  ! tmp link ID 
        CHARACTER(LEN=RWTLEN3) CRWT  ! tmp roadway type
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp source classification code
        CHARACTER(LEN=VTPLEN3) VTYPE ! tmp vehicle type

        CHARACTER*16 :: PROGNAME = 'RDIDAMB'   ! program name

C***********************************************************************
C   Begin body of subroutine RDIDAMB

C.........  Firstime routine is called, get the number of warnings
        IF( FIRSTIME ) THEN
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
            FIRSTIME = .FALSE.
        END IF

C.........  Reinitialize for multiple subroutine calls
        ICC   = -9
	IDROP = 0
        INY   = 0
        NVAR  = 0

C.........  Create formats
        WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3
        WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3

C.........  Make sure the file is at the beginning
        REWIND( FDEV )

C........................................................................
C.............  Head of the main read loop  .............................
C........................................................................

        SS   = NSRCSAV
        ES   = NSRCVAR
        IREC = 0
        TPF  = MTPRFAC * WKSET
        CLNK = ' '
        DO

            INVALID = .FALSE.

            READ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE
            IREC = IREC + 1
 
            IF ( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading inventory '//
     &              'file at line', IREC
                CALL M3MESG( MESG )

            END IF
	
            L = LEN_TRIM( LINE )  ! store width of line and check

C.............  Skip blank lines
            IF( L .EQ. 0 ) CYCLE

C.............  Scan for header lines and check to ensure all are set 
C               properly
            CALL GETHDR( MXVARFIL, MXIDAT, .TRUE., .TRUE., .FALSE., 
     &                   INVDNAM, LINE, ICC, INY, NVAR, IOS )

C.............  Interpret error status
            IF( IOS .EQ. 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXVARFIL=', MXVARFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  If NVAR was not set, then we are using the old format
C.............  Set the header parameters that we'd be setting otherwise
            IF( NVAR .EQ. 0 ) THEN

                NVAR = 1
                ALLOCATE( TMPNAM( NVAR ), STAT=IOS )
	        CALL CHECKMEM( IOS, 'TMPNAM', PROGNAME )
                ALLOCATE( POLPOS( NVAR ), STAT=IOS )
	        CALL CHECKMEM( IOS, 'POLPOS', PROGNAME )
                ALLOCATE( INVDUNT( 1,NMBPPOL3 ), STAT=IOS )
	        CALL CHECKMEM( IOS, 'INVDUNT', PROGNAME )
                ALLOCATE( INVDCNV( 1,NMBPPOL3 ), STAT=IOS )
	        CALL CHECKMEM( IOS, 'INVDCNV', PROGNAME )
                TMPNAM( 1 ) = 'VMT'                
                POLPOS( 1 ) = 1
                INVDUNT( 1,1 ) = 'mi/yr'
                INVDUNT( 1,2 ) = 'mi/day'
                INVDCNV      = 1000000     ! array
                ISINIT( 1 ) = 16
                IEINIT( 1 ) = 28
                OLDFMT = .TRUE.

            END IF

C.............  Check state/county codes, error for missing
            IF( .NOT. CHKINT( LINE( 1:2 ) ) .OR. 
     &          .NOT. CHKINT( LINE( 3:5 ) ) .OR.
     &          LINE( 1:2 ) .EQ. ' '        .OR.
     &          LINE( 3:5 ) .EQ. ' '             ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: State and/or county ' //
     &                 'code is non-integer or missing at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Initialize start and end positions
            IS = ISINIT - MBOTWIDE  ! array
            IE = IEINIT - MBOTWIDE  ! array

C.............  Make sure that all of the needed real values are real...

C.............  Emissions and associated data
            DO V = 1, NVAR

C.................  Update start and end positions
                DO K = 1, NMBPPOL3
                    IS( K ) = IS( K ) + MBOTWIDE
                    IE( K ) = IE( K ) + MBOTWIDE
                END DO

                IF( .NOT. CHKREAL( LINE( IS(1):IE(1) ) ) .OR.
     &              .NOT. CHKREAL( LINE( IS(2):IE(2) ) )      ) THEN

                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'ERROR: Inventory data for "' //
     &                     TMPNAM( V )( 1:L ) // '" are not a number '//
     &                     'or have bad formatting at line', IREC
                    CALL M3MESG( MESG )

                END IF

                IF( NWARN .LT. MXWARN .AND.
     &              LINE( IS(1):IE(1) ) .EQ. ' ' .AND.
     &              LINE( IS(2):IE(2) ) .EQ. ' '       ) THEN

                    L = LEN_TRIM( TMPNAM( V ) )
                    WRITE( MESG,94010 ) 'WARNING: All emissions ' //
     &                     'data for ' // TMPNAM( V )( 1:L ) //  
     &                     ' are missing at line', IREC
                    CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                    LINE( IS(1):IE(1) ) = '0.'
                    LINE( IS(2):IE(2) ) = '0.'

                END IF

            END DO

C.............  If there has been an error, do not try to store any of the
C               records.  Instead  go to next line of file.
            IF( EFLAG ) CYCLE
       
C.............  Now use the file format definition to parse the LINE into
C               the various data fields...

            FIP  = ICC * 100000 + 1000 * STR2INT( LINE( 1:2 ) ) +
     &             STR2INT( LINE( 3:5 ) )

            TSCC = ADJUSTL ( LINE( 6:15 ) )

C.............  Extract vehicle type and speeds
            IF( OLDFMT ) THEN
                VTYPE    = ADJUSTL ( LINE( 29:33 ) )
                SPEEDCHR = ADJUSTL ( LINE( 34:35 ) )
            ELSE
                VTYPE    = ADJUSTL ( LINE( 16:20 ) )
                SPEEDCHR = ADJUSTL ( LINE( 21:22 ) )
            END IF

C.............    Determine if vehicle type is valid
	    J = FINDC( VTYPE, NVTYPE, CVTYPLST )

            IF ( J .LE. 0 ) THEN 
                EFLAG = .TRUE.
	        INVALID = .TRUE.
                WRITE( MESG, 94010 )
     &	              'ERROR: Vehicle type "' // VTYPE // 
     &                '" not found in list of valid types at line', IREC
                CALL M3MESG( MESG )
            ELSE
                IVT = IVTIDLST( J )

            END IF

C.............    Determine if road class is valid
            RWT = STR2INT( TSCC( 8:10 ) )
            J = FIND1( RWT, NRCLAS, AMSRDCLS )

            IF ( J .LE. 0 ) THEN
                EFLAG = .TRUE.
	        INVALID = .TRUE.
                WRITE( MESG, 94010 )
     &	              'ERROR: Road class "' // TSCC( 8:10 ) // 
     &                '" not found in list of valid types at line', IREC
                CALL M3MESG( MESG )
            ELSE
                RWT = RDWAYTYP( J )

            END IF

C.............  Determine if speed is valid - blank or zero will prompt
C               message that speeds from another file will be used

            IF ( NWARN .LT. MXWARN .AND. 
     &         ( SPEEDCHR .EQ. ' ' .OR.  
     &           SPEEDCHR .EQ. '0'      ) ) THEN
                SFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &	               'WARNING: Blank or zero speed found ' //
     &                 'on line ', IREC
                CALL M3MESG( MESG )
                NWARN = NWARN + 1
                ISPEED = 0
            ELSE

                ISPEED = STR2INT( SPEEDCHR )

            END IF

C.............  Determine if SCC is valid

            SCCLEN = LEN_TRIM( TSCC )

            IF ( SCCLEN .NE. 10 ) THEN   ! check if SCC is 10 char. wide
                EFLAG = .TRUE.
	        INVALID = .TRUE.
                WRITE( MESG, 94010 )
     &	                   'SCC not 10 characters wide on line ', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Initialize start and end positions
            IS = ISINIT - MBOTWIDE  ! array
            IE = IEINIT - MBOTWIDE  ! array

C.............  Sum emissions for invalid records
            IF ( INVALID ) THEN

                IDROP = IDROP + 1
                DO V = 1, NVAR
                    VDROP( V ) = VDROP( V ) + 
     &                           STR2REAL( LINE( IS(1):IE(1) ) )
                END DO
                CYCLE

            END IF
	
C.............  Create string source characteristics. Pad with zeros,
C               if needed
C.............  Make adjustments to pad with zeros, if needed
            WRITE( CFIP,94120 ) FIP
            CALL PADZERO( CFIP )
            CALL PADZERO( TSCC )
            WRITE( CRWT,RWTFMT ) RWT
            WRITE( CIVT,VIDFMT ) IVT

C.............  Store source characteristics if dimension is okay
            SS = SS + 1

            IF( SS .LE. NRAWIN ) THEN
                IFIPA  ( SS ) = FIP
                IRCLASA( SS ) = RWT
                IVTYPEA( SS ) = IVT
                CLINKA ( SS ) = CLNK
                CSCCA  ( SS ) = TSCC
                CVTYPEA( SS ) = VTYPE
                SPEEDA ( SS ) = REAL( ISPEED )
                TPFLGA ( SS ) = TPF
                INVYRA ( SS ) = INY
                XLOC1A ( SS ) = BADVAL3
                YLOC1A ( SS ) = BADVAL3
                XLOC2A ( SS ) = BADVAL3
                YLOC2A ( SS ) = BADVAL3
            END IF 

                
            DO V = 1, NVAR

C.................  Update start and end positions
                DO K = 1, NMBPPOL3
                    IS( K ) = IS( K ) + MBOTWIDE
                    IE( K ) = IE( K ) + MBOTWIDE
                END DO

                DANN = STR2REAL( LINE( IS(1):IE(1) ) )
                IF( .NOT. OLDFMT ) 
     &              DOZN = STR2REAL( LINE( IS(2):IE(2) ) )
                 
                ES = ES + 1

                IF ( ES .LE. NRAWBV ) THEN

                    I = POLPOS( V )
                    POLVLA ( ES,NEM ) = INVDCNV( I,1 ) * DANN
                    POLVLA ( ES,NOZ ) = INVDCNV( I,2 ) * DOZN

                    WRITE( CCOD,94125 ) I

                    CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC, ' ', 
     &                            ' ', CCOD, CSOURCA( ES ) )

                END IF

            END DO     ! loop through data variables
	
        END DO             ! end of loop for reading input file
              
12      CONTINUE           ! end of read on input file

C.........  Update saved cumulative counts
        NSRCSAV = SS       !  source
        NSRCVAR = ES       !  source*pollutant

	NDROP   = NDROP + IDROP

C.........  Write message if overflow occurred
        IF( NSRCSAV .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source memory allocation ' //
     &             'insufficient for IDA inventory'
            CALL M3MSG2( MESG )

        END IF

        IF( NSRCVAR .GT. NRAWBV ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of valid src x variables ' //
     &        'encountered: ', NRAWOUT, CRLF() // BLANK5 //
     &        'Maximum number of raw records allowed: ', NRAWBV

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        ELSE
            NRAWOUT = NSRCVAR

        END IF
	
C.........  Deallocate local allocatable arrays 
        DEALLOCATE( TMPNAM, POLPOS )

C.........  Return from subroutine 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDIDAMB
