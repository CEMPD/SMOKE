
        SUBROUTINE RDEMSMB( EDEV, INY, NRAWIN, NRAWBP, WKSET, 
     &                      NRAWOUT, IOS, IREC, ERFILDSC, EFLAG, 
     &                      NDROP, DDROP )

C***********************************************************************
C  subroutine body starts at line 176
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 format for the mobile source formatted
C      files. It can be called multiple times for multiple files.
C
C  PRECONDITIONS REQUIRED:
C      Files must be opened and their unit numbers stored in EDEV() in the
C      order listed in the description.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDCSRC, CHECKMEM 
C      Functions: I/O API functions, YR2DAY
C
C  REVISION  HISTORY:
C      Copied from rdemsar.f by M. Houyoux (2/2000)
C
C****************************************************************************
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

C.........  MODULES for public variables
C.........  This module is the inventory arrays
        USE MODSOURC

C.........  This module is for mobile-specific data
        USE MODMOBIL

C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        LOGICAL                ENVYN
        INTEGER                FIND1
        INTEGER                FINDC
        INTEGER                GETVMIX
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY  

        EXTERNAL   CRLF, ENVYN, FIND1, FINDC, GETVMIX, INDEX1, STR2INT, 
     &             STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: EDEV      ! unit no. for input file
        INTEGER     , INTENT (IN) :: INY       ! inv year for this set of files
        INTEGER     , INTENT (IN) :: NRAWIN    ! total raw record-count
        INTEGER     , INTENT (IN) :: NRAWBP    ! total raw record times pols
        INTEGER     , INTENT (IN) :: WKSET     ! weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT   ! valid raw record-count
        INTEGER     , INTENT(OUT) :: IOS       ! I/O status
        INTEGER     , INTENT(OUT) :: IREC      ! line number
        CHARACTER(*), INTENT(OUT) :: ERFILDSC  ! file desc of file in error
        LOGICAL     , INTENT(OUT) :: EFLAG     ! error flag 
        INTEGER   , INTENT(INOUT) :: NDROP     ! number of records dropped
        REAL      , INTENT(INOUT) :: DDROP( MXIDAT ) ! emis dropped per pol/act

C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables
        INTEGER, PARAMETER :: CRVLEN3 = VTPLEN3+RWTLEN3

C.........  Local allocatable arrays
        CHARACTER*20, ALLOCATABLE :: SEGMENT( : )  ! Segments of parsed lines

C...........   Local arrays
        REAL             INVAL  ( MXDATFIL ) !  temporary data values
        REAL             DATAVAL( MXDATFIL ) !  temporary data values

C...........   Other local variables

        INTEGER          I, J, K, K1, K2, L, V       ! counters and indices

        INTEGER          COD         !  temporary pollutant code number
        INTEGER          CYID        !  tmp county ID
        INTEGER          ES          !  counter for rec x data vars
        INTEGER          FIP         !  temporary fip
        INTEGER          FMTCASE     !  code for format case
        INTEGER       :: ICC  = 0    !  country code, def = 0
        INTEGER          IVT         !  tmp vehicle type number
        INTEGER          MXTCOL      !  maximum number of table columns
        INTEGER       :: NPOA = 0    !  number of input data variables
        INTEGER          NPRECOL     !  no. src char columns for list-directed
        INTEGER, SAVE :: NSRCDAT = 0 !  cumulative source x pollutants count
        INTEGER, SAVE :: NSRCSAV = 0 !  cumulative source count
        INTEGER          ROAD        !  temporary road class code
        INTEGER          RWT         !  temporary roadway type
        INTEGER          SS          !  counter for rec (no data vars)
        INTEGER          STID        !  tmp state ID
        INTEGER          TPF         !  temporary temporal ID
        INTEGER          ZONE        !  tmp time zone

        REAL             DAY2YR    !  local, leap-year-able, DAY to YEAR factor
        REAL             X1, X2    !  x-dir link end point coordinates
        REAL             XLOC      !  tmp x-coordinate
        REAL             Y1, Y2    !  y-dir link end point coordinates 
        REAL             YLOC      !  tmp y-coordinate

        LOGICAL, SAVE :: FFLAG    = .TRUE.  ! true: fixed-column file format
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL       :: LFLAG    = .FALSE. ! true: link file
        LOGICAL, SAVE :: SFLAG    = .TRUE.  ! true: speeds data are available
        LOGICAL, SAVE :: WFLAG    = .TRUE.  ! true: convert coords to west-hemi

        CHARACTER*1            ARTP     ! tmp area type (urban/rural)
        CHARACTER*4            FCTP     ! tmp facility type
        CHARACTER*10 , SAVE :: FIPFMT   ! formt to write co/st/cy to string
        CHARACTER*10 , SAVE :: RWTFMT   ! formt to write rdway type to string
        CHARACTER*10 , SAVE :: VIDFMT   ! formt to write vehicle type to string
        CHARACTER*10           HDRBUF   ! tmp header buffer
        CHARACTER*20           EVNAME   ! tmp environment variable name
        CHARACTER*40           FMTDESC  ! tmp format description
        CHARACTER*300          MESG     ! message buffer
        CHARACTER*600          LINE     ! input line from inventory file

        CHARACTER(LEN=IOVLEN3) CDAT     ! tmp data variable code
        CHARACTER(LEN=POLLEN3) CCOD     ! character data var index to INVDNAM
        CHARACTER(LEN=FIPLEN3) CFIP     ! character FIP code
        CHARACTER(LEN=RWTLEN3) CROAD    ! tmp road class no.
        CHARACTER(LEN=RWTLEN3) CRWT     ! tmp roadway type no.
        CHARACTER(LEN=CRVLEN3) CRVC     ! tmp roadway // vehicle type
        CHARACTER(LEN=LNKLEN3) CLNK     ! tmp link ID
        CHARACTER(LEN=VIDLEN3) CIVT     ! tmp character vehicle ID
        CHARACTER(LEN=SCCLEN3) TSCC     ! tmp character SCC
        CHARACTER(LEN=VTPLEN3) VTYPE    ! tmp vehicle type

        CHARACTER*16 :: PROGNAME = 'RDEMSMB' ! Program name

C***********************************************************************
C   begin body of subroutine RDEMSMB

C.........  For first time routine called
        IF( FIRSTIME ) THEN

C.............  Get environment variable values for this routine
            EVNAME = 'SMK_EMS95_FIXFMT'
            MESG = 'Indicator for fixed-column EMS-95 format'
            FFLAG = ENVYN( EVNAME, MESG, FFLAG, IOS )

            EVNAME = 'WEST_HSPHERE'
            MESG = 'Western hemisphere flag'
            WFLAG = ENVYN( EVNAME, MESG, .TRUE., IOS )

C.............  Ensure that vehicle mix data are available
            IF( .NOT. ALLOCATED( VMTMIXA ) ) THEN

                MESG = 'Mobile VMT mix data are required for import ' //
     &                 'of EMS-95 mobile format'// CRLF() // BLANK10 //
     &                 'Set the IMPORT_VMTMIX_YN environment ' //
     &                 'variable to Y and try again.'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Set flag for whether speeds data are available or not - give a
C               warning message if speeds will not be used.
            IF( .NOT. ALLOCATED( SPDTBLA ) ) THEN
                SFLAG = .FALSE.
                MESG= 'WARNING: mobile speeds data are not available '//
     &                'during import of EMS-95 '// CRLF() // BLANK10 //
     &                'mobile format.  Only speeds from the MPREF ' //
     &                'file will be used in ' // CRLF() // BLANK10 //
     &                'subsequent processing.'
                CALL M3MSG2( MESG )
            END IF

C.............  Set up formats
            WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3
            WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3
            WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3

            FIRSTIME = .FALSE.

       END IF

C.........  Determine whether this file is a non-link or link data file
        DO

C.............  Read a line of file and check input status
            READ( EDEV, 93000, END=399, IOSTAT=IOS ) LINE

            IF( IOS .NE. 0 ) CYCLE
            
            IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN

                L = LEN_TRIM( LINE )
                LINE = ADJUSTL( LINE( 2:L ) )
                L = INDEX( LINE, ' ' )
                HDRBUF = LINE( 1:L )

                IF( HDRBUF .EQ. 'LINK' ) THEN
                    LFLAG = .TRUE.
                    EXIT
                ELSE IF( HDRBUF .EQ. 'NONLINK' ) THEN
                    LFLAG = .FALSE.
                    EXIT
                END IF

                IF( HDRBUF .EQ. 'UNITS' ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: #UNITS header is not ' //
     &                     'supported in EMS-95 mobile format.'
                    CALL M3MSG2( MESG )
                END IF

            END IF
             
        END DO

C.........  Set the description of the file type and format
        IF( .NOT. LFLAG .AND. FFLAG ) THEN
            FMTDESC = 'non-link EMS-95 fixed format'
            FMTCASE = 1

        ELSE IF( .NOT. LFLAG ) THEN
            FMTDESC = 'non-link EMS-95 list-directed format'
            FMTCASE = 2
            NPRECOL = 2

        ELSE IF( LFLAG .AND. FFLAG ) THEN
            FMTDESC = 'link EMS-95 fixed format'
            FMTCASE = 3

        ELSE IF( LFLAG ) THEN
            FMTDESC = 'link EMS-95 list-directed format'
            FMTCASE = 4
            NPRECOL = 8

        ELSE
            MESG = 'INTERNAL ERROR: mobile EMS-95 file format '//
     &             'is not recognized'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Make sure the file is at the beginning
        REWIND( EDEV )

C.........  Initialize units as converstion from day to year, in case units
C           are not provided by file
        INVDCNV = 1. / YR2DAY( INY )    ! Array

C.........  Initialize before loop
        SS   = NSRCSAV
        ES   = NSRCDAT
        IREC = 0
        ERFILDSC = 'EMS-95 mobile data'
        TPF      = WKSET     ! data expected as average day

C.........  Head of file read loop
        DO 

C.............  Read a line of file and check input status
            READ( EDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check read i/o status
            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC,
     &              CRLF() // BLANK10 // 'using ' // FMTDESC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            IF ( LINE .EQ. ' ' ) CYCLE      ! skip if line is blank

C.............  Scan for header lines and check to ensure all are set 
C               properly
            CALL GETHDR( MXDATFIL, .FALSE., .FALSE., .TRUE., 
     &                   LINE, ICC, INY, NPOA, IOS )

C.............  Interpret error status
            IF( IOS .EQ. 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXDATFIL=', MXDATFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Set the number of table columns and allocate memory
            IF( .NOT. ALLOCATED( SEGMENT ) ) THEN

                MXTCOL = NPRECOL + NPOA
                ALLOCATE( SEGMENT( MXTCOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

            END IF

C.............  Define day to year conversion factor and real type for integer 
C               missing value
            DAY2YR  = 1. / YR2DAY( INY )

C.............  Initialize link information in case record is nonlink         
            CLNK = ' '
            X1   = BADVAL3
            Y1   = BADVAL3
            X2   = BADVAL3
            Y2   = BADVAL3
            ZONE = 0

C.............  Set temporary variables for data from this file, depending
C               on the format of the file...
C.............  For non-link fixed format...
            SELECT CASE( FMTCASE )
            CASE( 1 )
                STID = STR2INT( LINE(  1:2  ) )
                CYID = STR2INT( LINE(  3:5  ) )
                ARTP = ADJUSTL( LINE(  6:6  ) )
                FCTP = ADJUSTL( LINE(  7:10 ) )
                INVAL( 1 ) = STR2REAL( LINE( 11:18 ) )
                ROAD = STR2INT( ARTP // FCTP ) 
                FIP = 1000*STID + CYID

C.................  Error if roadway type will be larger than 3 characters
                IF( ROAD .GT. 999 ) THEN
                    MESG = 'INTERNAL ERROR: internal roadway field ' //
     &                     'is not large enough to handle input file.'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
                END IF

C.............  For non-link list-directed format
            CASE( 2 )

                CALL PARSLINE( LINE, MXTCOL, SEGMENT )

                FIP   = STR2INT ( SEGMENT( 1 ) )
                ROAD  = STR2INT ( SEGMENT( 2 ) )

                DO I = 1, NPOA
                    K = NPRECOL + I
                    INVAL( I ) = STR2REAL( SEGMENT( K ) )
                END DO

C.............  For link fixed format...
            CASE( 3 )

                MESG = 'Link EMS-95 fixed format is not yet supported'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.............  For link list-directed format
            CASE( 4 )

                CALL PARSLINE( LINE, MXTCOL, SEGMENT )

                FIP   = STR2INT ( SEGMENT( 1 ) )
                ROAD  = STR2INT ( SEGMENT( 2 ) )
                CLNK  = ADJUSTL ( SEGMENT( 3 ) )
                X1    = STR2REAL( SEGMENT( 4 ) )
                Y1    = STR2REAL( SEGMENT( 5 ) )
                X2    = STR2REAL( SEGMENT( 6 ) )
                Y2    = STR2REAL( SEGMENT( 7 ) )
                ZONE  = STR2INT ( SEGMENT( 8 ) )

                DO I = 1, NPOA
                    K = NPRECOL + I
                    INVAL( I ) = STR2REAL( SEGMENT( K ) )
                END DO

            END SELECT 

C.............  Include country in country/state/county code
            IF( FIP .LT. 100000 ) THEN
                FIP = ICC * 100000 + FIP
            END IF

C.............  Create character string version of integer fields
            WRITE( CFIP ,FIPFMT ) FIP
            WRITE( CROAD,RWTFMT ) ROAD

C.............  Convert link coordinates from UTM to lat-lon
            IF( FMTCASE .GE. 3 .AND. ZONE .GT. 0 ) THEN
                XLOC = X1
                YLOC = Y1
                CALL UTM2LL( XLOC, YLOC, ZONE, X1, Y1 )

                XLOC = X2
                YLOC = Y2
                CALL UTM2LL( XLOC, YLOC, ZONE, X2, Y2 )

C.............  If needed, convert lat-lon coordinates to western hemisphere
            ELSE IF( FMTCASE .GE. 3 .AND. WFLAG ) THEN
                IF( X1 .GT. 0 ) X1 = -X1
                IF( X2 .GT. 0 ) X2 = -X2

            END IF

c note: issues - 
c    n: match up with SCC
c    n: match up with speed data file

C.............  Match up with unsorted vehicle mix table
            K1 = GETVMIX( CFIP, CROAD, CLNK )

C.............  Set error flag if return status is zero or less
C.............  Message are written out in GETVMIX
            IF( K1 .LE. 0 ) THEN
                EFLAG = .TRUE.
		NDROP = NDROP + 1
                CYCLE
            END IF

C.................  Convert road type number, and set character value also
            K2 = FIND1( ROAD, NRCLAS, AMSRDCLS )

            IF ( K2 .LE. 0 ) THEN   ! determine if road type is blank
                 EFLAG = .TRUE.
                 WRITE( MESG, 94010 )
     &                'ERROR: Road class "' // CROAD //
     &                '" not found in list of valid types at line', IREC
                 CALL M3MESG( MESG )

            ELSE
                RWT = RDWAYTYP( K2 )
                WRITE( CRWT,RWTFMT ) RWT

            END IF

C.............  Loop through vehicle types
            DO J = 1, NVTYPE

C.................  Set vehicle type code and name
                IVT   = IVTIDLST( J )
                VTYPE = CVTYPLST( J )
                WRITE( CIVT,VIDFMT ) IVT

C.................  Search for SCC
                CRVC = CROAD // VTYPE
                K2 = FINDC( CRVC, NSCCTBL, SCCRVC )

C.................  Assign SCC
                IF( K2 .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not find SCC for ' //
     &                     CRLF() // BLANK10 //
     &                     'Road Type: '// CROAD// ' Vtype: '// VTYPE
                    CYCLE

                ELSE
                    TSCC = SCCTBL( K2 )

                END IF

C.................  Match up with speed tables if they exist
c                SPD = SPDTBL( K2 )
c note: add later - SPEED is to be treated as an activity like VMT

C.................  Compute vehicle-mix-adjusted data values
                DATAVAL( 1:NPOA ) = VMTMIXA( J,K1 ) * INVAL( 1:NPOA )
C.................  Time to store data in unsorted list
                SS = SS + 1

                IF ( SS .LE. NRAWIN ) THEN

                    IFIPA  ( SS ) = FIP
                    IRCLASA( SS ) = RWT
                    IVTYPEA( SS ) = IVT
                    CLINKA ( SS ) = CLNK
                    CVTYPEA( SS ) = VTYPE
                    TPFLGA ( SS ) = TPF
                    INVYRA ( SS ) = INY
                    CSCCA  ( SS ) = TSCC
                    XLOC1A ( SS ) = X1
                    YLOC1A ( SS ) = Y1
                    XLOC2A ( SS ) = X2
                    YLOC2A ( SS ) = Y2
                END IF

C.................  Store data variable values
                DO V = 1, NPOA

                    ES = ES + 1
                    COD  = DATPOS( V )

                    IF( ES .LE. NRAWBP ) THEN
                        POLVLA ( ES,1 ) = INVDCNV( COD ) * DATAVAL( V )

                        WRITE( CCOD,94125 ) COD

                        CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC, ' ', 
     &                                ' ', CCOD, CSOURCA( ES ) )
                    END IF

                END DO
            END DO          !  vehicle types
        END DO              !  to head of file read loop

199     CONTINUE        !  end of the file read loop

        CLOSE( EDEV )

        WRITE( MESG,94010 ) 
     &         'DATA FILE processed:'  // CRLF() // BLANK10 //
     &              'This-file source-count', SS - NSRCSAV,
     &         CRLF() // BLANK10 //
     &              'Cumulative source-count', SS,
     &         CRLF() // BLANK10 //
     &              'This-file source*data-count', ES - NSRCDAT,
     &         CRLF() // BLANK10 //
     &              'Cumulative source*data-count', ES

        CALL M3MSG2( MESG )

        NSRCDAT = ES        !  cumulative records * data variables
        NSRCSAV = SS        !  cumulative records 

        IF( NSRCDAT .GT. NRAWBP ) THEN

            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Source time data memory ' //
     &             'allocation insufficient for inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCDAT

        END IF		!  if overflow or if errors

C.........  Deallocate local memory
        DEALLOCATE( SEGMENT )

        RETURN

399     MESG = 'File does not contain "#LINK" or "#NONLINK" ' //
     &         'header line with the # sign in ' // CRLF() // BLANK10//
     &         'the first column. Edit file and try again.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDEMSMB
