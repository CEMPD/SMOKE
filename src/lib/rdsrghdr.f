
        SUBROUTINE RDSRGHDR( FDEV, SRGFMT, GDNAM, GDESC, XCENT, YCENT, 
     &                       XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the header of the spatial surrogates file.
C     
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures
        INCLUDE 'FLTERR.EXT'    !  error filter statement function


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2    CRLF
        LOGICAL        DSCM3GRD
        INTEGER        FIND1
        INTEGER        STR2INT
        REAL           STR2REAL
        LOGICAL        CHKINT
	LOGICAL        CHKREAL

        EXTERNAL       CRLF, DSCM3GRD, FIND1, STR2INT, STR2REAL, 
     &                 CHKINT, CHKREAL

C...........   Subroutine arguments

        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        CHARACTER(*) , INTENT  (OUT):: SRGFMT     ! Format of surrogates file
        CHARACTER(*) , INTENT  (OUT):: GDNAM      ! Grid name
        CHARACTER(*) , INTENT  (OUT):: GDESC      ! Grid description
        REAL         , INTENT  (OUT):: XCENT      ! Center of coordinate system
        REAL         , INTENT  (OUT):: YCENT      ! Center of coordinate system
        REAL         , INTENT  (OUT):: XORIG      ! X origin
        REAL         , INTENT  (OUT):: YORIG      ! Y origin
        REAL         , INTENT  (OUT):: XCELL      ! Cell size, X direction
        REAL         , INTENT  (OUT):: YCELL      ! Cell size, Y direction
        INTEGER      , INTENT  (OUT):: NCOLS      ! # cells in X direction
        INTEGER      , INTENT  (OUT):: NROWS      ! # cells in Y direction

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 16          ! # of potential line segments
        INTEGER, PARAMETER :: MXGRDTYP = 5

C...........   Grid types and names arrays
        INTEGER      :: GRDTYPES( MXGRDTYP ) = ( / LATGRD3
     &                                           , LAMGRD3
     &                                           , MERGRD3
     &                                           , STEGRD3
     &                                           , UTMGRD3 / )

        CHARACTER*15 :: GRDNAMES( MXGRDTYP ) = ( / 'LAT-LON        '
     &                                           , 'LAMBERT        '
     &                                           , 'MERCATOR       '
     &                                           , 'STEREOGRAPHIC  '
     &                                           , 'UTM            ' / )

C...........   Other arrays

        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables

        INTEGER       I, J       ! counters and indices
        INTEGER       IOS        ! i/o status
        INTEGER       IREC       ! record counter
        INTEGER       NTHIK      ! boundary thickness

        REAL          P_ALP      ! First map projection description parameter
        REAL          P_BET      ! Second map projection description parameter
        REAL          P_GAM      ! Third map projection description parameter

        LOGICAL    :: EFLAG = .FALSE.           ! Error flag

        CHARACTER*7            PROJUNIT   ! Projection units
        CHARACTER*80           MSGEND     ! tmp end of message
        CHARACTER*200          LINE       ! Read buffer for a line
        CHARACTER*300          MESG       ! Message buffer
        CHARACTER(LEN=IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(LEN=IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(LEN=IOVLEN3) PROJTYPE   ! projection type

        CHARACTER*16 :: PROGNAME = 'RDSRGHDR'     ! Program name

C***********************************************************************
C   Begin body of subroutine RDSRGHDR

C.........  Read the Models-3 grid information file
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IREC  = 0
        EFLAG = .FALSE.

        REWIND( FDEV )

C.........  Loop through lines of file until the header line is encountered
        DO

            READ ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &            'I/O error', IOS, 'reading gridding surrogates '//
     &            'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL UPCASE( LINE )

            IF ( LINE .EQ. ' ' ) CYCLE  ! skip all blank lines
  
            I = INDEX( LINE, '#GRID' ) ! determine if current line is the header
 
            IF ( I .GT. 0 ) THEN
                SRGFMT = 'MODELS3'  ! set format to 'MODELS3'
                IF( IREC .NE. 1 ) THEN   
                    MESG = 'Line 1 of the spatial surrogates file '//
     &                     'did not contain the header line.'
                    CALL M3WARN( PROGNAME, 0, 0, MESG )
                END IF

                EXIT                ! exit read loop with header found

            END IF 

        END DO

        SELECT CASE( SRGFMT )

        CASE ( 'MODELS3' )

C.............  Parse the line of data into segments based on the rules
C               for "list-formatted" in fortran, but not requiring 
C               quotes around the text strings

            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Make sure appropriate segments of header line appear to
C               be either INTEGER or REAL as expected 

            IF ( .NOT. CHKREAL( SEGMENT( 3 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 4 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 5 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 6 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 7 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 8 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 9 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 12 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 13 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 14 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 15 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 16 ) ) ) EFLAG = .TRUE.

	    IF ( EFLAG ) THEN

                WRITE( MESG, 94010 )
     &            'Unexpected data type encountered in header of '//
     &            'spatial surrogates file'
                CALL M3MESG( MESG )

            ELSE

                GDNAM    =           SEGMENT ( 2 )
                XORIG    = STR2REAL( SEGMENT ( 3 ) )
                YORIG    = STR2REAL( SEGMENT ( 4 ) )
                XCELL    = STR2REAL( SEGMENT ( 5 ) )
                YCELL    = STR2REAL( SEGMENT ( 6 ) )
                NCOLS    = STR2INT ( SEGMENT ( 7 ) )
                NROWS    = STR2INT ( SEGMENT ( 8 ) )
                NTHIK    = STR2INT ( SEGMENT ( 9 ) )
                PROJTYPE =           SEGMENT ( 10 )
                PROJUNIT =           SEGMENT ( 11 )
                P_ALP    = STR2REAL( SEGMENT ( 12 ) )
                P_BET    = STR2REAL( SEGMENT ( 13 ) )
                P_GAM    = STR2REAL( SEGMENT ( 14 ) )
                XCENT    = STR2REAL( SEGMENT ( 15 ) )
                YCENT    = STR2REAL( SEGMENT ( 16 ) )

            END IF

        END SELECT

C.........  Now check header of surrogates file to ensure that this file
C           is consistent with the grid definition

        MSGEND = 'in surrogates does not match the grid description'

        IF( GDNAM .NE. GDNAM3D ) THEN
            EFLAG = .TRUE.
            MESG = 'Grid name "' // GDNAM( 1:LEN_TRIM( GDNAM ) ) //
     &             '" ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( XORIG, SNGL( XORIG3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'X-origin ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( YORIG, SNGL( YORIG3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Y-origin ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( XCELL, SNGL( XCELL3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'X cell size ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( YCELL, SNGL( YCELL3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Y cell size ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( REAL( NCOLS ), REAL( NCOLS3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Number of columns ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( REAL( NROWS ), REAL( NROWS3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Number of rows ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( REAL( NTHIK ), REAL( NTHIK3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Number of boundary cells ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        J =  FIND1( GDTYP3D, MXGRDTYP, GRDTYPES )

        IF( J .LE. 0 ) THEN

            WRITE( MESG,94010 ) 'INTERNAL ERROR: Could not find ' //
     &             'grid type', GDTYP3D, 'in subroutine grid type array'
     &             
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE

            IF( PROJTYPE .NE. GRDNAMES( J ) ) THEN
                EFLAG = .TRUE.
                MESG = 'Projection type "' // 
     &                 PROJTYPE( 1:LEN_TRIM( PROJTYPE ) ) //
     &                 '" ' // MSGEND
                CALL M3MSG2( MESG )
            END IF

        END IF

        IF( FLTERR( P_ALP, SNGL( P_ALP3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = '1st map projection parameter ' //
     &             CRLF() // BLANK10 // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( P_BET, SNGL( P_BET3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = '2nd map projection parameter ' //
     &             CRLF() // BLANK10 // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( P_GAM, SNGL( P_GAM3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = '3rd map projection parameter ' //
     &             CRLF() // BLANK10 // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( XCENT, SNGL( XCENT3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Projection X center ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( FLTERR( YCENT, SNGL( YCENT3D ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Projection Y center ' // MSGEND
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG ) THEN
            MESG = 'Problem with surrogates header.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C.........  If get here, then the header was never found!
999     MESG = 'Surrogates file is missing a header line'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSRGHDR

        



